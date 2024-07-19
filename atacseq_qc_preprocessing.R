library(tidyverse)
library(rjson)
library(parallel)
library(Rsubread)
# 1. ATAC-seq QC: Gather QC metrics from ENCODE ATAC output -----------
qc_json_dir <- "atac/qc_report/qc_json/"

qc_dir <- data.frame(dir = read.table("atac/qc_json_paths.txt")) %>%
  rename("dir" = "V1") %>%
  mutate(dir = paste0("atac/", dir)) %>%
  mutate(Folder = str_extract(dir, "(?<=qc_json/).*(?=/call-qc_report)"))

qc_metrics <- lapply(qc_dir$dir, function(dir) {
  unlist(fromJSON(file = dir))}) %>%
  do.call(cbind, .) %>%
  as.data.frame()%>%
  `colnames<-`(.["general.description", ]) %>%
  select(str_sort(colnames(.), numeric = TRUE)) %>%
  tibble::rownames_to_column("QC_metrics")

write.table(qc_metrics, "QC_metrics/QC_metrics_all.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)

qc_metrics_for_pre <- c(
  "align.frac_mito.rep1.non_mito_reads",
  "align.frag_len_stat.rep1.frac_reads_in_nfr_qc_pass",
  "align.frag_len_stat.rep1.nfr_peak_exists",
  "align.frag_len_stat.rep1.mono_nuc_peak_exists",
  "lib_complexity.lib_complexity.rep1.NRF",
  "lib_complexity.lib_complexity.rep1.PBC1",
  "lib_complexity.lib_complexity.rep1.PBC2",
  "replication.num_peaks.rep1.num_peaks",
  "align_enrich.tss_enrich.rep1.tss_enrich",
  "peak_enrich.frac_reads_in_peaks.macs2.rep1.frip",
  "align.samstat.rep1.pct_mapped_reads",
  "align.frac_mito.rep1.frac_mito_reads",
  "peak_stat.peak_region_size.rep1.25_pct",
  "peak_stat.peak_region_size.rep1.50_pct",
  "peak_stat.peak_region_size.rep1.75_pct")

qc_sub <- qc_metrics %>%
  filter(QC_metrics %in% qc_metrics_for_pre) %>%
  mutate(QC_metrics = factor(QC_metrics, qc_metrics_for_pre)) %>%
  arrange(QC_metrics) %>%
  mutate(QC_metrics = as.character(QC_metrics)) %>%
  mutate(QC_metrics = str_remove(QC_metrics, ".*\\."))

write.csv(qc_sub, "ATAC_QC.csv", row.names = FALSE)

# 2. ATAC-seq Preprocessing ------------
## 2.1. After running ENCODE ATAC pipeline, Gather .narrowPeak files from “call_peak” outputs --------
system("mkdir narrowPeaks")
system("cp samples/**/call-call_peak/**/**/*.bfilt.narrowPeak.gz narrowPeaks")
system("gunzip ./narrowPeaks/*.narrowPeak.gz")

## 2.2 Sort descending by q-values and take top 100K peaks for each sample ------
system("mkdir narrowPeaks_filtered")

narrowPeaks_dir <- dir(path = "narrowPeaks", pattern="*.narrowPeak", full.names = TRUE)
sample_names <- str_extract(narrowPeaks_dir, "(A|B)[0-9]+")

narrowPeaks <- mclapply(narrowPeaks_dir, import, mc.cores = ncore)
names(narrowPeaks) <- sample_names

narrowPeaks_sort <- lapply(narrowPeaks, function(gr) {
  sort(gr, by = ~ qValue, decreasing = TRUE)})

narrowPeaks_sort_top100K <- lapply(narrowPeaks_sort, function(gr) {
  gr[1:1E5, ]})

lapply(names(narrowPeaks_sort_top100K), function(x) {
  export.bed(sort(narrowPeaks_sort_top100K[[x]]), 
             con = paste0("narrowPeaks_filtered/", x, 
                          ".narrowPeak.top_100k_peaks_by_qval.bed"))
})

## 2.3 Merge filtered peaks across all samples -------
system("mkdir merged_filtered_peaks")
system("use BEDTools")
system(paste0("cat narrowPeaks_filtered/*narrowPeak.top_100k_peaks_by_qval.bed | sort -k1,1 -k2,2n | ",
              "bedtools merge -i - > merged_filtered_peaks/merged_filtered_peaks.bed"))


## 2.4 Run FeatureCounts on merged peaks -----------
annot <- read.table("merged_filtered_peaks/merged_filtered_peaks.bed",
                    sep = "\t", header = FALSE) %>%
  `colnames<-`(c("Chr", "Start", "End")) %>%
  mutate(GeneID = paste0("peak_", format(row_number() - 1, trim = TRUE, 
                                         scientific = FALSE)), Strand = "-") %>%
  select(GeneID, everything())

bamfiles <- dir(path = "BAMs", pattern="*.bam$", full.names = TRUE)
ncore <- min(length(bamfiles), detectCores())

cnts <- featureCounts(
  files = bamfiles,
  isGTFAnnotationFile = FALSE,
  isPairedEnd = TRUE,
  annot.ext = annot, 
  nthreads = ncore
)
save(cnts, file = "atac_TCF1cKO_counts.RData")

cnts_samplenames <- str_extract(colnames(cnts$counts), "^.*(?=_S[0-9])") 
cnts_stat_samplenames <- str_extract(colnames(cnts$stat)[-1], "^.*(?=_S[0-9])")   

cnts_df <- data.frame(cnts$annotation, cnts$counts, stringsAsFactors = FALSE) %>%
  `colnames<-`(c(colnames(cnts$annotation), cnts_samplenames))%>%
  select(colnames(cnts$annotation), str_sort(cnts_samplenames, numeric = TRUE))

cnts_df_250 <- subset(cnts_df, Length > 250)

cnts_stat <- cnts$stat %>%
  `colnames<-`(c("Status", cnts_stat_samplenames))%>%
  select(Status, str_sort(cnts_stat_samplenames, numeric = TRUE))


write.table(cnts_df, file = "atac_TCF1cKO_counts.txt", 
            row.names = FALSE, sep = "\t")
write.table(cnts_df_250, file = "atac_TCF1cKO_counts_min_peak_length_250bp.txt",
            row.names = FALSE, sep = "\t")
write.table(cnts_stat, "atac_TCF1cKO_featureCounts_stat.txt", sep = "\t",
            row.names = FALSE)
