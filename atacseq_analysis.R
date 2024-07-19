library(tidyverse)
library(DESeq2)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(ComplexHeatmap)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# 1. Load Data ----------------------
counts <- read.table("data/ATACseq/atac_TCF1cKO_counts_min_peak_length_250bp.txt",
                     sep = "\t", header = TRUE)
metadata <- read.csv("data/ATACseq/atac_metadata.csv") 

counts_tbl <- counts %>% dplyr::select(contains("dLN")) %>% `rownames<-`(counts$GeneID)
metadata_tbl <- metadata %>% column_to_rownames("Sample_ID")

all(metadata$Sample_ID == colnames(counts_tbl))

# 2. Peak Annotation using ChIPseeker::annotatePeak ---------------
counts_GR <- makeGRangesFromDataFrame(
  counts[, 1:6],
  keep.extra.columns = TRUE,
  ignore.strand = FALSE,
  seqnames.field = "Chr",
  start.field = "Start",
  end.field = "End",
  strand.field = "Strand"
)

peakAnno_all <- annotatePeak(counts_GR, tssRegion = c(-3000, 3000), 
                             TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, 
                             overlap = "all", annoDb = "org.Mm.eg.db")

peakAnno_all_sub <- peakAnno_all@anno %>%
  as.data.frame() %>%
  dplyr::select(GeneID, seqnames, start, end, SYMBOL, annotation, distanceToTSS) %>%
  mutate(anno_short = ifelse(
    str_detect(annotation, "Intron|Exon"), 
    str_extract(annotation, "Intron|Exon"), annotation))

save(peakAnno_all, file = "R_objects/PeakAnno.RData")

# 3. Run DESeq2 ----------------------
atacDDS <- DESeqDataSetFromMatrix(counts_tbl, metadata_tbl, design = ~ genotype)
atacDDS <- DESeq(atacDDS)

atacRES <- results(atacDDS, contrast = c("genotype", "KO", "WT"))
atacRES_ordered_all <- atacRES %>%
  as.data.frame() %>%
  cbind(GeneID = rownames(.), .) %>%
  arrange(padj) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  mutate(direction = factor(direction, c("up", "down"))) %>%
  left_join(peakAnno_all_sub, by = "GeneID")

save(atacDDS, atacRES, atacRES_ordered_all,
     file = "R_objects/Tcf7_Th17_ATACseq_DA_res.RData")

counts_tbl_norm <- DESeq2::counts(atacDDS, normalized = TRUE)

# 3. Prepare Figures ----------------------
## Function for generating heatmaps
htmap <- function(RES, counts_tbl, width = 6, height = 9, highlight_genes = "", fig_num = "",
                  top = NULL, suffix = "", use_row_annot = FALSE, gap = c(2.5, 1.5),
                  most_DA_only = TRUE, promoter_only = FALSE) {
  
  DA_peaks <- RES %>%
    filter(padj < 0.05) %>%
    dplyr::select(GeneID, direction, padj, SYMBOL, annotation)
  
  htmap_df <- as.data.frame(counts_tbl) %>%
    filter(rowSums(.) != 0) %>%
    rownames_to_column("GeneID") %>%
    inner_join(DA_peaks, by = "GeneID") %>%
    mutate(direction = factor(direction, c("up", "down")))
  if(promoter_only) { htmap_df <- htmap_df %>% 
    filter(str_detect(annotation, "Promoter")) 
  }
  
  if(most_DA_only) {
    htmap_df <- htmap_df %>% 
      group_by(SYMBOL) %>%
      arrange(padj) %>% dplyr::slice(1) %>% ungroup()
    row_split <- htmap_df$direction
    cluster_row_slices <- FALSE
  }
  else {
    row_split = htmap_df$SYMBOL
    cluster_row_slices <- TRUE
   }
  
  htmap_df_scaled <- htmap_df %>%
    dplyr::select(contains("dLN"))%>% 
    t() %>% scale() %>% t()
  
  col_labels <- str_extract(colnames(htmap_df_scaled), "[0-9]")
  col_split <- str_extract(colnames(htmap_df_scaled), "WT|KO") %>%
    factor(levels = c("WT", "KO"))
   
  htmap_col <- circlize::colorRamp2(
    c(min(htmap_df_scaled, na.rm=T), mean(htmap_df_scaled, na.rm=T),
      max(htmap_df_scaled, na.rm=T)),
    c("navy", "white", "firebrick3"), transparency = 0.05)
  
  if(use_row_annot) {
    row_annot <- rowAnnotation(
      link = anno_mark(at = which(htmap_df$SYMBOL %in% highlight_genes),
                       labels = htmap_df$SYMBOL[htmap_df$SYMBOL %in% highlight_genes],
                       labels_gp = gpar(fontsize = 9, fontface = "italic"), 
                       padding = unit(1, "mm")))
    show_row_names <- FALSE; row_labels <- NULL
    }
  else {
    row_annot <- NULL
    show_row_names <- TRUE
    row_labels <- htmap_df$SYMBOL
    row_labels[!row_labels %in% highlight_genes] <- ""
  }
  
  col_annot <- columnAnnotation(
    genotype = col_split,
    col = list(genotype = setNames(c("black", "maroon"), c("WT", "KO"))),
    simple_anno_size = unit(4, "mm"),
    annotation_name_gp = gpar(fontsize = 0))
  
  p <- htmap_df_scaled %>%
    Heatmap(
      col = htmap_col,
      row_labels = row_labels,
      row_split = row_split,
      row_names_gp = gpar(fontsize = 9, fontface = "italic"),
      row_title_gp = gpar(fontsize = 0),
      show_row_names = show_row_names,
      column_labels = col_labels,
      column_title_gp = gpar(fontsize = 10, fontface = "bold"),
      column_names_gp = gpar(fontsize = 9),
      column_names_centered = TRUE,
      column_names_rot = 0,
      column_split = col_split,
      cluster_row_slices = FALSE,
      cluster_columns = TRUE,
      cluster_column_slices = FALSE,
      right_annotation = row_annot,
      top_annotation = col_annot,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      name = "Z-score",
      row_gap = unit(gap[1], "mm"),
      column_gap = unit(gap[2], "mm")) 
  
  png(paste0("figures/Fig.", fig_num, "_DA", suffix, ".png"), 
      width = width, height = height, res = 400, units = "in")
  draw(p)
  dev.off()
  pdf(paste0("figures/Fig.", fig_num, "_DA", suffix, ".pdf"), 
      width = width, height = height)
  draw(p)
  dev.off()
}

## Function for signature score computation, testing and visualization (boxplot)
## Signature score = average of DESeq2 normalized counts
sig_score_boxplot <- function(counts_tbl, meta, anno, sig, min_FDR = FALSE, fig_num = "", 
                              title = "", suffix = "", nudge_x = 0) {
  counts_tbl <- as.data.frame(counts_tbl)
  cnt_sub <- anno %>%
    as.data.frame() %>%
    dplyr::select(GeneID, SYMBOL) %>%
    left_join(counts_tbl %>% rownames_to_column("GeneID"), by = "GeneID") %>%
    filter(str_to_upper(SYMBOL) %in% str_to_upper(sig)) 
  
  check_genes <- setdiff(str_to_upper(sig), str_to_upper(cnt_sub$SYMBOL))
  
  if(min_FDR) {
    cnt_sub <- cnt_sub %>%
      left_join(atacRES_ordered %>% dplyr::select(GeneID, padj), by = "GeneID") %>%
      group_by(SYMBOL) %>%
      arrange(padj) %>%
      dplyr::slice(1) %>% dplyr::select(-padj)
    
    suffix <- paste0(suffix, "_min_FDR")
  }
  
  score_df <- cnt_sub %>%
    pivot_longer(3:ncol(.), names_to = "Sample_ID", values_to = "expr")  %>%
    group_by(Sample_ID) %>%
    summarise(sig_score = mean(expr)) %>%
    left_join(meta, by = "Sample_ID") %>%
    mutate(genotype = factor(genotype, c("WT", "KO")))
  
  sig_score_WT <- subset(score_df, genotype == "WT")$sig_score
  sig_score_KO <- subset(score_df, genotype == "KO")$sig_score
  
  t_test_pval <- t.test(sig_score_WT, sig_score_KO)$p.value
  
  score_df %>%
    ggplot(aes(x = genotype, y = log2(sig_score), color = genotype, shape = genotype)) +
    geom_boxplot(width = 0.3, outlier.alpha = 0, size = 0.45, fill = NA) +
    geom_point(size = 1.5, position = position_nudge(x = nudge_x, y = 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_color_manual(values = c("black", "maroon")) +
    labs(x = NULL, y = "Log2 Signature Score", title = title) +
    guides(color = guide_legend(keywidth = 0.8, keyheight = 0.8)) +
    theme_cowplot(font_size = 10) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank())
  
  file_path <- paste0("figures/Fig.", fig_num, "_sig_score_ATACseq_", 
                      str_replace_all(title, " ", "_"), suffix)
  
  ggsave(paste0(file_path, ".png"), width = 4, height = 3.2)
  ggsave(paste0(file_path, ".pdf"), width = 4, height = 3.2)
  
  return(list(check_genes = check_genes, t_test_pval = round(t_test_pval, 5)))
}


## Fig 6F: Heatmap of npTh17 ATACseq DA peaks ---------
gene_list_Fig.6F <- c("Rorc", "Il17a", "Il1r1", "Csf2", "Il23r", "Icos", 
                      "Gzmb", "Ccl3", "Il7r", "Lag3", "Cxcr6", 
                      "Il10", "Il6st", "Ccr7", "Il9", "Ikzf3", "Slamf6", "Sell")

htmap(atacRES_ordered_all, counts_tbl = counts_tbl_norm,
      highlight_genes = gene_list_Fig.6F, fig_num = "6F", 
      suffix = "_htmap_npTh17_ATACseq_DA_peaks_all", use_row_annot = TRUE)


## Fig 6G: Pathogenic Th17 signature score boxplot, ATAC-seq ---------
pathogenic_sig <- read.table("data/pathogenic_th17_signature.txt")$V1

sig_score_boxplot(counts_tbl = counts_tbl_norm, meta = metadata, 
                  anno = peakAnno_all@anno, sig = pathogenic_sig, fig_num = "6G", 
                  min_FDR = FALSE, title = "Pathogenic Th17 Signature",
                  suffix = "_norm_all")
sig_score_boxplot(counts_tbl = counts_tbl_norm, meta = metadata, 
                  anno = peakAnno_all@anno, sig = pathogenic_sig, fig_num = "6G", 
                  min_FDR = TRUE, title = "Pathogenic Th17 Signature",
                  suffix = "_norm_all")


## Fig 7E: RORgt-TCF1 shared genes signature score boxplot, ATAC-seq ---------
Rorgt_TCF1 <- c("Bcl11b", "Bcl2", "Il17a", "Il17f", "Il2ra", 
                "Il7r", "Pparg", "Rara", "Rora", "Rorc")

sig_score_boxplot(counts_tbl = counts_tbl_norm, meta = metadata,
                  anno = peakAnno_all@anno, sig = Rorgt_TCF1, fig_num = "7E",
                  min_FDR = FALSE, title = "Rort-TCF1 Bound Genes",
                  suffix = "_norm_all")
sig_score_boxplot(counts_tbl = counts_tbl_norm, meta = metadata,
                  anno = peakAnno_all@anno, sig = Rorgt_TCF1, fig_num = "7E",
                  min_FDR = TRUE, title = "Rort-TCF1 bound Genes",
                  suffix = "_norm_all")
