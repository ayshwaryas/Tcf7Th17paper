library(tidyverse)
library(DESeq2)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(ComplexHeatmap)
library(ChIPseeker)
library(data.table)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
gene_list_Fig.6E <- c("Gzmb", "Rora", "Ccr2", "Il23r", "Icos", "Pparg",
               "Ccl1", "Cxcr6", "Gpr65", "Ccr3", "Ccr9", "Csf2",
               "Ccl3", "Ccr8", "Ccr6", "Il6ra", "Il17a", "Rorc",
               "Ccr1", "Tigit", "Sell", "Slamf6", "Bach2", "S1pr1",
               "Foxo1", "Il6st", "Il10", "Ccr7")

# 1. Load Data ----------------------
counts <- read.table("data/ATACseq/atac_TCF1cKO_counts_min_peak_length_250bp.txt",
                     sep = "\t", header = TRUE)
metadata <- read.csv("data/ATACseq/atac_metadata.csv") 

counts_tbl <- counts %>% dplyr::select(contains("dLN")) %>% `rownames<-`(counts$GeneID)
metadata_tbl <- metadata %>% column_to_rownames("Sample_ID")

all(metadata$Sample_ID == colnames(counts_tbl))

load("R_objects/Tcf7_Th17_ATACseq_DA_res.RData")
load("R_objects/PeakAnno.RData")
counts_tbl_norm <- DESeq2::counts(atacDDS, normalized = TRUE)

peaks_CXCR6_vs_SLAMF6 <- read.csv("data/ATACseq/CXCR6_vs_SLAMF6_all_peaks.csv") %>%
  mutate(name = paste0("peak.pub_", row_number())) %>%
  mutate(seqnames = str_extract(X, "chr.*(?=:)"),
         start = str_extract(X, "(?<=:)[0-9]+(?=\\-)"),
         end = str_extract(X, "(?<=\\-)[0-9]+$")) %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>%
  mutate(length = end - start + 1) %>%
  dplyr::select(-X) %>%
  dplyr::select(name, seqnames, start, end, length, everything()) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  mutate(direction = factor(direction, c("up", "down")))

# 2. TCF1cKO vs WT ATAC-seq DA peaks heatmap (FDR < 0.05) ------------------------------
## Direction of the overlapping peak in SLAMF6+ Th17 vs CXCR6+ Th16 ATAC-seq is annotated on the right side
## (yellow: up-regulated in SLAMF6+ Th17; black: down-regulated in SLAMF6+ Th17)

peaks_CXCR6_vs_SLAMF6.05 <- subset(peaks_CXCR6_vs_SLAMF6, padj < 0.05) %>%
  dplyr::select(name, seqnames, start, end, padj, log2FoldChange, direction) %>% 
  as.data.table %>% setkey(seqnames, start, end)

peaks_TCF1.05_minFDR <- subset(atacRES_ordered_all, padj < 0.05) %>%
  dplyr::select(GeneID, seqnames, start, end, padj, log2FoldChange, direction, SYMBOL) %>% 
  group_by(SYMBOL) %>%
  arrange(padj) %>% dplyr::slice(1) %>% ungroup %>%
  as.data.table %>% setkey(seqnames, start, end)

ATAC.overlap.05_minFDR <- foverlaps(
  peaks_TCF1.05_minFDR, peaks_CXCR6_vs_SLAMF6.05,
  type = "any", mult = "all", nomatch = NA) %>%
  group_by(GeneID) %>%
  arrange(padj) %>%
  dplyr::slice(1) %>%
  mutate(dir = paste(i.direction, direction, sep = ".")) %>%
  mutate(dir = factor(dir, c("up.up", "up.down", "up.NA", "down.NA", "down.up", "down.down")))

htmap_df <- counts_tbl_norm[ATAC.overlap.05_minFDR$GeneID, ] %>%
  as.data.frame() 
all(rownames(htmap_df) == ATAC.overlap.05_minFDR$GeneID)

htmap_df_scaled <- htmap_df %>%
  t() %>% scale() %>% t()

## Row annotation ----------------------------------------

row_annot <- rowAnnotation(
  `Cxcr6+ vs Slamf6+` = ATAC.overlap.05_minFDR$direction,
  link = anno_mark(
    at = which(ATAC.overlap.05_minFDR$SYMBOL %in% gene_list_Fig.6E),
    labels = ATAC.overlap.05_minFDR$SYMBOL[ATAC.overlap.05_minFDR$SYMBOL %in% gene_list_Fig.6E],
    labels_gp = gpar(fontsize = 10, fontface = 3),
    padding = unit(0.5, "mm"),
    link_width = unit(8, "mm")),
  col = list(`Cxcr6+ vs Slamf6+` = setNames(c("gold", "grey15"), c("up", "down"))),
  simple_anno_size = unit(3.5, "mm"),
  gap = unit(1, "points"),
  annotation_legend_param = list(
    `Cxcr6+ vs Slamf6+` = list(title = gt_render("CXCR6<sup>+</sup> vs SLAMF6<sup>+</sup>"),
                               title_gp = gpar(fontsize = 10, fontface = 2, lineheight = 1.15),
                               labels = gt_render(c("Up in CXCR6<sup>+</sup>", "Up in SLAMF6<sup>+</sup>")))),
  show_annotation_name = FALSE)

## Column annotation -------------------------------------
col_labels <- str_extract(colnames(htmap_df_scaled), "[0-9]")
col_split <- str_extract(colnames(htmap_df_scaled), "WT|KO") %>%
  factor(levels = c("WT", "KO"), labels = c("TCF1 WT", "TCF1cKO"))

col_annot <- columnAnnotation(
  genotype = col_split,
  col = list(genotype = setNames(c("black", "maroon"), c("TCF1 WT", "TCF1cKO"))),
  simple_anno_size = unit(4, "mm"),
  show_annotation_name = FALSE,
  show_legend = FALSE)

## Heatmap palette -----------------------------------------
htmap_col <- circlize::colorRamp2(
  c(min(htmap_df_scaled, na.rm=T), mean(htmap_df_scaled, na.rm=T),
    max(htmap_df_scaled, na.rm=T)),
  c("navy", "white", "firebrick3"), transparency = 0.05)

## Draw heatmap --------------------------------------------
p <- htmap_df_scaled %>%
  Heatmap(
    col = htmap_col,
    border_gp = gpar(color = "grey15"),
    row_split = ATAC.overlap.05_minFDR$dir,
    row_names_gp = gpar(fontsize = 9, fontface = "italic"),
    show_row_names = FALSE,
    row_title = NULL,
    column_labels = col_labels,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    show_column_names = FALSE,
    column_split = col_split,
    cluster_row_slices = FALSE,
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    right_annotation = row_annot,
    top_annotation = col_annot,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    name = "Z-score",
    row_gap = unit(c(1.2, 1.2, 2, 1.2, 1.2), "mm"),
    column_gap = unit(2, "mm")) 

filename <- "figures/revisions/heatmap_TCF1cKOvsWT_DApeaks_minFDR"
png(paste0(filename, ".png"), width = 5, height = 6.5, res = 400, units = "in")
draw(p, merge_legend = TRUE, padding = unit(c(2, 2, 2, -5), "mm"))
dev.off()
pdf(paste0(filename, ".pdf"), width = 5, height = 6.5)
draw(p, merge_legend = TRUE, padding = unit(c(2, 2, 2, -5), "mm"))
dev.off()

## Saving the order of the genes associated with DA peaks -------------------------------
p <- draw(p)
row_order <- row_order(p)

save(row_order, file = "results/2024-01-29.row_order_heatmap_TCF1cKOvsWT_DApeaks_minFDR.RData")
lapply(row_order, function(x) {
  data.frame(GeneID = rownames(htmap_df_scaled)[x])}) %>%
  data.table::rbindlist(idcol = "block") %>%
  separate(block, sep = "\\.", into = c("TCF1cKO_vs_WT", "CXCR6_vs_SLAMF6")) %>% 
  mutate(TCF1cKO_vs_WT = paste(TCF1cKO_vs_WT, "in TCF1cKO")) %>%
  mutate(CXCR6_vs_SLAMF6 = ifelse(is.na(CXCR6_vs_SLAMF6),  "Not DA", paste(CXCR6_vs_SLAMF6, "in CXCR6"))) %>%
  left_join(ATAC.overlap.05_minFDR %>% select(GeneID, SYMBOL), by = "GeneID") %>%
  write.csv("results/2024-01-29.heatmap_TCF1cKOvsWT_DApeaks_minFDR.csv", row.names = FALSE)
