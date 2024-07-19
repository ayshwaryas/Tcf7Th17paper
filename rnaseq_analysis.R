library(tidyverse)
library(DESeq2)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(ComplexHeatmap)

# 1. Load Data ----------------------
counts <- read.csv("data/RNAseq/Tcf7_Th17_RNAseq_count.csv")
tpm <- read.csv("data/RNAseq/Tcf7_Th17_RNAseq_tpm.csv")
metadata <- read.csv("data/RNAseq/Tcf7_Th17_RNAseq_metadata.csv") %>%
  mutate(genotype = factor(genotype, c("WT", "KO")))

counts_tbl <- counts[, -2] %>% column_to_rownames("gene_id") 
tpm_tbl <- tpm[, -2] %>% column_to_rownames("gene_id")
metadata_tbl <- metadata %>% column_to_rownames("Sample_ID")

all(metadata$Sample_ID == colnames(counts_tbl))
all(metadata$Sample_ID == colnames(tpm_tbl))


# 2. Differential gene expression analysis using DESeq2 ----------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_tbl), colData = metadata_tbl,
  design = ~ genotype)
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

res <- results(dds, contrast = c("genotype", "KO", "WT"))
res_ordered <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  left_join(counts[, c("gene_id", "gene_name")], by = "gene_id") %>%
  dplyr::select(gene_id, gene_name, everything()) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up in KO", "down in KO")) %>%
  mutate(direction = factor(direction, c("up in KO", "down in KO"))) %>%
  arrange(direction, padj)

counts_tbl_norm <- DESeq2::counts(dds, normalized = TRUE)

# 3. Prepare Figures ------------------------------------------
## Function for generating heatmaps
htmap <- function(DGE_res, highlight_genes = NULL, width = 5.5, height = 8, gap = c(2.5, 1.5),
                  fig.num = "", suffix = "", use_row_annot = FALSE, thres = 0.05) {
  
  DEG <- DGE_res %>% filter(padj < thres) 
  
  htmap_df <- tpm %>%
    filter(rowSums(.[, -(1:2)]) != 0) %>%
    inner_join(DEG %>% dplyr::select(gene_id, direction), by = "gene_id") %>%
    mutate(direction = factor(direction, c("up in KO", "down in KO"))) %>%
    arrange(direction)
  
  htmap_df_scaled <- htmap_df %>%
    dplyr::select(contains("dLN")) %>%
    t() %>% scale() %>% t()
  
  samples <- colnames(htmap_df_scaled)
  col_labels <- str_extract(samples, "[0-9]+")
  col_split <- str_extract(samples, "WT|KO") %>% 
    factor(levels = c("WT", "KO"))
  
  htmap_col <- circlize::colorRamp2(
    c(min(htmap_df_scaled, na.rm=T), mean(htmap_df_scaled, na.rm=T),
      max(htmap_df_scaled, na.rm=T)),
    c("navy", "white", "firebrick3"), transparency = 0.05)
  
  if(use_row_annot) {
    row_annot <- rowAnnotation(
      link = anno_mark(at = which(htmap_df$gene_name %in% highlight_genes),
                       labels = htmap_df$gene_name[htmap_df$gene_name %in% highlight_genes],
                       labels_gp = gpar(fontsize = 9, fontface = "italic"), 
                       padding = unit(1, "mm")))
    show_row_names <- FALSE; row_labels <- NULL
  }
  else {
    row_annot <- NULL
    show_row_names <- TRUE
    row_labels <- htmap_df$gene_name
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
      row_split = htmap_df$direction,
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
  
  png(paste0("figures/Fig.", fig.num, suffix, "_", thres, ".png"), 
      width = width, height = height, res = 500, units = "in")
  draw(p)
  dev.off()
  
  pdf(paste0("figures/Fig.", fig.num, suffix, "_", thres, ".pdf"), 
      width = width, height = height)
  draw(p)
  dev.off()
}

## Function for signature score computation, testing and visualization (boxplot)
## Signature score = average of DESeq2 normalized counts
sig_score_boxplot <- function(expr_mat, meta, sig, fig_num, 
                              title = "", suffix = "", nudge_x = 0) {
  check_genes <- setdiff(str_to_upper(sig), str_to_upper(tpm$gene_name))
  
  score_df <- data.frame(expr_mat) %>%
    rownames_to_column("gene_id") %>%
    left_join(tpm[, c("gene_id", "gene_name")], by = "gene_id") %>%
    dplyr::select(-gene_id) %>%
    dplyr::select(gene_name, everything()) %>%
    filter(str_to_upper(gene_name) %in% str_to_upper(sig)) %>%
    pivot_longer(2:ncol(.), names_to = "Sample_ID", values_to = "TPM")  %>%
    group_by(Sample_ID) %>%
    summarise(sig_score = mean(TPM)) %>%
    left_join(meta, by = "Sample_ID") %>%
    mutate(genotype = factor(genotype, c("WT", "KO"))) 
  
  sig_score_WT <- subset(score_df, genotype == "WT")$sig_score
  sig_score_KO <- subset(score_df, genotype == "KO")$sig_score
  
  t_test_pval <- t.test(sig_score_WT, sig_score_KO)$p.value
  
  score_df %>%
    ggplot(aes(x = genotype, y = log2(sig_score), 
               color = genotype, shape = genotype)) +
    geom_boxplot(width = 0.3, outlier.alpha = 0, size = 0.45, fill = NA) +
    geom_point(size = 1.5, position = position_nudge(x = nudge_x, y = 0)) +
    scale_color_manual(values = c("black", "maroon")) +
    labs(x = NULL, y = "Log2 Signature Score", title = title) +
    theme_cowplot(font_size = 10) +
    guides(color = guide_legend(keywidth = 0.8, keyheight = 0.8)) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank())
  
  file_path <- paste0("figures/Fig.", fig_num, "_sig_score_RNAseq_", 
                      str_replace_all(title, " ", "_"), suffix)
  
  ggsave(paste0(file_path, ".png"), width = 4, height = 3.2)
  ggsave(paste0(file_path, ".pdf"), width = 4, height = 3.2)
  
  return(list(check_genes = check_genes, t_test_pval = round(t_test_pval, 5)))
}


## Fig 6A: Heatmap of npTh17 RNAseq DEGs ---------
gene_list_Fig.6A <- c("Rorc", "Il23r", "Il7r", "Cxcr6", "Ccl5", "Icos", "Lag3",
                      "Il6st", "Slamf6", "Ccr7")

htmap(res_ordered, highlight_genes = gene_list_Fig.6A, fig.num = "6A", 
      suffix = "_htmap_npTh17_RNAseq_DEGs", thres = 0.05)

## Fig 6B: Pathogenic Th17 signature score boxplot, RNA-seq ---------
pathogenic_sig <- read.table("data/pathogenic_th17_signature.txt")$V1
sig_score_boxplot(counts_tbl_norm, metadata, sig = pathogenic_sig, fig_num = "6B",
                  title = "Pathogenic Th17 Signature", suffix = "_norm")

## Fig 7B: Heatmap of npTh17 RNAseq DEGs bounded to Rorgt -------
Rorgt_ChIPseq <- read.csv("data/Rorgt Ch-IP Seq peak regions.csv", check.names = FALSE) %>%
  dplyr::select("5' proximal", "Gene body + 5'", "3' proximal") %>%
  apply(2, function(x) { paste(x[x!= "-1"], collapse = ";")  }) %>%
  lapply(str_split, pattern = ";") %>%
  lapply(unlist)

Rorgt_ChIPseq_all <- unlist(Rorgt_ChIPseq) %>% unique()

overlap_DEG_Rorgt_ChIP <- res_ordered %>%
  mutate(gene_name_upper = str_to_upper(gene_name)) %>%
  filter(gene_name_upper %in% Rorgt_ChIPseq_all)

gene_list_Fig.7B <- c("Rora", "Pparg", "Il23r", "Rorc", "Foxo1")

htmap(overlap_DEG_Rorgt_ChIP, fig.num = "7B", highlight_genes = gene_list_Fig.7B, 
      suffix = "_htmap_npTh17_RNAseq_DEGs_Rorgt_bound_genes")


## Fig 7C: Rorgt core gene signature score boxplot, RNAseq -------
Rorgt_core_gene <- c("Il17a", "Il17f", "Il23r", "Ccl20", "Il1r1", "Ltb4r1", 
                     "2310007L2Rik", "Furin", "Fam124b", "Tmem176a", "Tmem176b")

sig_score_boxplot(counts_tbl_norm, metadata, sig = Rorgt_core_gene, fig_num = "7C",
                  title = "Rorgt core gene signature", suffix = "_norm")
