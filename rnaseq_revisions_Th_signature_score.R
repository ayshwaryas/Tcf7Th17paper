library(tidyverse)
library(DESeq2)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(ComplexHeatmap)
library(ggtext)
library(readxl)

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


# 2. Run DESeq2 ----------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_tbl), colData = metadata_tbl,
  design = ~ genotype)
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

res <- lfcShrink(dds, coef = "genotype_KO_vs_WT", type = "apeglm")

res_ordered <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  left_join(counts[, c("gene_id", "gene_name")], by = "gene_id") %>%
  dplyr::select(gene_id, gene_name, everything()) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up in KO", "down in KO")) %>%
  mutate(direction = factor(direction, c("up in KO", "down in KO"))) %>%
  arrange(direction, padj)

counts_tbl_norm <- DESeq2::counts(dds, normalized = TRUE)


# 3. Comparing Th signature score between TCF1cKO and WT (signature score volcano plot) ---------------------------
## (1) Th1 (2) Th2 (3) pathogenic Th17 (4) non-pathogenic Th17  (5) SLAMF6+ Th17 (6) CXCR6+ Th17

## (1) Th1 and (2) Th2 signatures ---------
Th_sig_names <- excel_sheets("data/Th_signatures.xlsx")
Th_sig_labels <- sapply(Th_sig_names, function(x) {
  x <- str_replace_all(x, "_", " ")
  if(str_detect(x, "[Uu]p$")) {str_remove(x, " [Uu]p$") %>% paste("Up in", .)}
  else {x}
})

Th_sig_raw <- sapply(Th_sig_names, function(x) {
  c(read_xlsx("data/Th_signatures.xlsx", sheet = x, col_names = FALSE)[, 1][[1]])},
  simplify = FALSE, USE.NAMES = TRUE) %>%
  `names<-`(Th_sig_labels)

mouse_human_genes_raw <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
convert_human_to_mouse <- function(gene_list){
  mouse_human_genes_raw %>%
    filter(DB.Class.Key %in% subset(mouse_human_genes_raw, Symbol %in% gene_list )$DB.Class.Key) %>% 
    filter(Common.Organism.Name == "mouse, laboratory") %>%
    mutate(Symbol = ifelse(is.na(Symbol), str_to_title(Symbol), Symbol)) %>%
    pull(Symbol) %>% unique
}

Th_sig <- lapply(Th_sig_raw, convert_human_to_mouse)


## (3) Pathogenic Th17 signature
pTh17_sig <- read.table("data/pathogenic_th17_signature.txt")$V1

## (4) Non-pathogenic Th17 signature
npTh17_sig <- c("Il6st",  "Il1rn", "Ikzl3", "Maf", "Ahr", "Il9", "Il10")

## (5) SLAMF6+ Th17 signature (top 100 up-regulated DEGs of CXCR6+ Th17 vs SLAMF6+ Th17, ordered by log2 FC)
## (6) Cxcr6+ Th17 signature (top 100 down-regulated DEGs of CXCR6+ Th17 vs SLAMF6+ Th17, ordered by log2 FC)
res_CXCR6_vs_SLAMF6 <- readxl::read_xlsx("data/CXCR6_vs_SLAMF6_Bulk RNAseq.xlsx") %>%
  `colnames<-`(c("gene_name", colnames(.)[-1]))

CXCR6_vs_SLAMF6_sig_logFC <- res_CXCR6_vs_SLAMF6 %>%
  filter(FDR < 0.05) %>%
  mutate(direction = ifelse(logFC > 0, "CXCR6+ Th17 signature", "SLAMF6+ Th17 signature")) %>%
  group_by(direction) %>%
  arrange(desc(abs(logFC))) %>%
  dplyr::slice(1:100) %>%
  split(f = .$direction) %>%
  lapply(pull, "gene_name")


## Compute signature score in TCF1cKO and WT
signatures <- c(list(`pTh17 signature` = pTh17_sig, `npTh17 signature` = npTh17_sig),
                Th_sig[1:2], CXCR6_vs_SLAMF6_sig_logFC)

Th_sig_score <- sapply(names(signatures), function(x) {
  sig_score_boxplot(counts_tbl_norm, metadata, sig = signatures[[x]], fig_num = "",
                    title = str_replace_all(x, "_", " "), suffix = "_norm", return_score = TRUE)}, simplify = FALSE, USE.NAMES = TRUE) %>%
  data.table::rbindlist(idcol = "signature")

## Compute log2 FC of signature score between TCF1cKO and WT samples, then calculate the
## BH-adjusted t-test p-values comparing the signature scores between the two genotypes
Th_sig_score_log2FC <- Th_sig_score %>%
  mutate(genotype = factor(genotype, c("KO", "WT"))) %>%
  split(f = .$signature) %>%
  lapply(function(df) {
    sig_score_WT <- subset(df, genotype == "WT")$sig_score
    sig_score_KO <- subset(df, genotype == "KO")$sig_score
    log2FC <- log2(mean(sig_score_KO) + 1) - log2(mean(sig_score_WT) + 1)
    ttest <- t.test(sig_score ~ genotype, data = df)
    data.frame(pval = ttest$p.value,
               statistic = ttest$statistic,
               log2FC = log2FC)
  }) %>%
  data.table::rbindlist(idcol = "signature")

Th_sig_score_log2FC$FDR <- p.adjust(Th_sig_score_log2FC$pval, method = "BH")

Th_sig_score_log2FC %>%
  ggplot(aes(x = log2FC, y = -log10(FDR))) +
  geom_vline(xintercept = 0, lty = 2, color = "grey60") + 
  geom_hline(yintercept = -log10(0.05), lty = 2, color = "grey60") + 
  geom_point(color = "red", size = 1.5) +
  geom_text_repel(aes(label = signature), lineheight = 0.8,
                  box.padding = 0.2, size = 4) +
  scale_y_continuous(sec.axis = dup_axis(name = NULL, breaks = -log10(0.05), label = "FDR=0.05")) +
  scale_color_identity() +
  labs(x = "Log2 Fold Change, TCF1cKO vs TCF1 WT", y = "-log<sub>10</sub>(FDR)") + 
  theme_cowplot(font_size = 12) +
  theme(axis.title.y = element_markdown(),
        axis.line.y.right = element_blank())
ggsave("figures/RNAseq_Th_signature_score_TCF1cKOvsWT.png", width = 5, height = 4)
ggsave("figures/RNAseq_Th_signature_score_TCF1cKOvsWT.pdf", width = 5, height = 4)
