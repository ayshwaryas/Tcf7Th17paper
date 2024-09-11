#### Network Analysis ####
## 1. Identify Targets:
##    (1) select genes associated with at least one DA peaks
##    (2) if multiple peaks were associated to the same Target, select the most significantly DA peak (smallest adjusted p-value)
##    (3) annotate genes by direction of DA and DE
## 2. Identify TFs:
##    (1) motif scan in selected peaks for Targets
## 3. Build and visualize network:
##    (1) 
##    (2) create igraph object
##    (3) identify communities and visualize as heatmap and network
##    (4) trim network

rm(list = ls())
setwd("/singerlab/linglin/Th17_single_cell_eae_ut/Davide_Tcf7/")
library(tidyverse)
library(GenomicRanges)
library(memes)
library(universalmotif)
library(BSgenome.Mmusculus.UCSC.mm10)
library(igraph)
library(pheatmap)
library(dendextend)
library(RColorBrewer)
library(openxlsx)
library(dynamicTreeCut)

direction_vec <- c("all", "KO_KO", "KO_n.s.", "KO_WT", "WT_KO", "WT_n.s.", "WT_WT")

#### 1. Identify Targets ####
# peak to gene association has been performed by Yufan and Ayshwarya
DA_df <- rbind(
  read.csv("DA_peaks_genes_overlap=all_UP_with_info.csv"),
  read.csv("DA_peaks_genes_overlap=all_DN_with_info.csv")
)
DA_df <- DA_df %>% 
  filter(!is.na(SYMBOL)) %>% 
  filter(abs(log2FoldChange) > log2(1.2)) %>% 
  group_by(SYMBOL) %>% 
  mutate(n_peaks = n(),
         is_min_padj_ATAC = seq_along(padj) == which.min(padj)) %>% 
  filter(is_min_padj_ATAC) # only keep the peak with smallest adjusted p-value
table(DA_df$n_peaks)

# DE has been performed by Yufan and Ayshwarya
DE_df <- read.csv("DGE_results_dLN_KO_vs_WT_ordered_by_padj.csv")
DE_df$direction <- "n.s."
DE_df$direction[(DE_df$padj < 0.05) & (DE_df$log2FoldChange > 0)] <- "up"
DE_df$direction[(DE_df$padj < 0.05) & (DE_df$log2FoldChange < 0)] <- "down"

target_df <- merge(DA_df, DE_df, 
                   by.x = "SYMBOL", by.y = "gene_symbol", 
                   all.x = TRUE, all.y = FALSE, # some genes were not in DE results (maybe not detected) 
                   suffixes = c("_ATAC", "_RNA")) %>% 
  # peak level direction of regulation
  mutate(direction_RNA_parsed = case_when(is.na(direction_RNA) ~ "n.s.",
                                          direction_RNA == "up" ~ "KO",
                                          direction_RNA == "down" ~ "WT",
                                          direction_RNA == "n.s." ~ "n.s."),
         direction_ATAC_parsed = case_when(direction_ATAC == "up" ~ "KO",
                                           direction_ATAC == "down" ~ "WT"),
         direction_ATAC_RNA = paste(direction_ATAC_parsed, direction_RNA_parsed, sep = "_")) %>% 
  mutate(direction_reg = case_when(direction_ATAC_RNA %in% c("WT_WT", "KO_KO") ~ "Activation",
                                   direction_ATAC_RNA %in% c("WT_KO", "KO_WT") ~ "Repression",
                                   direction_ATAC_RNA %in% c("WT_n.s.", "KO_n.s.") ~ "Poised"))
target_gr <- GRanges(target_df)

table(target_df$direction_ATAC, target_df$direction_RNA)

#### 2. Identify TFs ####
# (1) motif scan using FIMO
# load motif database
meme_db_path <- "~/utils/meme-5.4.1/motif_databases/MOUSE/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme"
meme_db <- read_meme(meme_db_path) %>% universalmotif::to_df()
motif_tf_map <- read.table("MOUSE_mono_motifs_20220517.tsv", sep = "\t", header = TRUE)
altname <- motif_tf_map$Transcription.factor[match(meme_db$name, motif_tf_map$Model)] %>% sub(".*:", "", .)
altname[is.na(altname)] <- sub("_.*", "", meme_db$name[is.na(altname)])
meme_db$altname <- altname

# filter motifs by expression: >=1TPM in >=1 samples
tpm <- read.csv("Tcf7_Th17_Davide_tpm.csv", as.is = TRUE)
expressed_genes <- unique(tpm$gene_name[rowMeans(tpm[,grep("dLN", colnames(tpm))]) >= 1])
setdiff(expressed_genes, DE_df$gene_symbol[!is.na(DE_df$padj)]) %>% intersect(meme_db$altname)
setdiff(DE_df$gene_symbol[!is.na(DE_df$padj)], expressed_genes) %>% intersect(meme_db$altname)
meme_db_sub <- meme_db %>% 
  filter(altname %in% expressed_genes)
meme_db_ls <- universalmotif::to_list(meme_db_sub, extrainfo = FALSE)

## scan for motifs
target_flank <- target_gr %>%
  `values<-`(value = NULL) %>%
  plyranges::anchor_center() %>%
  plyranges::mutate(width = 200)
seqlevelsStyle(target_flank) <- "UCSC"
saveRDS(target_flank, file = "target_flank.rds")
fimo_res <- target_flank %>%
  get_sequence(BSgenome.Mmusculus.UCSC.mm10) %>%
  runFimo(motifs = meme_db_ls, thresh = 1e-4, parse_genomic_coord = F, meme_path = "~/meme/bin/")
saveRDS(fimo_res, file = "fimo_res.rds")

## merge FIMO results with DA_DE results
fimo_res_df <- data.frame(fimo_res)
colnames(fimo_res_df) <- paste0(colnames(fimo_res_df), "_FIMO")
target_df$peak_flank <- with(data.frame(target_flank), paste0(seqnames, ":", start, "-", end))
tf_target_df <- merge(target_df, fimo_res_df,
                      by.x = "peak_flank", by.y = "seqnames_FIMO")

#### 3. Build network ####
# clean and parse DA, DE and motif scan results
graph_df <- tf_target_df %>%
  relocate(motif_alt_id_FIMO, SYMBOL) %>%
  # in each treatment, remove duplicated hits of the same TF in the same peak (keep the on with highest FIMO score)
  group_by(SYMBOL, motif_alt_id_FIMO, direction_ATAC_RNA, direction_reg) %>%
  summarise(n_hits = n()) %>% 
  mutate(direction_reg_num = ifelse(direction_reg == "Activation", 1, 
                                    ifelse(direction_reg == "Poised", 0, -1))) %>% 
  ungroup() %>% 
  filter(SYMBOL != "Tcf7") ## Tcf7 was knocked out
saveRDS(graph_df, file = "complete_graph_df.rds")
write.csv(graph_df, file = "complete_graph_df.csv", quote = FALSE)

# make heat maps to help identify clusters of TFs and targets
adj_mtx <- graph_df %>%
  dplyr::select(SYMBOL, direction_reg_num, motif_alt_id_FIMO) %>%
  spread(key = motif_alt_id_FIMO, value = direction_reg_num, fill = NA) %>% 
  column_to_rownames(var = "SYMBOL")

ann_tf <- data.frame(
  direction_ATAC = DA_df$direction[match(colnames(adj_mtx), DA_df$SYMBOL)],
  direction_RNA = DE_df$direction[match(colnames(adj_mtx), DE_df$gene_symbol)]
) %>% 
  mutate(direction_ATAC = factor(replace(direction_ATAC, is.na(direction_ATAC), "n.s."), 
                                 levels = c("down", "n.s.", "up"), labels = c("WT", "n.s.", "KO")),
         direction_RNA = factor(replace(direction_RNA, is.na(direction_RNA), "n.s."), 
                                levels = c("down", "n.s.", "up"), labels = c("WT", "n.s.", "KO"))) %>% 
  `rownames<-`(colnames(adj_mtx))
ann_target <- target_df %>% 
  filter(SYMBOL %in% rownames(adj_mtx)) %>% 
  dplyr::select(direction_ATAC, direction_RNA, SYMBOL) %>% 
  mutate(direction_ATAC = factor(direction_ATAC, levels = c("down", "n.s.", "up"), labels = c("WT", "n.s.", "KO")),
         direction_RNA = factor(replace(direction_RNA, is.na(direction_RNA), "n.s."), 
                                levels = c("down", "n.s.", "up"), labels = c("WT", "n.s.", "KO"))) %>% 
  column_to_rownames(var = "SYMBOL")
ann_colors <- list(
  "direction_ATAC" = c("WT" = "blue", "n.s." = "grey80", "KO" = "maroon"),
  "direction_RNA" = c("WT" = "blue", "n.s." = "grey80", "KO" = "maroon")
)
my_colors <- colorRampPalette(colors = c("blue", "gold", "red"))(3)
target_cluster_ls <- list()
wb <- createWorkbook()
for (i in direction_vec) {
  cat(i, "\n")
  if (i == "all") { # TF clusters
    idx <- 1:nrow(adj_mtx)
    pdf_height <- 6
    my_breaks <- NA
    show_rownames <- FALSE
    tmp <- -is.na(adj_mtx)
    cluster_tf <- hclust(dist(t(tmp)), method = "ward.D")
    k <- 10
    dend <- as.dendrogram(cluster_tf)
    pdf("dend_TF.pdf", width = 4, height = 16)
    dend %>% set("labels_col", value = brewer.pal(k, "Paired"), k = k) %>% 
      set("labels_cex", value = 0.4) %>% 
      set("branches_k_color", value = brewer.pal(k, "Paired"), k = k) %>% 
      plot(main = paste0("TF clusters (", k, " clusters)"), horiz = TRUE)
    dev.off()
    cluster_tf_label <- cutree(cluster_tf, k = k)
    col_gaps <- table(cluster_tf_label)[unique(cluster_tf_label[cluster_tf$order])] %>% cumsum
  } else {
    d_ATAC <- sub("_.*", "", i)
    d_RNA <- sub(".*_", "", i)
    ann_target_sub <- ann_target %>% filter(direction_ATAC == d_ATAC, direction_RNA == d_RNA)
    idx <- which(rownames(adj_mtx) %in% rownames(ann_target_sub))
    tmp <- -is.na(adj_mtx[idx,])
    pdf_height <- 1 + length(idx) / 10
    my_breaks <- seq(-1, 1, length.out = 3)
    show_rownames <- TRUE
    addWorksheet(wb, i)
  }
  
  distM <- dist(tmp)
  cluster_target <- hclust(distM, method = "ward.D2")
  cluster_target_label <- cutreeDynamic(cluster_target, minClusterSize = 3, method = "hybrid", distM = as.matrix(distM))
  names(cluster_target_label) <- cluster_target$labels
  
  target_cluster_ls[[i]] <- cluster_target_label[cluster_target$order]
  row_gaps <- table(cluster_target_label)[unique(cluster_target_label[cluster_target$order])] %>% cumsum
  
  cat("n =", length(idx), ", k =", k, ", height =", pdf_height, "\n")
  pdf(paste0("network_heatmap_", i, ".pdf"), width = 14, height = pdf_height)
  pheatmap(adj_mtx[idx[cluster_target$order], cluster_tf$order],
           color = my_colors,
           na_col = "grey90",
           breaks = my_breaks,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           gaps_col = col_gaps,
           gaps_row = row_gaps,
           annotation_row = ann_target,
           annotation_col = ann_tf,
           annotation_colors = ann_colors,
           show_rownames = show_rownames,
           fontsize_row = 4,
           fontsize_col = 4,
           legend = FALSE)
  dev.off()
  if (pdf_height > 20) {
    pdf(paste0("network_heatmap_", i, "_compressed.pdf"), width = 14, height = 6)
    pheatmap(adj_mtx[idx[cluster_target$order], cluster_tf$order],
             color = my_colors,
             na_col = "grey90",
             breaks = my_breaks,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             gaps_col = col_gaps,
             gaps_row = row_gaps,
             annotation_row = ann_target,
             annotation_col = ann_tf,
             annotation_colors = ann_colors,
             show_rownames = FALSE,
             fontsize_row = 4,
             fontsize_col = 4,
             legend = FALSE)
    dev.off()
  }
  if (i != "all") {
    writeData(wb, i, 
              data.frame(ann_target[rownames(adj_mtx)[idx[cluster_target$order]],],
                         cluster = target_cluster_ls[[i]],
                         0 + !is.na(adj_mtx[idx[cluster_target$order], cluster_tf$order])), 
              rowNames = TRUE, colNames = TRUE)
  }
}
saveRDS(target_cluster_ls, file = "target_cluster_ls.rds")
saveWorkbook(wb, file = "adj_matrix.xlsx", overwrite = TRUE)

#### 4. Make representative (summarized) network and prep data for visualization in Cytoscape ####
graph_df <- readRDS(file = "complete_graph_df.rds")
motif_tf_map <- read.table("MOUSE_mono_motifs_20220517.tsv", sep = "\t", header = TRUE)
## prune network based on selected targets
selected_target <- read.xlsx("TCF7_Th17_TF network analysis_new.xlsx", sheet = "Target_processed")
stopifnot(length(setdiff(selected_target$gene, graph_df$SYMBOL)) == 0)
stopifnot(all(graph_df$direction_ATAC_RNA[match(selected_target$gene, graph_df$SYMBOL)] == selected_target$ATAC_RNA))
## prune network based on selected TFs
selected_tf <- read.xlsx("TCF7_Th17_TF network analysis.xlsx", sheet = "TF_processed")
motif_tf_map <- motif_tf_map %>%
  mutate(Model_short = sub("_.*", "", Model),
         TF_short = sub(".*:", "", Transcription.factor))
selected_tf$altname <- motif_tf_map$TF_short[match(selected_tf$TF, motif_tf_map$Model_short)]
selected_tf$altname[selected_tf$TF == "TCF1"] <- "Tcf7"
selected_tf$altname[selected_tf$TF == "HEB"] <- "Tcf12"
stopifnot(all(selected_tf$altname == selected_tf$symbol_mgi))


## summarize module connectivity
adj_summ_ls <- list()
target_cluster_anno_ls <- list()
thresh_prop_hits <- 0.3
for (i in direction_vec) {
  if (i %in% c("all", "WT_KO", "KO_WT")) next
  cluster_adj <- read.xlsx("adj_matrix.xlsx", sheet = i, rowNames = TRUE)
  adj <- cluster_adj[,-(1:3)]
  tf_vec <- intersect(selected_tf$altname, colnames(adj))
  adj_sub <- adj[,tf_vec]
  
  n_total <- nrow(adj)
  n_total_by_cluster <- table(paste0("Cluster_", cluster_adj$cluster))
  n_hits <- colSums(adj_sub)
  n_hits_by_cluster <- rowsum(adj_sub, paste0("Cluster_", cluster_adj$cluster))
  
  p <- sapply(rownames(n_hits_by_cluster), function(a) {
    sapply(colnames(n_hits_by_cluster), function(b) {
      phyper(n_hits_by_cluster[a, b], n_total_by_cluster[a], n_total - n_total_by_cluster[a], n_hits[b], lower.tail = FALSE)
    })
  }) %>% t %>% as.data.frame(row.names = rownames(n_hits_by_cluster), col.names = colnames(n_hits_by_cluster))
  
  adj_summ <- (n_hits_by_cluster / n_total_by_cluster)
  adj_summ[(p > 0.05) & (adj_summ < thresh_prop_hits)] <- NA
  adj_summ_ls[[i]] <- adj_summ
  
  cur_targets <- selected_target$gene[selected_target$direction == i]
  stopifnot(all(cur_targets %in% rownames(cluster_adj)))
  target_cluster_name <- cluster_adj %>% 
    dplyr::select(cluster) %>% 
    rownames_to_column(var = "target") %>% 
    filter(target %in% cur_targets) %>% 
    group_by(cluster) %>% 
    summarise(cluster_name = paste(target, collapse = ";")) %>% 
    mutate(cluster = paste0("Cluster_", cluster))
  target_cluster_anno_ls[[i]] <- target_cluster_name
} 

graph_summ <- lapply(names(adj_summ_ls), function(i) {
  adj_summ_ls[[i]] %>% 
    data.frame() %>% 
    rownames_to_column(var = "target_cluster") %>% 
    gather(key = "TF_altname", value = "phits", -target_cluster) %>% 
    filter(!is.na(phits)) %>% 
    mutate(TF = selected_tf$TF[match(TF_altname, selected_tf$altname)],
           direction_target = i, 
           target = target_cluster_anno_ls[[i]]$cluster_name[match(target_cluster, target_cluster_anno_ls[[i]]$cluster)],
           direction_tf = selected_tf$direction[match(TF_altname, selected_tf$altname)],
           direction_reg = ifelse(direction_target %in% c("WT_WT", "KO_KO"), "Activation", "Poised")) %>% 
    filter(!is.na(target)) %>% 
    relocate(TF, target)
}) %>% Reduce(rbind, .)
graph_summ <- graph_summ %>% 
  group_by(TF) %>% 
  add_tally(name = "degree_TF") %>% 
  ungroup() %>% 
  group_by(target) %>% 
  add_tally(name = "degree_target") %>% 
  ungroup()
write.csv(graph_summ, file = paste0("graph_target_summarized_", thresh_prop_hits, "_", Sys.Date(), ".csv"), quote = FALSE, row.names = FALSE)
