############################################################
# RNA-seq Analysis Pipeline (ERα / Nalm6 project)
# Description:
#   Salmon gene counts → QC → DE → GSEA → TF activity →
#   Stemness scores → Stemness marker expression vs Veh →
############################################################

## =========================
## 0. LIBRARIES
## =========================
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(tidyverse)
  library(pheatmap)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
  library(UpSetR)
  library(viper)
  library(dorothea)
  library(igraph)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggrepel)
})

# Force base/dplyr to win namespace conflicts from igraph/viper/IRanges
select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
intersect <- base::intersect
union     <- base::union
setdiff   <- base::setdiff

set.seed(1)

## =========================
## 1. PATHS
## =========================
# Update counts_file to point to your local salmon.merged.gene_counts.tsv
counts_file <- "data/salmon.merged.gene_counts.tsv"

outdir <- "results"
figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")

dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

## =========================
## 2. METADATA
## =========================
meta <- data.frame(
  sample    = c("N6_Veh_r1","N6_Veh_r2",
                "N6_E2_r1","N6_E2_r2",
                "N6_Tam_r1","N6_Tam_r2",
                "N6_TamE2_r1","N6_TamE2_r2"),
  condition = factor(c("Veh","Veh","E2","E2","Tam","Tam","TamE2","TamE2"),
                     levels = c("Veh","E2","Tam","TamE2")),
  rep       = factor(c("r1","r2","r1","r2","r1","r2","r1","r2")),
  batch     = factor(c("Batch1","Batch2","Batch1","Batch2",
                       "Batch1","Batch2","Batch1","Batch2"))
)

## =========================
## 3. LOAD COUNTS
## =========================
counts <- read.delim(counts_file, row.names = 1)
counts <- counts[, meta$sample]

## =========================
## 4. edgeR MODEL
## =========================
design <- model.matrix(~ rep + condition, data = meta)

y <- DGEList(counts = counts)
y <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust = TRUE)
fit    <- glmQLFit(y, design, robust = TRUE)
logcpm <- cpm(y, log = TRUE, prior.count = 2)

## =========================
## 5. GENE SYMBOL MAPPING
## =========================
# Salmon produces versioned ENSEMBL IDs (e.g. ENSG00000000419.14)
clean_ids <- sub("\\.\\d+$", "", rownames(logcpm))

symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys      = clean_ids,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

message(sprintf("Mapped %d / %d genes (%.1f%%)",
  sum(!is.na(symbols)), length(symbols), 100 * mean(!is.na(symbols))))

# Lookup table for annotating DE tables
id_map <- data.frame(
  ensembl_versioned = rownames(logcpm),
  ensembl           = clean_ids,
  symbol            = unname(symbols),
  stringsAsFactors  = FALSE
)

# Apply symbols; drop unmapped and duplicates
rownames(logcpm) <- unname(symbols)
logcpm <- logcpm[!is.na(rownames(logcpm)), ]
logcpm <- logcpm[!duplicated(rownames(logcpm)), ]

message(sprintf("logcpm: %d genes x %d samples", nrow(logcpm), ncol(logcpm)))

## =========================
## 5B. FILTER RIBOSOMAL & MITOCHONDRIAL GENES
## =========================
ribo_genes <- grep("^RP[LS]", rownames(logcpm), value = TRUE)
mito_genes <- grep("^MT-",    rownames(logcpm), value = TRUE)

message(sprintf("Removing %d ribosomal and %d mitochondrial genes",
                length(ribo_genes), length(mito_genes)))

genes_to_remove <- base::union(ribo_genes, mito_genes)
logcpm <- logcpm[!rownames(logcpm) %in% genes_to_remove, ]

message(sprintf("logcpm after ribo/mito removal: %d genes", nrow(logcpm)))

## =========================
## 6. BATCH CORRECTION (visualization only)
## =========================
logcpm_bc <- removeBatchEffect(
  logcpm,
  batch  = meta$batch,
  design = model.matrix(~ condition, data = meta)
)

logcpm_bc <- logcpm_bc[!rownames(logcpm_bc) %in% genes_to_remove, ]
message("Batch correction applied — logcpm_bc ready for visualization")

# Shared column annotation for pheatmap
ann_col <- data.frame(condition = meta$condition,
                      rep       = meta$rep,
                      row.names = meta$sample)

## =========================
## 7. QC FIGURES (Fig 1)
## =========================

plot_pca <- function(mat, title) {
  pca <- prcomp(t(mat), scale. = FALSE)
  pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  df  <- data.frame(pca$x, meta)
  ggplot(df, aes(PC1, PC2, color = condition, shape = rep)) +
    geom_point(size = 4) +
    labs(title = title,
         x = paste0("PC1 (", pct[1], "%)"),
         y = paste0("PC2 (", pct[2], "%)")) +
    theme_classic()
}

# A) Library size
lib_df <- data.frame(
  sample    = colnames(counts),
  libsize   = colSums(counts),
  condition = meta$condition
)

p_lib <- ggplot(lib_df, aes(sample, libsize, fill = condition)) +
  geom_col() +
  scale_y_continuous(labels = scales::comma) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(y = "Library size (raw counts)", x = NULL)

ggsave(file.path(figdir, "Fig1A_LibrarySize.png"), p_lib, width = 8, height = 4, dpi = 300)

# B) PCA before batch correction
p_pca_raw <- plot_pca(logcpm, "PCA — uncorrected")
ggsave(file.path(figdir, "Fig1B_PCA_uncorrected.png"), p_pca_raw, width = 6, height = 5, dpi = 300)

# C) PCA after batch correction
p_pca_bc <- plot_pca(logcpm_bc, "PCA — batch-corrected")
ggsave(file.path(figdir, "Fig1C_PCA_corrected.png"), p_pca_bc, width = 6, height = 5, dpi = 300)

# D) Sample correlation heatmap
png(file.path(figdir, "Fig1D_Correlation.png"), width = 800, height = 600)
pheatmap(cor(logcpm_bc), annotation_col = ann_col, main = "Sample correlation")
dev.off()

## =========================
## 8. DIFFERENTIAL EXPRESSION
## =========================
contrasts <- list(
  E2_vs_Veh    = makeContrasts(conditionE2,    levels = design),
  Tam_vs_Veh   = makeContrasts(conditionTam,   levels = design),
  TamE2_vs_Veh = makeContrasts(conditionTamE2, levels = design)
)

DE_list <- lapply(names(contrasts), function(nm) {
  tab <- topTags(glmQLFTest(fit, contrast = contrasts[[nm]]), n = Inf)$table
  tab$ensembl <- sub("\\.\\d+$", "", rownames(tab))
  tab$symbol  <- id_map$symbol[match(tab$ensembl, id_map$ensembl)]
  write.csv(tab, file.path(tabdir, paste0("DE_", nm, ".csv")))
  tab
})
names(DE_list) <- names(contrasts)

tab_E2 <- DE_list[["E2_vs_Veh"]]

## =========================
## 9. VOLCANO & EFFECT SIZE (Fig 2)
## =========================
tab_E2_plot <- tab_E2 %>%
  mutate(
    sig   = FDR < 0.05 & abs(logFC) > 1,
    label = ifelse(rank(-abs(logFC)) <= 15 & sig, symbol, NA)
  )

p_vol <- ggplot(tab_E2_plot, aes(logFC, -log10(FDR), color = sig)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_label_repel(aes(label = label), size = 3,
                   max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c("grey60","firebrick")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_classic() +
  labs(title = "E2 vs Veh", color = "FDR<0.05 & |logFC|>1")

ggsave(file.path(figdir, "Fig2A_Volcano.png"), p_vol, width = 7, height = 6, dpi = 300)

p_fc <- ggplot(tab_E2, aes(logFC)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  theme_classic() +
  labs(title = "Effect size distribution (E2 vs Veh)")

ggsave(file.path(figdir, "Fig2B_EffectSize.png"), p_fc, dpi = 300)

## =========================
## 10. HEATMAP — TOP 100 DE GENES (Fig 3)
## =========================
top_genes <- base::intersect(na.omit(tab_E2$symbol[1:100]), rownames(logcpm_bc))
message("Genes in heatmap: ", length(top_genes))

mat_z <- t(scale(t(logcpm_bc[top_genes, ])))

png(file.path(figdir, "Fig3_Heatmap.png"), width = 800, height = 1000)
pheatmap(mat_z,
         annotation_col = ann_col,
         show_rownames  = length(top_genes) <= 100,
         fontsize_row   = 7,
         main           = "Top 100 DE genes (E2 vs Veh)")
dev.off()

## =========================
## 11. GSEA — HALLMARK (Fig 4)
## =========================
msig     <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(msig$gene_symbol, msig$gs_name)

make_ranks <- function(tab) {
  df <- tab %>%
    dplyr::filter(!is.na(symbol)) %>%
    arrange(desc(abs(logFC))) %>%
    distinct(symbol, .keep_all = TRUE)
  r        <- df$logFC
  names(r) <- df$symbol
  sort(r, decreasing = TRUE)
}

gsea_list <- lapply(names(DE_list), function(nm) {
  ranks <- make_ranks(DE_list[[nm]])
  fg    <- fgsea(pathways, ranks, minSize = 15, maxSize = 500)
  fg    <- fg[order(fg$padj), ]
  write.csv(fg[, !sapply(fg, is.list)],
            file.path(tabdir, paste0("fgsea_", nm, ".csv")))
  fg
})
names(gsea_list) <- names(DE_list)

p_bar <- gsea_list[["E2_vs_Veh"]] %>%
  slice_min(padj, n = 20) %>%
  ggplot(aes(reorder(pathway, NES), NES, fill = NES > 0)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("steelblue","firebrick")) +
  coord_flip() +
  theme_classic() +
  labs(x = NULL, title = "GSEA Hallmark (E2 vs Veh)")

ggsave(file.path(figdir, "Fig4_GSEA_barplot.png"), p_bar, width = 9, height = 6, dpi = 300)

## =========================
## 12. UpSet — MULTI-CONTRAST DE OVERLAP (Fig 5)
## =========================
sets <- lapply(DE_list, function(tab)
  tab$symbol[tab$FDR < 0.05 & !is.na(tab$symbol)])

png(file.path(figdir, "Fig5_UpSet.png"), width = 900, height = 700)
upset(fromList(sets), order.by = "freq")
dev.off()

## =========================
## 13. TF ACTIVITY — DoRothEA + VIPER (Figs 6–7)
## =========================
data(dorothea_hs, package = "dorothea")

regulon_df <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A","B","C"))

regulon <- split(regulon_df, regulon_df$tf)
regulon <- lapply(regulon, function(x) {
  tfmode     <- setNames(x$mor, x$target)
  likelihood <- rep(1, nrow(x))
  list(tfmode = tfmode, likelihood = likelihood)
})

tf_activity <- viper(logcpm_bc, regulon, minsize = 5, eset.filter = FALSE)

png(file.path(figdir, "Fig6_TF_activity.png"), width = 800, height = 1000)
pheatmap(tf_activity[1:50, ], scale = "row", annotation_col = ann_col,
         main = "Top 50 TF activities (VIPER)")
dev.off()

top_tf <- names(sort(apply(tf_activity, 1, var), decreasing = TRUE))[1:30]

png(file.path(figdir, "Fig7_TF_heatmap.png"), width = 800, height = 1000)
pheatmap(tf_activity[top_tf, ], scale = "row", annotation_col = ann_col,
         main = "Top 30 variable TFs (VIPER)")
dev.off()

tf_df <- as.data.frame(t(tf_activity))
tf_df$sample <- rownames(tf_df)
tf_df <- dplyr::left_join(tf_df, meta, by = "sample")

tf_summary <- tf_df %>%
  group_by(condition) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

write.csv(tf_summary, file.path(tabdir, "TF_activity_summary.csv"), row.names = FALSE)

## =========================
## 14. NETWORK CENTRALITY (Fig 8)
## =========================
g    <- graph_from_data_frame(regulon_df %>% dplyr::select(tf, target))
cent <- data.frame(gene = names(degree(g)), degree = degree(g))

p_net <- cent %>%
  arrange(desc(degree)) %>%
  head(20) %>%
  ggplot(aes(reorder(gene, degree), degree)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_classic() +
  labs(x = NULL, y = "Degree centrality",
       title = "Top 20 hub TFs (DoRothEA network)")

ggsave(file.path(figdir, "Fig8_Network.png"), p_net, width = 6, height = 5, dpi = 300)

## =========================
## 15. HELPER: SCORE GENE SET
## =========================
score_set <- function(mat, genes) {
  g <- base::intersect(genes, rownames(mat))
  if (length(g) < 3) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(mat[g, ])))
  z[is.na(z)] <- 0
  colMeans(z)
}

## =========================
## 16. STEMNESS SCORES (Fig 9)
## =========================
StemSets <- list(
  STEM_CORE      = c("CD34","PROM1","KIT","THY1","BMI1","MYC"),
  SELF_RENEWAL   = c("KLF4","SOX2","NANOG","MYC"),
  HEMATOPOIETIC  = c("MEIS1","HOXA9","PBX3","RUNX1"),
  QUIESCENT_STEM = c("FOXO3","CDKN1A","GADD45A","BTG1")
)

stem_scores <- as.data.frame(
  sapply(names(StemSets), function(nm) score_set(logcpm_bc, StemSets[[nm]]))
)
stem_scores$sample   <- colnames(logcpm_bc)
stem_scores          <- dplyr::left_join(stem_scores, meta, by = "sample")
stem_scores$Stemness <- rowMeans(stem_scores[, names(StemSets)], na.rm = TRUE)

write.csv(stem_scores, file.path(tabdir, "Stemness_scores.csv"), row.names = FALSE)

p_stem <- ggplot(stem_scores, aes(condition, Stemness, color = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_point(aes(shape = rep), size = 3,
             position = position_jitter(width = 0.08, seed = 1)) +
  theme_classic() +
  labs(title = "Composite stemness score", y = "Score", x = NULL)

ggsave(file.path(figdir, "Fig9_Stemness.png"), p_stem, width = 6, height = 4, dpi = 300)

## =========================
## 17. STEMNESS MARKER EXPRESSION vs VEH (Figs 10–13)
## =========================
all_stem_markers  <- unique(unlist(StemSets))
available_markers <- base::intersect(all_stem_markers, rownames(logcpm_bc))

message(sprintf("Stemness markers available: %d / %d",
  length(available_markers), length(all_stem_markers)))

# 17A. Heatmap
stem_mat_z <- t(scale(t(logcpm_bc[available_markers, ])))

row_ann <- data.frame(
  GeneSet = sapply(available_markers, function(g) {
    paste(names(StemSets)[sapply(StemSets, function(s) g %in% s)], collapse = "/")
  }),
  row.names = available_markers
)

png(file.path(figdir, "Fig10_StemMarker_Heatmap.png"), width = 900, height = 700)
pheatmap(stem_mat_z, annotation_col = ann_col, annotation_row = row_ann,
         fontsize_row = 9, main = "Stemness marker expression")
dev.off()

# 17B. Δ log-CPM dotplot vs Veh
stem_long <- logcpm_bc[available_markers, ] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "logCPM") %>%
  dplyr::left_join(meta, by = "sample")

veh_means <- stem_long %>%
  dplyr::filter(condition == "Veh") %>%
  group_by(gene) %>%
  summarise(veh_mean = mean(logCPM), .groups = "drop")

stem_fc <- stem_long %>%
  group_by(gene, condition) %>%
  summarise(mean_logCPM = mean(logCPM), .groups = "drop") %>%
  dplyr::left_join(veh_means, by = "gene") %>%
  mutate(delta_logCPM = mean_logCPM - veh_mean)

write.csv(stem_fc, file.path(tabdir, "StemMarker_FC_vs_Veh.csv"), row.names = FALSE)

p_dot <- stem_fc %>%
  dplyr::filter(condition != "Veh") %>%
  ggplot(aes(x = condition,
             y = reorder(gene, delta_logCPM),
             color = delta_logCPM,
             size  = abs(delta_logCPM))) +
  geom_point() +
  scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                        midpoint = 0, name = "Δ log-CPM\nvs Veh") +
  scale_size_continuous(range = c(2, 8), name = "|Δ log-CPM|") +
  theme_classic() +
  labs(x = NULL, y = NULL, title = "Stemness marker expression vs Vehicle") +
  theme(axis.text.y = element_text(size = 9))

ggsave(file.path(figdir, "Fig11_StemMarker_DotPlot.png"), p_dot,
       width = 7, height = 6, dpi = 300)

# 17C. Per gene-set barplot
set_membership <- stack(lapply(StemSets, function(g)
  base::intersect(g, rownames(logcpm_bc))))
colnames(set_membership) <- c("gene","GeneSet")

stem_fc_sets <- stem_fc %>%
  dplyr::filter(condition != "Veh") %>%
  dplyr::left_join(set_membership, by = "gene", relationship = "many-to-many") %>%
  group_by(GeneSet, condition) %>%
  summarise(mean_delta = mean(delta_logCPM, na.rm = TRUE), .groups = "drop")

p_bar_stem <- ggplot(stem_fc_sets, aes(condition, mean_delta, fill = mean_delta > 0)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("steelblue","firebrick")) +
  facet_wrap(~ GeneSet, scales = "free_y") +
  theme_classic() +
  labs(x = NULL, y = "Mean Δ log-CPM vs Veh",
       title = "Stemness gene set activity vs Vehicle")

ggsave(file.path(figdir, "Fig12_StemSet_Barplot.png"), p_bar_stem,
       width = 8, height = 6, dpi = 300)

# 17D. Per-gene boxplots, paginated
marker_chunks <- split(available_markers,
                       ceiling(seq_along(available_markers) / 9))

for (i in seq_along(marker_chunks)) {
  p_box <- stem_long %>%
    dplyr::filter(gene %in% marker_chunks[[i]]) %>%
    ggplot(aes(condition, logCPM, color = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    geom_point(aes(shape = rep), size = 2,
               position = position_jitter(width = 0.08, seed = 1)) +
    facet_wrap(~ gene, scales = "free_y", ncol = 3) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(x = NULL, y = "log-CPM")

  ggsave(file.path(figdir, sprintf("Fig13_%02d_StemMarker_Boxplots.png", i)),
         p_box, width = 9, height = 9, dpi = 300)
}

## =========================
## 18. CELL STATE AXIS (Fig 14)
## =========================
ProlifSet  <- c("MKI67","PCNA","E2F1","CCNB1","CDK1","MYC")
PersistSet <- c("CDKN1A","FOXO3","GADD45A","HMOX1","NQO1","BCL2")

axis_df <- data.frame(
  sample        = colnames(logcpm_bc),
  Proliferation = score_set(logcpm_bc, ProlifSet),
  Persistence   = score_set(logcpm_bc, PersistSet)
) %>%
  dplyr::left_join(meta, by = "sample") %>%
  mutate(State = Persistence - Proliferation)

write.csv(axis_df, file.path(tabdir, "Prolif_vs_Persist.csv"), row.names = FALSE)

p_axis <- ggplot(axis_df, aes(Proliferation, Persistence,
                               color = condition, label = sample)) +
  geom_point(size = 4) +
  geom_label_repel(size = 3, show.legend = FALSE) +
  theme_classic() +
  labs(title = "Cell state axis")

ggsave(file.path(figdir, "Fig14_CellState.png"), p_axis, width = 6, height = 5, dpi = 300)

## =========================
## 19. INTEGRATED ANALYSIS (Figs 15–16)
## =========================
integrated <- stem_scores %>%
  dplyr::select(sample, Stemness) %>%
  dplyr::left_join(
    axis_df %>% dplyr::select(sample, Proliferation, Persistence),
    by = "sample"
  ) %>%
  dplyr::left_join(meta, by = "sample")

cor_mat <- cor(integrated[, c("Stemness","Proliferation","Persistence")],
               use = "complete.obs")

png(file.path(figdir, "Fig15_ProgramCorrelation.png"), width = 600, height = 500)
pheatmap(cor_mat, main = "Program correlations",
         display_numbers = TRUE, number_format = "%.2f")
dev.off()

p_corr <- ggplot(integrated, aes(Stemness, Persistence, color = condition)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  theme_classic() +
  labs(title = "Stemness vs Persistence")

ggsave(file.path(figdir, "Fig16_Stem_vs_Persist.png"), p_corr, width = 6, height = 5, dpi = 300)

message("Pipeline complete.")
message("Figures : ", figdir)
message("Tables  : ", tabdir)

############################################################
# END OF PIPELINE
############################################################
