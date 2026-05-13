required_packages <- c("DESeq2", "ggplot2", "pheatmap")
missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  stop(sprintf("Не установлены пакеты: %s", paste(missing_required, collapse = ", ")))
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
  script_dir <- dirname(script_path)
} else {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
counts_root <- file.path(project_root, "Quantification_results")
output_dir <- file.path(project_root, "EDA_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

count_files <- list.files(counts_root, pattern = "counts\\.txt$", recursive = TRUE, full.names = TRUE)
if (length(count_files) == 0) {
  stop("Файлы counts.txt не найдены в Quantification_results.")
}

sample_ids <- basename(dirname(count_files))
names(count_files) <- sample_ids

read_count_file <- function(path, sample_id) {
  tab <- read.delim(path, sep = "\t", header = TRUE, comment.char = "#", check.names = FALSE)
  if (!all(c("Geneid") %in% colnames(tab))) {
    stop(sprintf("В файле %s отсутствует колонка Geneid.", path))
  }
  count_col <- colnames(tab)[ncol(tab)]
  out <- data.frame(Geneid = tab$Geneid, count = as.numeric(tab[[count_col]]), stringsAsFactors = FALSE)
  colnames(out)[2] <- sample_id
  out
}

count_list <- mapply(read_count_file, count_files, names(count_files), SIMPLIFY = FALSE)
count_table <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_list)
count_table[is.na(count_table)] <- 0
rownames(count_table) <- count_table$Geneid
count_matrix <- as.matrix(count_table[, -1, drop = FALSE])
storage.mode(count_matrix) <- "integer"

metadata <- data.frame(sample_id = colnames(count_matrix), stringsAsFactors = FALSE)
metadata$genotype <- NA_character_
metadata$condition <- NA_character_

metadata$genotype[metadata$sample_id %in% c("SRR10408671", "SRR10408672", "SRR10408673", "SRR10408683", "SRR10408684", "SRR10408685")] <- "WT"
metadata$genotype[metadata$sample_id %in% c("SRR10408674", "SRR10408675", "SRR10408676", "SRR10408686", "SRR10408687", "SRR10408688")] <- "DoubleMutant"
metadata$condition[metadata$sample_id %in% c("SRR10408671", "SRR10408672", "SRR10408673", "SRR10408674", "SRR10408675", "SRR10408676")] <- "PN"
metadata$condition[metadata$sample_id %in% c("SRR10408683", "SRR10408684", "SRR10408685", "SRR10408686", "SRR10408687", "SRR10408688")] <- "NPN"

if (any(is.na(metadata$genotype)) || any(is.na(metadata$condition))) {
  stop("Для части образцов не удалось определить genotype/condition по Variant 4.")
}

metadata$genotype <- factor(metadata$genotype, levels = c("WT", "DoubleMutant"))
metadata$condition <- factor(metadata$condition, levels = c("PN", "NPN"))
metadata$group <- factor(paste(metadata$genotype, metadata$condition, sep = "_"))
rownames(metadata) <- metadata$sample_id
metadata <- metadata[colnames(count_matrix), , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~ genotype + condition + genotype:condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
vsd_matrix <- assay(vsd)

pca <- prcomp(t(vsd_matrix), scale. = TRUE)
explained_variance <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
variance_table <- data.frame(
  PC = paste0("PC", seq_along(explained_variance)),
  explained_variance = explained_variance,
  explained_variance_percent = explained_variance * 100,
  stringsAsFactors = FALSE
)
write.table(variance_table, file.path(output_dir, "pca_explained_variance.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)
pca_df <- merge(pca_df, metadata, by.x = "sample_id", by.y = "sample_id", all.x = TRUE, sort = FALSE)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = genotype, label = sample_id)) +
  geom_point(size = 3.3) +
  geom_text(vjust = -0.8, size = 3.1) +
  labs(
    title = "PCA: проверка batch effects",
    x = sprintf("PC1 (%.2f%%)", explained_variance[1] * 100),
    y = sprintf("PC2 (%.2f%%)", explained_variance[2] * 100),
    color = "Condition",
    shape = "Genotype"
  ) +
  theme_bw(base_size = 12)
ggsave(filename = file.path(output_dir, "PCA_batch_check.png"), plot = pca_plot, width = 9, height = 7, dpi = 300)

pc1_lm <- lm(PC1 ~ genotype + condition + genotype:condition, data = pca_df)
pc2_lm <- lm(PC2 ~ genotype + condition + genotype:condition, data = pca_df)
batch_effect_report <- data.frame(
  component = c("PC1", "PC2"),
  r_squared = c(summary(pc1_lm)$r.squared, summary(pc2_lm)$r.squared),
  p_condition = c(anova(pc1_lm)["condition", "Pr(>F)"], anova(pc2_lm)["condition", "Pr(>F)"]),
  p_genotype = c(anova(pc1_lm)["genotype", "Pr(>F)"], anova(pc2_lm)["genotype", "Pr(>F)"]),
  p_interaction = c(anova(pc1_lm)["genotype:condition", "Pr(>F)"], anova(pc2_lm)["genotype:condition", "Pr(>F)"])
)
write.table(batch_effect_report, file.path(output_dir, "batch_effect_model.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

euclidean_dist <- as.matrix(dist(t(vsd_matrix), method = "euclidean"))
write.table(euclidean_dist, file.path(output_dir, "euclidean_distance_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)

annotation_df <- metadata[, c("genotype", "condition"), drop = FALSE]
colnames(annotation_df) <- c("Genotype", "Condition")

pheatmap(
  euclidean_dist,
  clustering_distance_rows = as.dist(euclidean_dist),
  clustering_distance_cols = as.dist(euclidean_dist),
  clustering_method = "ward.D2",
  annotation_col = annotation_df,
  annotation_row = annotation_df,
  fontsize = 10,
  filename = file.path(output_dir, "heatmap_euclidean.png"),
  width = 9,
  height = 8
)

poisson_status <- data.frame(method = "PoissonDistance", available = FALSE, note = "Пакет PoiClaClu не установлен", stringsAsFactors = FALSE)
if (requireNamespace("PoiClaClu", quietly = TRUE)) {
  poisson_status <- tryCatch({
    poisson_dist <- as.matrix(PoiClaClu::PoissonDistance(t(counts(dds, normalized = TRUE)))$dd)
    expected_n <- ncol(dds)
    if (nrow(poisson_dist) != expected_n || ncol(poisson_dist) != expected_n) {
      poisson_dist <- as.matrix(PoiClaClu::PoissonDistance(counts(dds, normalized = TRUE))$dd)
    }
    sample_order <- colnames(dds)
    rownames(poisson_dist) <- sample_order
    colnames(poisson_dist) <- sample_order
    write.table(poisson_dist, file.path(output_dir, "poisson_distance_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
    pheatmap(
      poisson_dist,
      clustering_distance_rows = as.dist(poisson_dist),
      clustering_distance_cols = as.dist(poisson_dist),
      clustering_method = "ward.D2",
      annotation_col = annotation_df,
      annotation_row = annotation_df,
      fontsize = 10,
      filename = file.path(output_dir, "heatmap_poisson.png"),
      width = 9,
      height = 8
    )
    data.frame(method = "PoissonDistance", available = TRUE, note = "Poisson heatmap построена", stringsAsFactors = FALSE)
  }, error = function(e) {
    data.frame(method = "PoissonDistance", available = FALSE, note = paste("Ошибка построения:", conditionMessage(e)), stringsAsFactors = FALSE)
  })
}
write.table(poisson_status, file.path(output_dir, "poisson_distance_status.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

hc <- hclust(as.dist(euclidean_dist), method = "ward.D2")
png(filename = file.path(output_dir, "hierarchical_clustering.png"), width = 1300, height = 900, res = 150)
plot(hc, main = "Hierarchical clustering of samples", xlab = "", sub = "", cex = 0.9)
dev.off()

clusters_k2 <- cutree(hc, k = 2)
cluster_table <- data.frame(
  sample = names(clusters_k2),
  cluster_k2 = as.integer(clusters_k2),
  condition = metadata[names(clusters_k2), "condition"],
  genotype = metadata[names(clusters_k2), "genotype"],
  stringsAsFactors = FALSE
)
write.table(cluster_table, file.path(output_dir, "hierarchical_clusters_k2.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cor_matrix <- cor(vsd_matrix, method = "pearson")
write.table(cor_matrix, file.path(output_dir, "replicate_correlation_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)

pairwise_correlation <- list()
for (grp in levels(metadata$group)) {
  grp_samples <- rownames(metadata)[metadata$group == grp]
  if (length(grp_samples) < 2) {
    next
  }
  cmb <- combn(grp_samples, 2)
  grp_df <- data.frame(
    group = grp,
    sample_1 = cmb[1, ],
    sample_2 = cmb[2, ],
    correlation = mapply(function(a, b) cor_matrix[a, b], cmb[1, ], cmb[2, ]),
    stringsAsFactors = FALSE
  )
  grp_df$pass_gt_0_90 <- grp_df$correlation > 0.90
  pairwise_correlation[[grp]] <- grp_df
}
pairwise_correlation_df <- do.call(rbind, pairwise_correlation)
write.table(pairwise_correlation_df, file.path(output_dir, "replicate_correlations_by_group.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

median_within_group_corr <- sapply(rownames(metadata), function(sample_name) {
  grp <- metadata[sample_name, "group"]
  grp_samples <- rownames(metadata)[metadata$group == grp]
  peer_samples <- setdiff(grp_samples, sample_name)
  if (length(peer_samples) == 0) {
    return(NA_real_)
  }
  median(cor_matrix[sample_name, peer_samples], na.rm = TRUE)
})

pc_coords <- pca_df[, c("sample_id", "PC1", "PC2", "group")]
group_centroid <- aggregate(cbind(PC1, PC2) ~ group, data = pc_coords, FUN = mean)
pc_coords <- merge(pc_coords, group_centroid, by = "group", suffixes = c("", "_centroid"), all.x = TRUE, sort = FALSE)
pc_coords$distance_to_group_centroid <- sqrt((pc_coords$PC1 - pc_coords$PC1_centroid) ^ 2 + (pc_coords$PC2 - pc_coords$PC2_centroid) ^ 2)

pc_coords$z_distance <- NA_real_
for (grp in unique(pc_coords$group)) {
  idx <- which(pc_coords$group == grp)
  grp_dist <- pc_coords$distance_to_group_centroid[idx]
  if (length(grp_dist) >= 2 && sd(grp_dist) > 0) {
    pc_coords$z_distance[idx] <- (grp_dist - mean(grp_dist)) / sd(grp_dist)
  } else {
    pc_coords$z_distance[idx] <- 0
  }
}

outlier_table <- data.frame(
  sample = pc_coords$sample_id,
  group = pc_coords$group,
  median_within_group_correlation = median_within_group_corr[pc_coords$sample_id],
  distance_to_group_centroid = pc_coords$distance_to_group_centroid,
  z_distance = pc_coords$z_distance,
  stringsAsFactors = FALSE
)
outlier_table$flag_corr_lt_0_90 <- outlier_table$median_within_group_correlation < 0.90
outlier_table$flag_z_distance_gt_3 <- outlier_table$z_distance > 3
outlier_table$exclude_recommended <- outlier_table$flag_corr_lt_0_90 & outlier_table$flag_z_distance_gt_3

outlier_table$formal_justification <- ifelse(
  outlier_table$exclude_recommended,
  "Исключение рекомендовано: низкая внутригрупповая корреляция (<0.90) и сильное отклонение в PCA (z>3).",
  ifelse(
    outlier_table$flag_corr_lt_0_90 | outlier_table$flag_z_distance_gt_3,
    "Пограничный случай: требуется ручная проверка FASTQC, картирования и биологического контекста.",
    "Формальных оснований для исключения образца нет."
  )
)
write.table(outlier_table, file.path(output_dir, "outlier_analysis.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

summary_report <- c(
  sprintf("Количество образцов: %d", ncol(count_matrix)),
  sprintf("Количество генов после фильтрации: %d", nrow(vsd_matrix)),
  sprintf("PC1 объясняет %.2f%% вариации", explained_variance[1] * 100),
  sprintf("PC2 объясняет %.2f%% вариации", explained_variance[2] * 100),
  sprintf("Доля пар репликатов с корреляцией >0.90: %.2f%%", 100 * mean(pairwise_correlation_df$pass_gt_0_90)),
  sprintf("Количество образцов с рекомендацией исключения: %d", sum(outlier_table$exclude_recommended))
)
writeLines(summary_report, con = file.path(output_dir, "EDA_report.txt"))

writeLines("EDA завершен. Результаты сохранены в папке EDA_results.")
