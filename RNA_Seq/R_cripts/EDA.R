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
plot_theme <- theme_bw(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

pca <- prcomp(t(vsd_matrix), scale. = TRUE)
explained_variance <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)

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
    title = "PCA",
    x = sprintf("PC1 (%.2f%%)", explained_variance[1] * 100),
    y = sprintf("PC2 (%.2f%%)", explained_variance[2] * 100),
    color = "Condition",
    shape = "Genotype"
  ) +
  plot_theme
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

annotation_df <- metadata[, c("genotype", "condition"), drop = FALSE]
colnames(annotation_df) <- c("Genotype", "Condition")

pheatmap(
  euclidean_dist,
  clustering_distance_rows = as.dist(euclidean_dist),
  clustering_distance_cols = as.dist(euclidean_dist),
  clustering_method = "ward.D2",
  annotation_col = annotation_df,
  annotation_row = annotation_df,
  main = "Евклидово расстояние между образцами",
  fontsize = 10,
  filename = file.path(output_dir, "heatmap_euclidean.png"),
  width = 9,
  height = 8
)

hc <- hclust(as.dist(euclidean_dist), method = "ward.D2")
png(filename = file.path(output_dir, "hierarchical_clustering.png"), width = 1300, height = 900, res = 150)
plot(hc, main = "Иерархическая кластеризация образцов", xlab = "", sub = "", cex = 0.9, font.main = 2)
dev.off()

cor_matrix <- cor(vsd_matrix, method = "pearson")
cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("sample_1", "sample_2", "correlation")
cor_plot <- ggplot(cor_df, aes(x = sample_1, y = sample_2, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "#00C2FF", mid = "#FF6B6B", high = "#FF1744", midpoint = 0.9, limits = c(0, 1)) +
  labs(
    title = "Корреляция между репликатами",
    x = "",
    y = "",
    fill = "Корреляция"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9)
  )
ggsave(filename = file.path(output_dir, "replicate_correlation_heatmap.png"), plot = cor_plot, width = 10, height = 8, dpi = 300)

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

outlier_plot <- ggplot(outlier_table, aes(x = median_within_group_correlation, y = z_distance, color = exclude_recommended, label = sample)) +
  geom_vline(xintercept = 0.90, linetype = "dashed", color = "#636363") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "#636363") +
  geom_point(size = 3.3) +
  geom_text(vjust = -0.8, size = 3.1) +
  scale_color_manual(values = c(`TRUE` = "#D7301F", `FALSE` = "#225EA8")) +
  labs(
    title = "Оценка выбросов",
    x = "Медианная внутригрупповая корреляция",
    y = "Z-оценка расстояния до центроида группы",
    color = "Рекомендовано исключить"
  ) +
  plot_theme
ggsave(filename = file.path(output_dir, "outlier_assessment.png"), plot = outlier_plot, width = 10, height = 7, dpi = 300)

writeLines("EDA завершен. Результаты сохранены в папке EDA_results.")
