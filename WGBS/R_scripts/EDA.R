required_packages <- c("ggplot2", "pheatmap")
missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  stop(sprintf("Не установлены пакеты: %s", paste(missing_required, collapse = ", ")))
}

suppressPackageStartupMessages({
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

wgbs_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
extract_root <- file.path(wgbs_root, "Extract_results")
output_dir <- file.path(wgbs_root, "EDA_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

contexts <- c("CpG", "CHG", "CHH")
min_depth <- 5L

find_bedgraph <- function(sample_id, context) {
  sample_dir <- file.path(extract_root, sample_id)
  candidates <- c(
    file.path(sample_dir, paste0(sample_id, "_", context, ".bedGraph")),
    file.path(sample_dir, paste0(sample_id, "_", context, ".bedGraph.gz"))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[1]
}

summarize_bedgraph <- function(path, min_depth = 5L) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "r") else file(path, "r")
  on.exit(close(con), add = TRUE)
  n_sites <- 0L
  n_methyl <- 0
  n_unmethyl <- 0
  chr_methyl <- list()
  chr_unmethyl <- list()
  while (length(line <- readLines(con, n = 1L)) > 0) {
    if (startsWith(line, "track")) {
      next
    }
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 6) {
      next
    }
    chr_name <- parts[1]
    methyl <- as.numeric(parts[5])
    unmethyl <- as.numeric(parts[6])
    if (is.na(methyl) || is.na(unmethyl)) {
      next
    }
    depth <- methyl + unmethyl
    if (depth < min_depth) {
      next
    }
    n_sites <- n_sites + 1L
    n_methyl <- n_methyl + methyl
    n_unmethyl <- n_unmethyl + unmethyl
    if (is.null(chr_methyl[[chr_name]])) {
      chr_methyl[[chr_name]] <- 0
      chr_unmethyl[[chr_name]] <- 0
    }
    chr_methyl[[chr_name]] <- chr_methyl[[chr_name]] + methyl
    chr_unmethyl[[chr_name]] <- chr_unmethyl[[chr_name]] + unmethyl
  }
  total_depth <- n_methyl + n_unmethyl
  weighted_meth <- if (total_depth > 0) 100 * n_methyl / total_depth else NA_real_
  chr_table <- data.frame(
    chr = names(chr_methyl),
    n_methyl = unlist(chr_methyl, use.names = FALSE),
    n_unmethyl = unlist(chr_unmethyl, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  chr_table$depth <- chr_table$n_methyl + chr_table$n_unmethyl
  chr_table$weighted_meth_pct <- ifelse(chr_table$depth > 0, 100 * chr_table$n_methyl / chr_table$depth, NA_real_)
  list(
    n_sites = n_sites,
    n_methyl = n_methyl,
    n_unmethyl = n_unmethyl,
    weighted_meth_pct = weighted_meth,
    chr_table = chr_table
  )
}

sample_dirs <- list.dirs(extract_root, recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[grepl("^SRR", sample_dirs)]
if (length(sample_dirs) == 0) {
  stop("В Extract_results не найдены папки образцов (SRR*). Сначала запустите Extract_script.sh.")
}
sample_dirs <- sort(sample_dirs)

sample_annotation <- data.frame(
  sample_id = c("SRR10382875", "SRR10382876", "SRR10382877", "SRR10382878"),
  genotype = c("DeltaDelta", "DeltaDelta", "WT", "WT"),
  condition = c("NPN", "PN", "NPN", "PN"),
  stringsAsFactors = FALSE
)

present_samples <- intersect(sample_annotation$sample_id, sample_dirs)
if (length(present_samples) < 2) {
  stop("Для EDA нужно минимум 2 образца из набора SRR10382875-78.")
}
sample_dirs <- sample_annotation$sample_id[sample_annotation$sample_id %in% present_samples]

context_summary_list <- list()
feature_profile_list <- list()

for (sample_id in sample_dirs) {
  for (context in contexts) {
    bed_path <- find_bedgraph(sample_id, context)
    if (is.na(bed_path)) {
      warning(sprintf("Пропуск %s (%s): bedGraph не найден.", sample_id, context))
      next
    }
    stats <- summarize_bedgraph(bed_path, min_depth = min_depth)
    context_summary_list[[length(context_summary_list) + 1]] <- data.frame(
      sample_id = sample_id,
      context = context,
      n_sites = stats$n_sites,
      n_methyl = stats$n_methyl,
      n_unmethyl = stats$n_unmethyl,
      depth = stats$n_methyl + stats$n_unmethyl,
      weighted_meth_pct = stats$weighted_meth_pct,
      stringsAsFactors = FALSE
    )
    if (nrow(stats$chr_table) > 0) {
      chr_df <- stats$chr_table
      chr_df$sample_id <- sample_id
      chr_df$context <- context
      chr_df$feature_id <- paste(chr_df$context, chr_df$chr, sep = "|")
      feature_profile_list[[length(feature_profile_list) + 1]] <- chr_df
    }
  }
}

if (length(context_summary_list) == 0) {
  stop("Не удалось прочитать ни одного bedGraph. Проверьте Extract_results.")
}

context_summary <- do.call(rbind, context_summary_list)
write.table(
  context_summary,
  file.path(output_dir, "context_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

global_by_sample <- aggregate(
  cbind(n_methyl, n_unmethyl, depth) ~ sample_id,
  data = context_summary,
  FUN = sum
)
global_by_sample$weighted_meth_pct <- ifelse(
  global_by_sample$depth > 0,
  100 * global_by_sample$n_methyl / global_by_sample$depth,
  NA_real_
)
write.table(
  global_by_sample,
  file.path(output_dir, "global_methylation_by_sample.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

context_summary$context <- factor(context_summary$context, levels = c("CHH", "CHG", "CpG"))
plot_theme <- theme_bw(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

meth_by_context_plot <- ggplot(context_summary, aes(x = context, y = weighted_meth_pct, fill = context)) +
  geom_col(width = 0.7) +
  facet_wrap(~ sample_id, scales = "free_y") +
  labs(
    title = "Средний уровень метилирования по контекстам",
    x = "Контекст",
    y = "Метилирование, % (взвешенное)"
  ) +
  plot_theme +
  theme(legend.position = "none")
ggsave(
  file.path(output_dir, "methylation_by_context.png"),
  meth_by_context_plot,
  width = 10,
  height = 7,
  dpi = 300
)

context_summary$context_fraction <- ave(
  context_summary$depth,
  context_summary$sample_id,
  FUN = function(x) if (sum(x) > 0) x / sum(x) else NA_real_
)

context_fraction_plot <- ggplot(context_summary, aes(x = sample_id, y = context_fraction, fill = context)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
  labs(
    title = "Покрытие по контекстам",
    x = "Образец",
    y = "Доля сайтов с покрытием в контексте",
    fill = "Контекст"
  ) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(output_dir, "context_coverage_fraction.png"),
  context_fraction_plot,
  width = 9,
  height = 6,
  dpi = 300
)

if (length(feature_profile_list) > 0) {
  feature_profiles <- do.call(rbind, feature_profile_list)
  write.table(
    feature_profiles,
    file.path(output_dir, "methylation_by_context_chromosome.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  feature_matrix <- reshape(
    feature_profiles[, c("sample_id", "feature_id", "weighted_meth_pct")],
    idvar = "feature_id",
    timevar = "sample_id",
    direction = "wide"
  )
  rownames(feature_matrix) <- feature_matrix$feature_id
  feature_matrix$feature_id <- NULL
  feature_matrix_num <- as.matrix(feature_matrix)
  mode(feature_matrix_num) <- "numeric"
  feature_matrix_num[is.na(feature_matrix_num)] <- 0
  row_var <- apply(feature_matrix_num, 1, var, na.rm = TRUE)
  feature_matrix_num <- feature_matrix_num[row_var > 0 & !is.na(row_var), , drop = FALSE]

  if (ncol(feature_matrix_num) >= 2 && nrow(feature_matrix_num) >= 2) {
    cor_samples <- cor(feature_matrix_num, method = "pearson", use = "pairwise.complete.obs")
    write.table(
      cor_samples,
      file.path(output_dir, "sample_correlation_context_chr.tsv"),
      sep = "\t",
      quote = FALSE
    )

    sample_ids <- sub("^weighted_meth_pct\\.", "", colnames(feature_matrix_num))
    colnames(feature_matrix_num) <- sample_ids
    sample_meta <- merge(
      data.frame(sample_id = sample_ids, stringsAsFactors = FALSE),
      sample_annotation,
      by = "sample_id",
      all.x = TRUE,
      sort = FALSE
    )
    if (any(is.na(sample_meta$genotype)) || any(is.na(sample_meta$condition))) {
      stop("Для части WGBS-образцов не удалось определить genotype/condition.")
    }
    sample_meta$genotype <- factor(sample_meta$genotype, levels = c("WT", "DeltaDelta"))
    sample_meta$condition <- factor(sample_meta$condition, levels = c("PN", "NPN"))
    rownames(sample_meta) <- sample_meta$sample_id

    cor_df <- as.data.frame(as.table(cor_samples))
    colnames(cor_df) <- c("sample_1", "sample_2", "correlation")
    cor_plot <- ggplot(cor_df, aes(x = sample_1, y = sample_2, fill = correlation)) +
      geom_tile() +
      scale_fill_gradient2(low = "#00C2FF", mid = "#FF6B6B", high = "#FF1744", midpoint = 0.9, limits = c(0, 1)) +
      labs(title = "Корреляция между репликатами", x = "", y = "", fill = "Корреляция") +
      plot_theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9)
      )
    ggsave(file.path(output_dir, "replicate_correlation_heatmap.png"), cor_plot, width = 10, height = 8, dpi = 300)

    euclidean_dist <- as.matrix(dist(t(feature_matrix_num), method = "euclidean"))
    pheatmap(
      euclidean_dist,
      clustering_distance_rows = as.dist(euclidean_dist),
      clustering_distance_cols = as.dist(euclidean_dist),
      clustering_method = "ward.D2",
      annotation_col = sample_meta[, c("genotype", "condition"), drop = FALSE],
      annotation_row = sample_meta[, c("genotype", "condition"), drop = FALSE],
      main = "Евклидово расстояние между WGBS-образцами",
      fontsize = 10,
      filename = file.path(output_dir, "heatmap_euclidean.png"),
      width = 9,
      height = 8
    )

    pca_chr <- prcomp(t(feature_matrix_num), scale. = TRUE)
    explained <- (pca_chr$sdev ^ 2) / sum(pca_chr$sdev ^ 2)
    pca_df <- data.frame(
      sample_id = rownames(pca_chr$x),
      PC1 = pca_chr$x[, 1],
      PC2 = pca_chr$x[, 2],
      stringsAsFactors = FALSE
    )
    pca_df <- merge(pca_df, sample_meta, by = "sample_id", all.x = TRUE, sort = FALSE)

    pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = genotype, label = sample_id)) +
      geom_point(size = 3.3) +
      geom_text(vjust = -0.8, size = 3.1) +
      labs(
        title = "PCA",
        x = sprintf("PC1 (%.2f%%)", explained[1] * 100),
        y = sprintf("PC2 (%.2f%%)", explained[2] * 100),
        color = "Condition",
        shape = "Genotype"
      ) +
      plot_theme
    ggsave(file.path(output_dir, "PCA_batch_check.png"), pca_plot, width = 9, height = 7, dpi = 300)

    write.table(
      pca_df,
      file.path(output_dir, "pca_coordinates.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    dist_chr <- dist(t(feature_matrix_num), method = "euclidean")
    hc <- hclust(dist_chr, method = "ward.D2")
    png(file.path(output_dir, "hierarchical_clustering.png"), width = 1300, height = 900, res = 150)
    plot(hc, main = "Иерархическая кластеризация образцов", xlab = "", sub = "", cex = 0.9, font.main = 2)
    dev.off()
  }
}

writeLines("EDA завершен. Результаты сохранены в папке WGBS/EDA_results.")
