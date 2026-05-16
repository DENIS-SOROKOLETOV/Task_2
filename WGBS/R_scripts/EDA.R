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

context_summary_list <- list()
chr_profile_list <- list()

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
    if (context == "CpG" && nrow(stats$chr_table) > 0) {
      chr_df <- stats$chr_table
      chr_df$sample_id <- sample_id
      chr_profile_list[[length(chr_profile_list) + 1]] <- chr_df
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

context_summary$context <- factor(context_summary$context, levels = contexts)
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
  geom_col(position = "stack", width = 0.75) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
  labs(
    title = "Доля покрытия по контекстам (CpG / CHG / CHH)",
    x = "Образец",
    y = "Доля сайтов с покрытием",
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

if (length(chr_profile_list) > 0) {
  chr_profiles <- do.call(rbind, chr_profile_list)
  write.table(
    chr_profiles,
    file.path(output_dir, "cpg_methylation_by_chromosome.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  chr_matrix <- reshape(
    chr_profiles[, c("sample_id", "chr", "weighted_meth_pct")],
    idvar = "chr",
    timevar = "sample_id",
    direction = "wide"
  )
  rownames(chr_matrix) <- chr_matrix$chr
  chr_matrix$chr <- NULL
  chr_matrix_num <- as.matrix(chr_matrix)
  mode(chr_matrix_num) <- "numeric"
  chr_matrix_num[is.na(chr_matrix_num)] <- 0
  row_var <- apply(chr_matrix_num, 1, var, na.rm = TRUE)
  chr_matrix_num <- chr_matrix_num[row_var > 0 & !is.na(row_var), , drop = FALSE]
  if (ncol(chr_matrix_num) >= 2 && nrow(chr_matrix_num) >= 2) {
    cor_chr <- cor(chr_matrix_num, method = "pearson", use = "pairwise.complete.obs")
    write.table(
      cor_chr,
      file.path(output_dir, "sample_correlation_chromosome_cpg.tsv"),
      sep = "\t",
      quote = FALSE
    )
    cor_df <- as.data.frame(as.table(cor_chr))
    colnames(cor_df) <- c("sample_1", "sample_2", "correlation")
    cor_plot <- ggplot(cor_df, aes(x = sample_1, y = sample_2, fill = correlation)) +
      geom_tile() +
      scale_fill_gradient2(low = "#00C2FF", mid = "white", high = "#FF1744", midpoint = 0.5, limits = c(-1, 1)) +
      labs(title = "Корреляция образцов (профиль CpG по хромосомам)", x = "", y = "", fill = "r") +
      plot_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(output_dir, "sample_correlation_heatmap.png"), cor_plot, width = 8, height = 7, dpi = 300)
    pheatmap(
      cor_chr,
      clustering_distance_rows = as.dist(1 - cor_chr),
      clustering_distance_cols = as.dist(1 - cor_chr),
      clustering_method = "ward.D2",
      main = "Сходство образцов (CpG по хромосомам)",
      display_numbers = TRUE,
      number_format = "%.2f",
      fontsize = 10,
      filename = file.path(output_dir, "sample_similarity_heatmap.png"),
      width = 8,
      height = 7
    )
    pca_chr <- prcomp(t(chr_matrix_num), scale. = TRUE)
    explained <- (pca_chr$sdev^2) / sum(pca_chr$sdev^2)
    pca_df <- data.frame(
      sample_id = rownames(pca_chr$x),
      PC1 = pca_chr$x[, 1],
      PC2 = pca_chr$x[, 2],
      stringsAsFactors = FALSE
    )
    pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = sample_id)) +
      geom_point(size = 3.5, color = "#225EA8") +
      geom_text(vjust = -0.8, size = 3.2) +
      labs(
        title = "PCA: профиль CpG-метилирования по хромосомам",
        x = sprintf("PC1 (%.1f%%)", explained[1] * 100),
        y = sprintf("PC2 (%.1f%%)", explained[2] * 100)
      ) +
      plot_theme
    ggsave(file.path(output_dir, "pca_chromosome_profile.png"), pca_plot, width = 8, height = 6, dpi = 300)
    write.table(
      pca_df,
      file.path(output_dir, "pca_coordinates.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    dist_chr <- dist(t(chr_matrix_num), method = "euclidean")
    hc <- hclust(dist_chr, method = "ward.D2")
    png(file.path(output_dir, "hierarchical_clustering.png"), width = 1200, height = 800, res = 150)
    plot(hc, main = "Иерархическая кластеризация образцов", xlab = "", sub = "")
    dev.off()
  }
}

writeLines("EDA завершен. Результаты сохранены в папке WGBS/EDA_results.")
