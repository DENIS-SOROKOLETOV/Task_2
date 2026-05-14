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
output_dir <- file.path(project_root, "DE_results")
gtf_path <- file.path(project_root, "Reference_genome", "Fusarium_graminearum_ph_1_gca_000240135.ASM24013v3.62.gtf")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(gtf_path)) {
  stop("Не найден референсный GTF: Fusarium_graminearum_ph_1_gca_000240135.ASM24013v3.62.gtf")
}

count_files <- list.files(counts_root, pattern = "counts\\.txt$", recursive = TRUE, full.names = TRUE)
if (length(count_files) == 0) {
  stop("Файлы counts.txt не найдены в Quantification_results.")
}

sample_ids <- basename(dirname(count_files))
names(count_files) <- sample_ids

read_count_file <- function(path, sample_id) {
  tab <- read.delim(path, sep = "\t", header = TRUE, comment.char = "#", check.names = FALSE)
  if (!("Geneid" %in% colnames(tab))) {
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

extract_gtf_attr <- function(attr_string, key) {
  pattern <- paste0(key, " \"([^\"]+)\"")
  m <- regexec(pattern, attr_string)
  res <- regmatches(attr_string, m)
  out <- vapply(res, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  out
}

gtf <- read.delim(
  gtf_path,
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "#",
  stringsAsFactors = FALSE
)
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
gtf_gene <- gtf[gtf$feature == "gene", , drop = FALSE]
if (nrow(gtf_gene) == 0) {
  gtf_gene <- gtf[gtf$feature == "transcript", , drop = FALSE]
}
gene_id_raw <- extract_gtf_attr(gtf_gene$attribute, "gene_id")
gene_name <- extract_gtf_attr(gtf_gene$attribute, "gene_name")
gene_biotype <- extract_gtf_attr(gtf_gene$attribute, "gene_biotype")
gene_name[is.na(gene_name) | gene_name == ""] <- gene_id_raw[is.na(gene_name) | gene_name == ""]
gene_annotation <- data.frame(
  Geneid = gene_id_raw,
  gene_id = gene_id_raw,
  gene_name = gene_name,
  gene_biotype = gene_biotype,
  stringsAsFactors = FALSE
)
gene_annotation <- gene_annotation[!is.na(gene_annotation$Geneid), , drop = FALSE]
gene_annotation <- gene_annotation[!duplicated(gene_annotation$Geneid), , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~ genotype + condition + genotype:condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)
coef_names <- resultsNames(dds)

interaction_coef <- grep("genotype.*condition|condition.*genotype", coef_names, value = TRUE)
if (length(interaction_coef) == 0) {
  stop("Не найден interaction-коэффициент в resultsNames(dds).")
}
interaction_coef <- interaction_coef[1]

contrast_specs <- list(
  list(
    id = "C1_WT_NPN_vs_WT_PN",
    title = "Эффект среды в WT: NPN vs PN",
    plot_title = "WT: NPN × PN",
    type = "contrast",
    value = c("condition", "NPN", "PN")
  ),
  list(
    id = "C2_DoubleMutant_NPN_vs_DoubleMutant_PN",
    title = "Эффект среды в ΔΔ: NPN vs PN",
    plot_title = "ΔΔ: NPN × PN",
    type = "contrast_list",
    value = list(c("condition_NPN_vs_PN", interaction_coef))
  ),
  list(
    id = "C3_DoubleMutant_PN_vs_WT_PN",
    title = "Эффект генотипа в PN: ΔΔ vs WT",
    plot_title = "PN: ΔΔ × WT",
    type = "contrast",
    value = c("genotype", "DoubleMutant", "WT")
  ),
  list(
    id = "C4_DoubleMutant_NPN_vs_WT_NPN",
    title = "Эффект генотипа в NPN: ΔΔ vs WT",
    plot_title = "NPN: ΔΔ × WT",
    type = "contrast_list",
    value = list(c("genotype_DoubleMutant_vs_WT", interaction_coef))
  ),
  list(
    id = "C5_Interaction_effect",
    title = "Интеракция genotype:condition",
    plot_title = "[WT/ΔΔ×PN/NPN]",
    type = "name",
    value = interaction_coef
  )
)

get_res <- function(dds_obj, spec) {
  if (spec$type == "contrast") {
    return(results(dds_obj, contrast = spec$value, alpha = 0.05))
  }
  if (spec$type == "contrast_list") {
    return(results(dds_obj, contrast = spec$value, alpha = 0.05))
  }
  if (spec$type == "name") {
    return(results(dds_obj, name = spec$value, alpha = 0.05))
  }
  stop(sprintf("Неизвестный тип контраста: %s", spec$type))
}

make_volcano <- function(df, title_text, output_path) {
  df$neg_log10_padj <- -log10(df$padj)
  df$neg_log10_padj[is.infinite(df$neg_log10_padj)] <- NA_real_
  p <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj, color = DEG_status)) +
    geom_point(alpha = 0.6, size = 1.2, na.rm = TRUE) +
    scale_color_manual(values = c(up = "#D7301F", down = "#225EA8", ns = "#808080", not_tested = "#BDBDBD")) +
    labs(
      title = title_text,
      x = "log2FoldChange",
      y = "-log10(adjusted p-value)",
      color = "DE status"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  ggsave(filename = output_path, plot = p, width = 8.5, height = 6.5, dpi = 300)
}

annotate_res_table <- function(df, annot_tbl) {
  merged <- merge(df, annot_tbl, by = "Geneid", all.x = TRUE, sort = FALSE)
  merged
}

get_contrast_folder <- function(contrast_id) {
  sub("^C[0-9]+_", "", contrast_id)
}

make_contrast_paths <- function(base_dir, spec) {
  contrast_dir <- file.path(base_dir, get_contrast_folder(spec$id))
  ma_dir <- file.path(contrast_dir, "MA_plot")
  volcano_dir <- file.path(contrast_dir, "Volcano_plot")
  heatmap_dir <- file.path(contrast_dir, "Heatmap")
  dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ma_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)
  list(
    contrast_dir = contrast_dir,
    all_genes = file.path(contrast_dir, paste0(spec$id, "_all_genes.tsv")),
    deg = file.path(contrast_dir, paste0(spec$id, "_DEG.tsv")),
    top = file.path(contrast_dir, paste0(spec$id, "_top20_by_padj.tsv")),
    ma_plot = file.path(ma_dir, "MA_plot.png"),
    volcano_plot = file.path(volcano_dir, "Volcano_plot.png"),
    heatmap_plot = file.path(heatmap_dir, "Heatmap.png")
  )
}

save_de_heatmap <- function(gene_ids, vsd_mat, metadata_tbl, annot_tbl, title_text, out_path) {
  gene_ids <- unique(gene_ids)
  gene_ids <- gene_ids[gene_ids %in% rownames(vsd_mat)]
  if (length(gene_ids) <= 1) {
    return(FALSE)
  }
  heat_mat <- vsd_mat[gene_ids, , drop = FALSE]
  heat_mat_scaled <- t(scale(t(heat_mat)))
  heat_mat_scaled[is.na(heat_mat_scaled)] <- 0
  row_annot <- annot_tbl[match(rownames(heat_mat_scaled), annot_tbl$Geneid), c("gene_name", "gene_biotype"), drop = FALSE]
  rownames(row_annot) <- rownames(heat_mat_scaled)
  if (ncol(row_annot) > 0) {
    keep_cols <- vapply(row_annot, function(x) any(!is.na(x) & x != ""), logical(1))
    row_annot <- row_annot[, keep_cols, drop = FALSE]
  }
  if (ncol(row_annot) == 0) {
    row_annot <- NULL
  }
  col_annot <- metadata_tbl[, c("genotype", "condition"), drop = FALSE]
  colnames(col_annot) <- c("Genotype", "Condition")
  pheatmap(
    heat_mat_scaled,
    annotation_col = col_annot,
    annotation_row = row_annot,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    show_rownames = FALSE,
    fontsize = 9,
    main = title_text,
    filename = out_path,
    width = 10,
    height = 11
  )
  TRUE
}

summary_rows <- list()
deg_master <- list()
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)
for (spec in contrast_specs) {
  paths <- make_contrast_paths(output_dir, spec)
  res <- get_res(dds, spec)
  res_df <- as.data.frame(res)
  res_df$Geneid <- rownames(res_df)
  res_df <- res_df[, c("Geneid", setdiff(colnames(res_df), "Geneid"))]
  res_df$DEG_status <- "ns"
  res_df$DEG_status[is.na(res_df$padj)] <- "not_tested"
  res_df$DEG_status[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >= 1] <- "up"
  res_df$DEG_status[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "down"
  res_df <- res_df[order(res_df$padj), ]
  res_df <- annotate_res_table(res_df, gene_annotation)

  write.table(res_df, paths$all_genes, sep = "\t", row.names = FALSE, quote = FALSE)

  deg_df <- res_df[res_df$DEG_status %in% c("up", "down"), ]
  write.table(deg_df, paths$deg, sep = "\t", row.names = FALSE, quote = FALSE)
  if (nrow(deg_df) > 0) {
    deg_master[[spec$id]] <- data.frame(
      contrast_id = spec$id,
      contrast_title = spec$title,
      Geneid = deg_df$Geneid,
      gene_id = deg_df$gene_id,
      gene_name = deg_df$gene_name,
      gene_biotype = deg_df$gene_biotype,
      log2FoldChange = deg_df$log2FoldChange,
      lfcSE = deg_df$lfcSE,
      pvalue = deg_df$pvalue,
      padj = deg_df$padj,
      DEG_status = deg_df$DEG_status,
      stringsAsFactors = FALSE
    )
  }

  top_n <- min(20, nrow(deg_df))
  top_df <- if (top_n > 0) deg_df[order(deg_df$padj, -abs(deg_df$log2FoldChange)), , drop = FALSE][seq_len(top_n), , drop = FALSE] else deg_df
  write.table(top_df, paths$top, sep = "\t", row.names = FALSE, quote = FALSE)

  png(filename = paths$ma_plot, width = 1200, height = 900, res = 150)
  plotMA(res, alpha = 0.05, main = spec$plot_title, ylim = c(-8, 8))
  dev.off()

  make_volcano(res_df, spec$plot_title, paths$volcano_plot)

  save_de_heatmap(
    gene_ids = top_df$Geneid,
    vsd_mat = vsd_mat,
    metadata_tbl = metadata,
    annot_tbl = gene_annotation,
    title_text = spec$plot_title,
    out_path = paths$heatmap_plot
  )

  summary_rows[[spec$id]] <- data.frame(
    contrast_id = spec$id,
    title = spec$title,
    tested_genes = sum(!is.na(res_df$pvalue)),
    padj_significant = sum(!is.na(res_df$padj) & res_df$padj < 0.05),
    deg_up = sum(res_df$DEG_status == "up"),
    deg_down = sum(res_df$DEG_status == "down"),
    stringsAsFactors = FALSE
  )
}

de_summary <- do.call(rbind, summary_rows)
write.table(de_summary, file.path(output_dir, "DEG_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
deg_master_table <- data.frame()
if (length(deg_master) > 0) {
  deg_master_table <- do.call(rbind, deg_master)
  write.table(deg_master_table, file.path(output_dir, "DEG_significant_all_contrasts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

top_gene_ids <- character(0)
if (nrow(deg_master_table) > 0) {
  ranked_deg <- deg_master_table[order(deg_master_table$padj, -abs(deg_master_table$log2FoldChange)), , drop = FALSE]
  ranked_deg <- ranked_deg[!duplicated(ranked_deg$Geneid), , drop = FALSE]
  top_gene_ids <- head(ranked_deg$Geneid, n = min(20, nrow(ranked_deg)))
}
top_gene_ids <- top_gene_ids[top_gene_ids %in% rownames(vsd_mat)]
if (save_de_heatmap(
  gene_ids = top_gene_ids,
  vsd_mat = vsd_mat,
  metadata_tbl = metadata,
  annot_tbl = gene_annotation,
  title_text = "Top-20 DEG across C1-C5",
  out_path = file.path(output_dir, "DEG_top_genes_heatmap.png")
)) {
  top_genes_table <- data.frame(
    Geneid = top_gene_ids,
    stringsAsFactors = FALSE
  )
  top_genes_table <- merge(top_genes_table, gene_annotation, by = "Geneid", all.x = TRUE, sort = FALSE)
  write.table(top_genes_table, file.path(output_dir, "DEG_top_genes_for_heatmap.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

writeLines("DE-анализ завершен. Результаты сохранены в папке DE_results.")
