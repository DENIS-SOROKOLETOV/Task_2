required_packages <- c("KEGGREST", "ggplot2")
missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  stop(sprintf("Не установлены пакеты: %s", paste(missing_required, collapse = ", ")))
}

suppressPackageStartupMessages({
  library(KEGGREST)
  library(ggplot2)
})
options(timeout = 300)

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
  script_dir <- dirname(script_path)
} else {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
de_dir <- file.path(project_root, "DE_results")
venn_dir <- file.path(project_root, "Venn_results")
output_dir <- file.path(project_root, "Enrichment_results", "KEGG")
venn_output_dir <- file.path(output_dir, "CommonUnique_KEGG")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(venn_output_dir, recursive = TRUE, showWarnings = FALSE)

deg_files <- list.files(de_dir, pattern = "_DEG\\.tsv$", recursive = TRUE, full.names = TRUE)
if (length(deg_files) == 0) {
  stop("В папке DE_results не найдены файлы *_DEG.tsv")
}
all_gene_files <- list.files(de_dir, pattern = "_all_genes\\.tsv$", recursive = TRUE, full.names = TRUE)
all_gene_lookup <- setNames(all_gene_files, sub("_all_genes\\.tsv$", "", basename(all_gene_files)))

format_contrast_label <- function(x) {
  out <- sub("^C[0-9]+_", "", x)
  out <- gsub("DoubleMutant", "ΔΔ", out, fixed = TRUE)
  out <- gsub("_vs_", " vs ", out, fixed = TRUE)
  out <- gsub("_", " ", out, fixed = TRUE)
  trimws(out)
}

format_set_label <- function(x) {
  if (identical(x, "venn_C1_WT_NPN_vs_WT_PN_C2_DoubleMutant_NPN_vs_DoubleMutant_PN_intersection")) {
    return("WT NPN vs WT PN ∩ ΔΔ NPN vs ΔΔ PN")
  }
  if (identical(x, "venn_C1_WT_NPN_vs_WT_PN_C2_DoubleMutant_NPN_vs_DoubleMutant_PN_only_WT_NPN_vs_PN")) {
    return("WT NPN vs WT PN ∖ ΔΔ NPN vs ΔΔ PN")
  }
  if (identical(x, "venn_C1_WT_NPN_vs_WT_PN_C2_DoubleMutant_NPN_vs_DoubleMutant_PN_only_DoubleMutant_NPN_vs_PN")) {
    return("ΔΔ NPN vs ΔΔ PN ∖ WT NPN vs WT PN")
  }
  if (identical(x, "unique_core_C1_WT_NPN_vs_WT_PN")) {
    return("C1 ∖ (C2 ∪ C3 ∪ C4)")
  }
  if (identical(x, "unique_core_C2_DoubleMutant_NPN_vs_DoubleMutant_PN")) {
    return("C2 ∖ (C1 ∪ C3 ∪ C4)")
  }
  out <- gsub("DoubleMutant", "ΔΔ", x, fixed = TRUE)
  out <- gsub("_", " ", out, fixed = TRUE)
  trimws(out)
}

normalize_gene_ids <- function(x) {
  x <- gsub("^gene:", "", x)
  x <- gsub("^fgr:", "", x)
  x <- trimws(x)
  x[!is.na(x) & nzchar(x)]
}

build_kegg_mapping <- function() {
  links <- keggLink("pathway", "fgr")
  term2gene <- unique(data.frame(
    term_id = gsub("^path:", "", unname(links)),
    gene_id = gsub("^fgr:", "", names(links)),
    stringsAsFactors = FALSE
  ))
  path_names <- keggList("pathway", "fgr")
  term2name <- unique(data.frame(
    term_id = gsub("^path:", "", names(path_names)),
    term_name = unname(path_names),
    stringsAsFactors = FALSE
  ))
  list(term2gene = term2gene, term2name = term2name)
}

run_enrichment <- function(gene_set, universe_set, term2gene, term2name, min_term_size = 5) {
  gene_set <- unique(gene_set)
  universe_set <- unique(universe_set)
  term2gene <- unique(term2gene)
  term2gene <- term2gene[term2gene$gene_id %in% universe_set, , drop = FALSE]
  gene_set <- gene_set[gene_set %in% universe_set]
  if (length(gene_set) == 0 || nrow(term2gene) == 0) {
    return(data.frame())
  }

  term_genes <- split(term2gene$gene_id, term2gene$term_id)
  term_genes <- lapply(term_genes, unique)
  term_sizes <- vapply(term_genes, length, integer(1))
  valid_terms <- names(term_sizes)[term_sizes >= min_term_size]
  if (length(valid_terms) == 0) {
    return(data.frame())
  }

  bg_size <- length(universe_set)
  set_size <- length(gene_set)
  out_list <- vector("list", length(valid_terms))
  for (i in seq_along(valid_terms)) {
    term_id <- valid_terms[i]
    tg <- term_genes[[term_id]]
    a <- length(intersect(gene_set, tg))
    b <- length(setdiff(tg, gene_set))
    c <- length(setdiff(gene_set, tg))
    d <- bg_size - a - b - c
    if (a == 0 || d < 0) {
      out_list[[i]] <- NULL
      next
    }
    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    p <- fisher.test(mat, alternative = "greater")$p.value
    out_list[[i]] <- data.frame(
      term_id = term_id,
      overlap = a,
      term_size = length(tg),
      set_size = set_size,
      bg_size = bg_size,
      pvalue = p,
      gene_ratio = a / set_size,
      bg_ratio = length(tg) / bg_size,
      overlap_genes = paste(sort(intersect(gene_set, tg)), collapse = ","),
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, out_list)
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame())
  }
  out$padj <- p.adjust(out$pvalue, method = "BH")
  out <- merge(out, term2name, by = "term_id", all.x = TRUE, sort = FALSE)
  out <- out[order(out$padj, out$pvalue, -out$overlap), ]
  out
}

save_top_plot <- function(tbl, output_path, title_text, top_n = 15) {
  if (nrow(tbl) == 0) {
    if (file.exists(output_path)) {
      unlink(output_path, force = TRUE)
    }
    return(invisible(NULL))
  }
  tbl <- tbl[!is.na(tbl$padj) & tbl$padj < 0.05, , drop = FALSE]
  if (nrow(tbl) == 0) {
    if (file.exists(output_path)) {
      unlink(output_path, force = TRUE)
    }
    return(invisible(NULL))
  }
  top_tbl <- head(tbl[order(tbl$padj, tbl$pvalue), ], n = min(top_n, nrow(tbl)))
  top_tbl$label <- ifelse(is.na(top_tbl$term_name) | top_tbl$term_name == "", top_tbl$term_id, top_tbl$term_name)
  top_tbl$score <- -log10(pmax(top_tbl$padj, 1e-300))
  p <- ggplot(top_tbl, aes(x = .data$score, y = reorder(.data$label, .data$score))) +
    geom_point(aes(size = .data$overlap, color = .data$gene_ratio), alpha = 0.9) +
    scale_color_gradient(low = "#377EB8", high = "#E41A1C") +
    labs(
      title = title_text,
      x = "-log10(FDR)",
      y = "",
      color = "GeneRatio",
      size = "Overlap"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  ggsave(filename = output_path, plot = p, width = 11, height = 7.5, dpi = 300)
}

get_de_universe <- function(de_root) {
  all_files <- list.files(de_root, pattern = "_all_genes\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(all_files) == 0) {
    return(character(0))
  }
  all_genes <- unique(unlist(lapply(all_files, function(f) {
    tbl <- read.delim(f, sep = "\t", header = TRUE, check.names = FALSE)
    if (!("Geneid" %in% colnames(tbl))) {
      return(character(0))
    }
    normalize_gene_ids(tbl$Geneid)
  })))
  unique(all_genes)
}

read_gene_set_file <- function(path) {
  tbl <- read.delim(path, sep = "\t", header = TRUE, check.names = FALSE)
  gene_col <- intersect(colnames(tbl), c("Geneid", "gene_id", "gene"))
  if (length(gene_col) == 0) {
    return(character(0))
  }
  normalize_gene_ids(tbl[[gene_col[1]]])
}

kegg_map <- build_kegg_mapping()

summary_rows <- list()
for (deg_file in deg_files) {
  contrast_id <- sub("_DEG\\.tsv$", "", basename(deg_file))
  all_file <- all_gene_lookup[[contrast_id]]
  if (is.null(all_file) || !file.exists(all_file)) {
    next
  }

  deg_tbl <- read.delim(deg_file, sep = "\t", header = TRUE, check.names = FALSE)
  all_tbl <- read.delim(all_file, sep = "\t", header = TRUE, check.names = FALSE)
  if (!("Geneid" %in% colnames(deg_tbl)) || !("Geneid" %in% colnames(all_tbl))) {
    next
  }

  deg_genes <- normalize_gene_ids(deg_tbl$Geneid)
  universe_genes <- normalize_gene_ids(all_tbl$Geneid)
  universe_genes <- unique(universe_genes)
  deg_genes <- unique(deg_genes[deg_genes %in% universe_genes])

  kegg_res <- run_enrichment(
    gene_set = deg_genes,
    universe_set = universe_genes,
    term2gene = kegg_map$term2gene,
    term2name = kegg_map$term2name,
    min_term_size = 5
  )

  kegg_out <- file.path(output_dir, paste0(contrast_id, "_KEGG_enrichment.tsv"))
  write.table(kegg_res, kegg_out, sep = "\t", row.names = FALSE, quote = FALSE)

  save_top_plot(
    tbl = kegg_res,
    output_path = file.path(output_dir, paste0(contrast_id, "_KEGG_top_terms.png")),
    title_text = paste0("KEGG-обогащение DEG: ", format_contrast_label(contrast_id))
  )

  summary_rows[[contrast_id]] <- data.frame(
    contrast_id = contrast_id,
    n_deg = length(deg_genes),
    n_background = length(universe_genes),
    n_kegg_enriched_fdr_0_05 = sum(kegg_res$padj < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

summary_tbl <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame()
write.table(summary_tbl, file.path(output_dir, "KEGG_enrichment.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

venn_summary_rows <- list()
de_universe <- get_de_universe(de_dir)
if (length(de_universe) > 0 && dir.exists(venn_dir)) {
  venn_set_files <- c(
    list.files(venn_dir, pattern = "_intersection\\.tsv$", full.names = TRUE),
    list.files(venn_dir, pattern = "_only_.*\\.tsv$", full.names = TRUE),
    file.path(venn_dir, "common_genes_all_core_contrasts.tsv")
  )
  venn_set_files <- unique(venn_set_files[file.exists(venn_set_files)])
  venn_set_files <- venn_set_files[grepl("C1_|C2_", basename(venn_set_files))]

  for (set_file in venn_set_files) {
    set_id <- sub("\\.tsv$", "", basename(set_file))
    set_genes <- unique(read_gene_set_file(set_file))
    set_genes <- set_genes[set_genes %in% de_universe]
    if (length(set_genes) == 0) {
      next
    }

    kegg_res <- run_enrichment(
      gene_set = set_genes,
      universe_set = de_universe,
      term2gene = kegg_map$term2gene,
      term2name = kegg_map$term2name,
      min_term_size = 5
    )

    write.table(
      kegg_res,
      file.path(venn_output_dir, paste0(set_id, "_KEGG_enrichment.tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    save_top_plot(
      tbl = kegg_res,
      output_path = file.path(venn_output_dir, paste0(set_id, "_KEGG_top_terms.png")),
      title_text = paste0("KEGG-обогащение: ", format_set_label(set_id))
    )

    venn_summary_rows[[set_id]] <- data.frame(
      set_id = set_id,
      n_genes_in_set = length(set_genes),
      n_background = length(de_universe),
      n_kegg_enriched_fdr_0_05 = sum(kegg_res$padj < 0.05, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  core_unique_path <- file.path(venn_dir, "unique_genes_by_core_contrast.tsv")
  if (file.exists(core_unique_path)) {
    core_unique_tbl <- read.delim(core_unique_path, sep = "\t", header = TRUE, check.names = FALSE)
    if (all(c("contrast_id", "Geneid") %in% colnames(core_unique_tbl))) {
      contrast_ids <- unique(core_unique_tbl$contrast_id)
      contrast_ids <- contrast_ids[grepl("^C[12]_", contrast_ids)]
      for (cid in contrast_ids) {
        cid_genes <- normalize_gene_ids(core_unique_tbl$Geneid[core_unique_tbl$contrast_id == cid])
        cid_genes <- unique(cid_genes[cid_genes %in% de_universe])
        if (length(cid_genes) == 0) {
          next
        }
        set_id <- paste0("unique_core_", cid)
        kegg_res <- run_enrichment(
          gene_set = cid_genes,
          universe_set = de_universe,
          term2gene = kegg_map$term2gene,
          term2name = kegg_map$term2name,
          min_term_size = 5
        )
        write.table(
          kegg_res,
          file.path(venn_output_dir, paste0(set_id, "_KEGG_enrichment.tsv")),
          sep = "\t",
          row.names = FALSE,
          quote = FALSE
        )
        save_top_plot(
          tbl = kegg_res,
          output_path = file.path(venn_output_dir, paste0(set_id, "_KEGG_top_terms.png")),
          title_text = paste0("KEGG-обогащение: ", format_set_label(set_id))
        )
        venn_summary_rows[[set_id]] <- data.frame(
          set_id = set_id,
          n_genes_in_set = length(cid_genes),
          n_background = length(de_universe),
          n_kegg_enriched_fdr_0_05 = sum(kegg_res$padj < 0.05, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

venn_summary_tbl <- if (length(venn_summary_rows) > 0) do.call(rbind, venn_summary_rows) else data.frame()
write.table(
  venn_summary_tbl,
  file.path(venn_output_dir, "venn_sets_enrichment_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

writeLines("Enrichment-анализ завершен. Результаты сохранены в папке Enrichment_results/KEGG.")
