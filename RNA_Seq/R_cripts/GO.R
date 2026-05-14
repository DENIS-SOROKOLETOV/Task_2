required_packages <- c("ggplot2", "jsonlite")
missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  stop(sprintf("Не установлены пакеты: %s", paste(missing_required, collapse = ", ")))
}

suppressPackageStartupMessages({
  library(ggplot2)
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
de_dir <- file.path(project_root, "DE_results")
reference_dir <- file.path(project_root, "Reference_genome")
output_dir <- file.path(project_root, "Enrichment_results", "GO")
venn_dir <- file.path(project_root, "Venn_results")
common_unique_output_dir <- file.path(output_dir, "CommonUnique_GO")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(common_unique_output_dir, recursive = TRUE, showWarnings = FALSE)
old_common_unique <- list.files(common_unique_output_dir, full.names = TRUE)
if (length(old_common_unique) > 0) {
  unlink(old_common_unique, recursive = TRUE, force = TRUE)
}

gff_path <- file.path(reference_dir, "GCF_000240135.3_ASM24013v3_genomic.gff")
gaf_path <- file.path(reference_dir, "GCF_000240135.3_ASM24013v3_gene_ontology.gaf")
if (!file.exists(gff_path)) {
  stop("Не найден файл аннотации генов GFF.")
}
if (!file.exists(gaf_path)) {
  stop("Не найден файл GO-аннотаций GAF.")
}

parse_attr_value <- function(attr, key) {
  pattern <- paste0(".*", key, "=([^;]+).*")
  hit <- grepl(paste0(key, "="), attr)
  out <- rep(NA_character_, length(attr))
  out[hit] <- sub(pattern, "\\1", attr[hit])
  out
}

extract_geneid_num <- function(attr) {
  hit <- grepl("GeneID:[0-9]+", attr)
  out <- rep(NA_character_, length(attr))
  out[hit] <- sub(".*GeneID:([0-9]+).*", "\\1", attr[hit])
  out
}

gff <- read.delim(
  gff_path,
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "#",
  stringsAsFactors = FALSE
)
colnames(gff) <- c("seqid", "source", "feature", "start", "end", "score", "strand", "phase", "attribute")
gff_gene <- gff[gff$feature == "gene", c("feature", "attribute"), drop = FALSE]
gene_symbol <- parse_attr_value(gff_gene$attribute, "Name")
geneid_num <- extract_geneid_num(gff_gene$attribute)
gene_map <- data.frame(
  gene_symbol = gene_symbol,
  geneid_num = geneid_num,
  stringsAsFactors = FALSE
)
gene_map <- gene_map[!is.na(gene_map$gene_symbol) & !is.na(gene_map$geneid_num), , drop = FALSE]
gene_map <- gene_map[!duplicated(gene_map$gene_symbol), , drop = FALSE]
gene_map$de_geneid <- paste0("gene:", gene_map$gene_symbol)

gaf <- read.delim(
  gaf_path,
  sep = "\t",
  header = FALSE,
  comment.char = "!",
  quote = "",
  stringsAsFactors = FALSE
)
if (ncol(gaf) < 9) {
  stop("Формат GAF не распознан: слишком мало колонок.")
}
colnames(gaf)[1:9] <- c("DB", "GeneID", "Symbol", "Qualifier", "GO_ID", "Reference", "Evidence_Code", "With_From", "Aspect")
gaf <- gaf[, c("GeneID", "GO_ID", "Aspect"), drop = FALSE]
gaf$GeneID <- as.character(gaf$GeneID)
gaf$GO_ID <- as.character(gaf$GO_ID)
gaf$Aspect <- as.character(gaf$Aspect)
gaf <- gaf[!is.na(gaf$GeneID) & !is.na(gaf$GO_ID) & !is.na(gaf$Aspect), , drop = FALSE]
gaf <- gaf[gaf$GO_ID != "" & gaf$GeneID != "" & gaf$Aspect != "", , drop = FALSE]
gaf <- gaf[!duplicated(gaf), , drop = FALSE]

aspect_to_ontology <- function(x) {
  out <- x
  out[x == "P"] <- "biological_process"
  out[x == "C"] <- "cellular_component"
  out[x == "F"] <- "molecular_function"
  out[!(x %in% c("P", "C", "F"))] <- "unknown"
  out
}
gaf$Ontology <- aspect_to_ontology(gaf$Aspect)

go2genes <- split(gaf$GeneID, gaf$GO_ID)
go2genes <- lapply(go2genes, unique)
go2ontology <- tapply(gaf$Ontology, gaf$GO_ID, function(x) unique(x)[1])

background_gene_ids <- sort(unique(gaf$GeneID))
if (length(background_gene_ids) == 0) {
  stop("Не удалось сформировать фоновый набор генов из GAF.")
}

deg_files <- list.files(de_dir, pattern = "_DEG\\.tsv$", recursive = TRUE, full.names = TRUE)
if (length(deg_files) == 0) {
  stop("Не найдены DEG-файлы в DE_results.")
}

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
  out <- gsub("DoubleMutant", "ΔΔ", x, fixed = TRUE)
  out <- gsub("_", " ", out, fixed = TRUE)
  trimws(out)
}

format_de_gene <- function(x) {
  y <- as.character(x)
  y <- sub("^gene:", "", y)
  y <- sub("^gene-", "", y)
  y
}

fetch_go_names_quickgo <- function(go_ids, batch_size = 200L) {
  go_ids <- unique(go_ids)
  out <- data.frame(GO_ID = character(0), term_name = character(0), stringsAsFactors = FALSE)
  if (length(go_ids) == 0) {
    return(out)
  }
  idx <- seq(1, length(go_ids), by = batch_size)
  for (i in idx) {
    chunk <- go_ids[i:min(i + batch_size - 1, length(go_ids))]
    url <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", paste(chunk, collapse = ","))
    batch_ok <- FALSE
    try({
      raw_txt <- readLines(url, warn = FALSE, encoding = "UTF-8")
      payload <- jsonlite::fromJSON(paste(raw_txt, collapse = ""))
      if (!is.null(payload$results) && nrow(payload$results) > 0) {
        tmp <- data.frame(
          GO_ID = as.character(payload$results$id),
          term_name = as.character(payload$results$name),
          stringsAsFactors = FALSE
        )
        out <- rbind(out, tmp)
        batch_ok <- TRUE
      }
    }, silent = TRUE)
    if (!batch_ok) {
      next
    }
    Sys.sleep(0.1)
  }
  out <- out[!duplicated(out$GO_ID), , drop = FALSE]
  out
}

go_name_cache_path <- file.path(reference_dir, "go_term_names_cache.tsv")
go_ids_all <- sort(unique(gaf$GO_ID))
cached_names <- data.frame(GO_ID = character(0), term_name = character(0), stringsAsFactors = FALSE)
if (file.exists(go_name_cache_path)) {
  cached_names <- read.delim(go_name_cache_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("GO_ID", "term_name") %in% colnames(cached_names))) {
    cached_names <- data.frame(GO_ID = character(0), term_name = character(0), stringsAsFactors = FALSE)
  } else {
    cached_names <- cached_names[, c("GO_ID", "term_name"), drop = FALSE]
  }
}
missing_go_ids <- setdiff(go_ids_all, cached_names$GO_ID)
if (length(missing_go_ids) > 0) {
  fetched_names <- fetch_go_names_quickgo(missing_go_ids, batch_size = 200L)
  if (nrow(fetched_names) > 0) {
    cached_names <- rbind(cached_names, fetched_names)
    cached_names <- cached_names[!duplicated(cached_names$GO_ID), , drop = FALSE]
    write.table(cached_names, go_name_cache_path, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}
go_name_table <- cached_names
go_name_map <- setNames(go_name_table$term_name, go_name_table$GO_ID)

run_go_ora <- function(foreground_gene_ids, background_gene_ids, go2genes, go2ontology, go_name_map, min_term_size = 5) {
  N <- length(background_gene_ids)
  n <- length(foreground_gene_ids)
  if (n == 0 || N == 0) {
    return(data.frame())
  }

  term_stats <- lapply(names(go2genes), function(go_id) {
    term_genes <- unique(go2genes[[go_id]])
    term_genes <- intersect(term_genes, background_gene_ids)
    K <- length(term_genes)
    if (K < min_term_size) {
      return(NULL)
    }
    overlap <- intersect(foreground_gene_ids, term_genes)
    k <- length(overlap)
    if (k == 0) {
      return(NULL)
    }
    p_value <- phyper(q = k - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
    term_name <- unname(go_name_map[go_id])
    if (length(term_name) == 0 || is.na(term_name) || term_name == "") {
      term_name <- go_id
    }
    data.frame(
      GO_ID = go_id,
      term_name = term_name,
      Ontology = as.character(go2ontology[[go_id]]),
      term_size = K,
      foreground_size = n,
      overlap_size = k,
      gene_ratio = k / n,
      background_ratio = K / N,
      pvalue = p_value,
      overlap_genes = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  term_stats <- term_stats[!vapply(term_stats, is.null, logical(1))]
  if (length(term_stats) == 0) {
    return(data.frame())
  }
  res <- do.call(rbind, term_stats)
  res$qvalue <- p.adjust(res$pvalue, method = "BH")
  res$log10_qvalue <- -log10(pmax(res$qvalue, .Machine$double.xmin))
  res <- res[order(res$qvalue, res$pvalue), , drop = FALSE]
  rownames(res) <- NULL
  res
}

make_go_plot <- function(res_tbl, title_text, out_path, top_n_per_ontology = 10) {
  if (nrow(res_tbl) == 0) {
    return(FALSE)
  }
  res_sig <- res_tbl[res_tbl$qvalue < 0.05, , drop = FALSE]
  if (nrow(res_sig) == 0) {
    return(FALSE)
  }

  split_by_ont <- split(res_sig, res_sig$Ontology)
  top_tbl <- do.call(rbind, lapply(split_by_ont, function(df) {
    df <- df[order(df$qvalue, -df$overlap_size), , drop = FALSE]
    head(df, n = min(top_n_per_ontology, nrow(df)))
  }))
  if (nrow(top_tbl) == 0) {
    return(FALSE)
  }

  top_tbl$term_label <- top_tbl$term_name
  top_tbl <- top_tbl[order(top_tbl$Ontology, top_tbl$log10_qvalue), , drop = FALSE]
  top_tbl$term_label <- factor(top_tbl$term_label, levels = top_tbl$term_label)

  p <- ggplot(top_tbl, aes(x = log10_qvalue, y = term_label, fill = Ontology)) +
    geom_col(width = 0.85) +
    scale_fill_manual(values = c(
      biological_process = "#F8766D",
      cellular_component = "#00BA38",
      molecular_function = "#619CFF",
      unknown = "#B3B3B3"
    )) +
    labs(
      title = title_text,
      x = expression(-log[10](qvalue)),
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
    )
  ggsave(out_path, p, width = 13, height = 9, dpi = 300)
  TRUE
}

safe_id <- function(x) {
  out <- gsub("[^A-Za-z0-9_]+", "_", x)
  out <- gsub("_+", "_", out)
  out <- gsub("^_|_$", "", out)
  out
}

analyze_go_for_gene_set <- function(gene_ids, set_id, title_text, out_dir) {
  input_ids <- unique(gene_ids[!is.na(gene_ids) & gene_ids != ""])
  input_symbols <- format_de_gene(input_ids)
  mapped_geneid_num <- gene_map$geneid_num[match(input_symbols, gene_map$gene_symbol)]
  foreground_gene_ids <- sort(unique(mapped_geneid_num[!is.na(mapped_geneid_num)]))
  foreground_gene_ids <- intersect(foreground_gene_ids, background_gene_ids)

  go_res <- run_go_ora(
    foreground_gene_ids = foreground_gene_ids,
    background_gene_ids = background_gene_ids,
    go2genes = go2genes,
    go2ontology = go2ontology,
    go_name_map = go_name_map,
    min_term_size = 5
  )

  file_prefix <- safe_id(set_id)
  out_table <- file.path(out_dir, paste0(file_prefix, "_GO_enrichment.tsv"))
  if (nrow(go_res) == 0) {
    write.table(data.frame(), out_table, sep = "\t", row.names = FALSE, quote = FALSE)
    return(data.frame(
      set_id = set_id,
      n_input_genes = length(input_ids),
      n_genes_with_go_mapping = length(foreground_gene_ids),
      n_significant_terms = 0,
      plot_created = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  write.table(go_res, out_table, sep = "\t", row.names = FALSE, quote = FALSE)
  plot_path <- file.path(out_dir, paste0(file_prefix, "_GO_barplot.png"))
  has_plot <- make_go_plot(
    res_tbl = go_res,
    title_text = title_text,
    out_path = plot_path,
    top_n_per_ontology = 10
  )

  data.frame(
    set_id = set_id,
    n_input_genes = length(input_ids),
    n_genes_with_go_mapping = length(foreground_gene_ids),
    n_significant_terms = sum(go_res$qvalue < 0.05, na.rm = TRUE),
    plot_created = has_plot,
    stringsAsFactors = FALSE
  )
}

summary_rows <- list()
for (deg_file in deg_files) {
  contrast_id <- sub("_DEG\\.tsv$", "", basename(deg_file))
  deg_tbl <- read.delim(deg_file, sep = "\t", header = TRUE, check.names = FALSE)
  if (!("Geneid" %in% colnames(deg_tbl))) {
    next
  }

  input_ids <- unique(deg_tbl$Geneid[!is.na(deg_tbl$Geneid) & deg_tbl$Geneid != ""])
  input_symbols <- format_de_gene(input_ids)
  mapped_geneid_num <- gene_map$geneid_num[match(input_symbols, gene_map$gene_symbol)]
  foreground_gene_ids <- sort(unique(mapped_geneid_num[!is.na(mapped_geneid_num)]))
  foreground_gene_ids <- intersect(foreground_gene_ids, background_gene_ids)

  go_res <- run_go_ora(
    foreground_gene_ids = foreground_gene_ids,
    background_gene_ids = background_gene_ids,
    go2genes = go2genes,
    go2ontology = go2ontology,
    go_name_map = go_name_map,
    min_term_size = 5
  )

  out_table <- file.path(output_dir, paste0(contrast_id, "_GO_enrichment.tsv"))
  if (nrow(go_res) == 0) {
    write.table(data.frame(), out_table, sep = "\t", row.names = FALSE, quote = FALSE)
    summary_rows[[contrast_id]] <- data.frame(
      contrast_id = contrast_id,
      n_input_deg = length(input_ids),
      n_deg_with_go_mapping = length(foreground_gene_ids),
      n_significant_terms = 0,
      stringsAsFactors = FALSE
    )
    next
  }

  write.table(go_res, out_table, sep = "\t", row.names = FALSE, quote = FALSE)

  plot_path <- file.path(output_dir, paste0(contrast_id, "_GO_barplot.png"))
  has_plot <- make_go_plot(
    res_tbl = go_res,
    title_text = paste0("GO-обогащение DEG: ", format_contrast_label(contrast_id)),
    out_path = plot_path,
    top_n_per_ontology = 10
  )

  summary_rows[[contrast_id]] <- data.frame(
    contrast_id = contrast_id,
    n_input_deg = length(input_ids),
    n_deg_with_go_mapping = length(foreground_gene_ids),
    n_significant_terms = sum(go_res$qvalue < 0.05, na.rm = TRUE),
    plot_created = has_plot,
    stringsAsFactors = FALSE
  )
}

summary_tbl <- do.call(rbind, summary_rows)
write.table(summary_tbl, file.path(output_dir, "GO_enrichment_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

common_unique_summary_rows <- list()
if (dir.exists(venn_dir)) {
  pair_set_files <- c(
    file.path(venn_dir, "venn_C1_WT_NPN_vs_WT_PN_C2_DoubleMutant_NPN_vs_DoubleMutant_PN_intersection.tsv"),
    file.path(venn_dir, "venn_C1_WT_NPN_vs_WT_PN_C2_DoubleMutant_NPN_vs_DoubleMutant_PN_only_DoubleMutant_NPN_vs_PN.tsv"),
    file.path(venn_dir, "venn_C1_WT_NPN_vs_WT_PN_C2_DoubleMutant_NPN_vs_DoubleMutant_PN_only_WT_NPN_vs_PN.tsv")
  )
  pair_set_files <- pair_set_files[file.exists(pair_set_files)]
  if (length(pair_set_files) > 0) {
    for (set_file in pair_set_files) {
      set_tbl <- read.delim(set_file, sep = "\t", header = TRUE, check.names = FALSE)
      if (!("Geneid" %in% colnames(set_tbl))) {
        next
      }
      set_name <- sub("\\.tsv$", "", basename(set_file))
      common_unique_summary_rows[[set_name]] <- analyze_go_for_gene_set(
        gene_ids = set_tbl$Geneid,
        set_id = set_name,
        title_text = paste0("GO-обогащение: ", format_set_label(set_name)),
        out_dir = common_unique_output_dir
      )
    }
  }

  common_files <- c()
  for (cf in common_files) {
    if (!file.exists(cf)) {
      next
    }
    common_tbl <- read.delim(cf, sep = "\t", header = TRUE, check.names = FALSE)
    if (!("Geneid" %in% colnames(common_tbl))) {
      next
    }
    set_name <- sub("\\.tsv$", "", basename(cf))
    common_unique_summary_rows[[set_name]] <- analyze_go_for_gene_set(
      gene_ids = common_tbl$Geneid,
      set_id = set_name,
      title_text = paste0("GO-обогащение: ", format_set_label(set_name)),
      out_dir = common_unique_output_dir
    )
  }

  unique_files <- c(
    file.path(venn_dir, "unique_genes_by_core_contrast.tsv"),
    file.path(venn_dir, "unique_genes_by_contrast.tsv")
  )
  for (uf in unique_files) {
    if (!file.exists(uf)) {
      next
    }
    unique_tbl <- read.delim(uf, sep = "\t", header = TRUE, check.names = FALSE)
    if (!all(c("contrast_id", "Geneid") %in% colnames(unique_tbl))) {
      next
    }
    unique_tbl <- unique_tbl[!is.na(unique_tbl$contrast_id) & !is.na(unique_tbl$Geneid), , drop = FALSE]
    split_sets <- split(unique_tbl$Geneid, unique_tbl$contrast_id)
    for (cid in names(split_sets)) {
      if (!(cid %in% c("C1_WT_NPN_vs_WT_PN", "C2_DoubleMutant_NPN_vs_DoubleMutant_PN"))) {
        next
      }
      set_name <- paste0("unique_", cid)
      common_unique_summary_rows[[set_name]] <- analyze_go_for_gene_set(
        gene_ids = split_sets[[cid]],
        set_id = set_name,
        title_text = paste0("GO-обогащение: ", format_set_label(set_name)),
        out_dir = common_unique_output_dir
      )
    }
  }
}

if (length(common_unique_summary_rows) > 0) {
  common_unique_summary_tbl <- do.call(rbind, common_unique_summary_rows)
} else {
  common_unique_summary_tbl <- data.frame(
    set_id = character(0),
    n_input_genes = integer(0),
    n_genes_with_go_mapping = integer(0),
    n_significant_terms = integer(0),
    plot_created = logical(0),
    stringsAsFactors = FALSE
  )
}
write.table(
  common_unique_summary_tbl,
  file.path(common_unique_output_dir, "GO_common_unique_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

writeLines("GO-анализ завершен. Результаты сохранены в папке Enrichment_results/GO.")
