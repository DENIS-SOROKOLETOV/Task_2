required_packages <- c("data.table", "DSS", "GenomicRanges", "GenomicFeatures")
missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  stop(sprintf("Не установлены пакеты: %s", paste(missing_required, collapse = ", ")))
}

suppressPackageStartupMessages({
  library(data.table)
  library(DSS)
  library(GenomicRanges)
  library(GenomicFeatures)
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
output_dir <- file.path(wgbs_root, "DMR_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

annotation_gff <- normalizePath(
  file.path(wgbs_root, "..", "RNA_Seq", "Reference_genome", "Fusarium_graminearum_ph_1_gca_000240135.ASM24013v3.62.gff3"),
  winslash = "/",
  mustWork = TRUE
)

sample_annotation <- data.frame(
  sample_id = c("SRR10382875", "SRR10382876", "SRR10382877", "SRR10382878"),
  genotype = c("DoubleMutant", "DoubleMutant", "WT", "WT"),
  condition = c("NPN", "PN", "NPN", "PN"),
  stringsAsFactors = FALSE
)

write.table(
  sample_annotation,
  file = file.path(output_dir, "sample_annotation.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

context_params <- list(
  CpG = list(min_coverage = 3L, p_threshold = 0.05, delta = 0.05, minlen = 30L, min_sites = 2L, merge_dist = 200L, pct_sig = 0.30),
  CHG = list(min_coverage = 3L, p_threshold = 0.05, delta = 0.05, minlen = 30L, min_sites = 2L, merge_dist = 200L, pct_sig = 0.30)
)

contrast_specs <- list(
  list(
    id = "C1_WT_NPN_vs_WT_PN",
    title = "Эффект среды в WT: NPN vs PN",
    method = "pairwise",
    group1 = "SRR10382877",
    group2 = "SRR10382878"
  ),
  list(
    id = "C2_DoubleMutant_NPN_vs_DoubleMutant_PN",
    title = "Эффект среды в ΔΔ: NPN vs PN",
    method = "pairwise",
    group1 = "SRR10382875",
    group2 = "SRR10382876"
  ),
  list(
    id = "C3_DoubleMutant_PN_vs_WT_PN",
    title = "Эффект генотипа в PN: ΔΔ vs WT",
    method = "pairwise",
    group1 = "SRR10382876",
    group2 = "SRR10382878"
  ),
  list(
    id = "C4_DoubleMutant_NPN_vs_WT_NPN",
    title = "Эффект генотипа в NPN: ΔΔ vs WT",
    method = "pairwise",
    group1 = "SRR10382875",
    group2 = "SRR10382877"
  ),
  list(
    id = "P1_pooled_genotype_DoubleMutant_vs_WT",
    title = "Pooled эффект генотипа: ΔΔ vs WT (PN+NPN)",
    method = "pooled",
    group1 = c("SRR10382875", "SRR10382876"),
    group2 = c("SRR10382877", "SRR10382878")
  ),
  list(
    id = "P2_pooled_condition_NPN_vs_PN",
    title = "Pooled эффект среды: NPN vs PN (WT+ΔΔ)",
    method = "pooled",
    group1 = c("SRR10382875", "SRR10382877"),
    group2 = c("SRR10382876", "SRR10382878")
  )
)

find_bedgraph <- function(sample_id, context) {
  candidates <- c(
    file.path(extract_root, sample_id, paste0(sample_id, "_", context, ".bedGraph")),
    file.path(extract_root, sample_id, paste0(sample_id, "_", context, ".bedGraph.gz"))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[1]
}

read_methylation_counts <- function(path, min_coverage) {
  dt <- fread(
    file = path,
    sep = "\t",
    header = FALSE,
    skip = "track",
    select = c(1, 2, 5, 6),
    col.names = c("chr", "start", "methylated", "unmethylated"),
    showProgress = FALSE
  )
  dt[, pos := as.integer(start) + 1L]
  dt[, N := as.numeric(methylated) + as.numeric(unmethylated)]
  dt[, X := as.numeric(methylated)]
  dt <- dt[!is.na(N) & !is.na(X) & N >= min_coverage, .(chr, pos, N, X)]
  as.data.frame(dt)
}

run_dmr_pooled <- function(context, group1, group2, params) {
  load_pooled <- function(sample_ids) {
    combined <- NULL
    for (sid in sample_ids) {
      path <- find_bedgraph(sid, context)
      if (is.na(path)) {
        stop(sprintf("Не найден bedGraph %s для образца %s", context, sid))
      }
      dt <- as.data.table(read_methylation_counts(path, min_coverage = params$min_coverage))
      combined <- rbind(combined, dt)
    }
    combined[, .(N = sum(N), X = sum(X)), by = .(chr, pos)]
  }
  dt1 <- load_pooled(group1)
  dt2 <- load_pooled(group2)
  setnames(dt1, c("N", "X"), c("N1", "X1"))
  setnames(dt2, c("N", "X"), c("N2", "X2"))
  merged <- merge(dt1, dt2, by = c("chr", "pos"), all = FALSE)
  if (nrow(merged) == 0) {
    return(data.frame())
  }
  merged[, mu1 := X1 / N1]
  merged[, mu2 := X2 / N2]
  merged[, diff := mu1 - mu2]
  merged <- merged[abs(diff) >= params$delta]
  if (nrow(merged) == 0) {
    return(data.frame())
  }
  merged[, p_pool := (X1 + X2) / (N1 + N2)]
  merged[, se := sqrt(pmax(p_pool * (1 - p_pool) * (1 / N1 + 1 / N2), .Machine$double.eps))]
  merged[, stat := diff / se]
  merged[, pval := 2 * pnorm(-abs(stat))]
  merged[, diff.se := sqrt(pmax(mu1 * (1 - mu1) / N1 + mu2 * (1 - mu2) / N2, .Machine$double.eps))]
  dml <- as.data.frame(merged[, .(chr, pos, stat, pval, diff, diff.se, mu1, mu2)])
  dmr <- callDMR(
    dml,
    p.threshold = params$p_threshold,
    delta = params$delta,
    minlen = params$minlen,
    minCG = params$min_sites,
    dis.merge = params$merge_dist,
    pct.sig = params$pct_sig
  )
  if (is.null(dmr) || nrow(dmr) == 0) {
    return(data.frame())
  }
  as.data.frame(dmr)
}

run_dmr_pairwise <- function(context, sample_g1, sample_g2, params) {
  path1 <- find_bedgraph(sample_g1, context)
  path2 <- find_bedgraph(sample_g2, context)
  if (is.na(path1) || is.na(path2)) {
    stop(sprintf("Не найден bedGraph для pairwise-сравнения %s vs %s (%s)", sample_g1, sample_g2, context))
  }
  dt1 <- as.data.table(read_methylation_counts(path1, min_coverage = params$min_coverage))
  dt2 <- as.data.table(read_methylation_counts(path2, min_coverage = params$min_coverage))
  setnames(dt1, c("N", "X"), c("N1", "X1"))
  setnames(dt2, c("N", "X"), c("N2", "X2"))
  merged <- merge(dt1, dt2, by = c("chr", "pos"), all = FALSE)
  if (nrow(merged) == 0) {
    return(data.frame())
  }
  merged[, mu1 := X1 / N1]
  merged[, mu2 := X2 / N2]
  merged[, diff := mu1 - mu2]
  merged <- merged[abs(diff) >= params$delta]
  if (nrow(merged) == 0) {
    return(data.frame())
  }
  merged[, p_pool := (X1 + X2) / (N1 + N2)]
  merged[, se := sqrt(pmax(p_pool * (1 - p_pool) * (1 / N1 + 1 / N2), .Machine$double.eps))]
  merged[, stat := diff / se]
  merged[, pval := 2 * pnorm(-abs(stat))]
  merged[, diff.se := sqrt(pmax(mu1 * (1 - mu1) / N1 + mu2 * (1 - mu2) / N2, .Machine$double.eps))]
  dml <- as.data.frame(merged[, .(chr, pos, stat, pval, diff, diff.se, mu1, mu2)])
  dmr <- callDMR(
    dml,
    p.threshold = params$p_threshold,
    delta = params$delta,
    minlen = params$minlen,
    minCG = params$min_sites,
    dis.merge = params$merge_dist,
    pct.sig = params$pct_sig
  )
  if (is.null(dmr) || nrow(dmr) == 0) {
    return(data.frame())
  }
  as.data.frame(dmr)
}

run_dmr_for_contrast <- function(context, spec, params) {
  if (spec$method == "pooled") {
    dmr <- run_dmr_pooled(context, spec$group1, spec$group2, params)
  } else if (spec$method == "pairwise") {
    dmr <- run_dmr_pairwise(context, spec$group1, spec$group2, params)
  } else {
    stop(sprintf("Неизвестный method: %s", spec$method))
  }
  if (nrow(dmr) == 0) {
    return(dmr)
  }
  dmr$context <- context
  dmr$contrast_id <- spec$id
  dmr$contrast_title <- spec$title
  dmr$method <- spec$method
  dmr$group1_samples <- paste(spec$group1, collapse = ";")
  dmr$group2_samples <- paste(spec$group2, collapse = ";")
  dmr$dmr_id <- paste0(context, "_", spec$id, "_DMR_", seq_len(nrow(dmr)))
  dmr
}

txdb <- makeTxDbFromGFF(annotation_gff, format = "gff3")
genes_gr <- genes(txdb)
gene_ids <- names(genes_gr)
if (is.null(gene_ids) || any(gene_ids == "")) {
  gene_ids <- as.character(seq_along(genes_gr))
}
mcols(genes_gr)$gene_id <- gene_ids

promoters_gr <- promoters(genes_gr, upstream = 2000, downstream = 500)
mcols(promoters_gr)$gene_id <- mcols(genes_gr)$gene_id

gene_body_gr <- genes_gr
mcols(gene_body_gr)$gene_id <- mcols(genes_gr)$gene_id

downstream_gr <- flank(genes_gr, width = 2000, start = FALSE, both = FALSE, ignore.strand = FALSE)
mcols(downstream_gr)$gene_id <- mcols(genes_gr)$gene_id

collapse_hits <- function(dmr_gr, feature_gr, region_name) {
  h <- findOverlaps(dmr_gr, feature_gr, ignore.strand = TRUE)
  if (length(h) == 0) {
    return(data.frame())
  }
  data.frame(
    dmr_index = queryHits(h),
    gene_id = as.character(mcols(feature_gr)$gene_id[subjectHits(h)]),
    region = region_name,
    stringsAsFactors = FALSE
  )
}

annotate_dmr_table <- function(dmr_table) {
  if (nrow(dmr_table) == 0) {
    return(dmr_table)
  }
  dmr_gr <- GRanges(
    seqnames = dmr_table$chr,
    ranges = IRanges(start = dmr_table$start, end = dmr_table$end)
  )
  hit_tables <- rbind(
    collapse_hits(dmr_gr, promoters_gr, "promoter"),
    collapse_hits(dmr_gr, gene_body_gr, "gene_body"),
    collapse_hits(dmr_gr, downstream_gr, "downstream_2kb")
  )
  if (nrow(hit_tables) == 0) {
    dmr_table$primary_region <- "intergenic"
    dmr_table$associated_genes <- ""
    dmr_table$n_associated_genes <- 0L
    return(dmr_table)
  }
  hit_dt <- as.data.table(hit_tables)
  hit_dt <- unique(hit_dt[gene_id != "", ])
  if (nrow(hit_dt) == 0) {
    dmr_table$primary_region <- "intergenic"
    dmr_table$associated_genes <- ""
    dmr_table$n_associated_genes <- 0L
    return(dmr_table)
  }
  hit_dt[, region_rank := fifelse(region == "promoter", 1L, fifelse(region == "gene_body", 2L, 3L))]
  primary_region <- hit_dt[order(region_rank), .SD[1], by = dmr_index][order(dmr_index)]
  genes_per_dmr <- hit_dt[, .(associated_genes = paste(sort(unique(gene_id)), collapse = ";"), n_associated_genes = uniqueN(gene_id)), by = dmr_index]
  dmr_dt <- as.data.table(dmr_table)
  dmr_dt[, dmr_index := .I]
  dmr_dt <- merge(dmr_dt, primary_region[, .(dmr_index, primary_region = region)], by = "dmr_index", all.x = TRUE)
  dmr_dt <- merge(dmr_dt, genes_per_dmr, by = "dmr_index", all.x = TRUE)
  dmr_dt[is.na(primary_region), primary_region := "intergenic"]
  dmr_dt[is.na(associated_genes), associated_genes := ""]
  dmr_dt[is.na(n_associated_genes), n_associated_genes := 0L]
  dmr_dt$dmr_index <- NULL
  as.data.frame(dmr_dt)
}

make_gene_priority <- function(dmr_annotated) {
  if (nrow(dmr_annotated) == 0) {
    return(data.frame())
  }
  dmr_dt <- as.data.table(dmr_annotated)
  dmr_dt <- dmr_dt[associated_genes != "" & primary_region != "intergenic"]
  if (nrow(dmr_dt) == 0) {
    return(data.frame())
  }
  diff_col <- if ("diff.Methy" %in% names(dmr_dt)) "diff.Methy" else "diff"
  dmr_dt[, abs_diff := abs(get(diff_col))]
  dmr_dt[, dmr_width := end - start + 1]
  dmr_dt[, region_score := fifelse(primary_region == "promoter", 3, fifelse(primary_region == "gene_body", 2, 1))]
  dmr_dt[, dmr_priority_score := region_score * 1000 + abs_diff * 100 + log10(dmr_width + 1) * 10]
  expanded <- dmr_dt[, .(
    context,
    contrast_id,
    dmr_id,
    primary_region,
    abs_diff,
    dmr_width,
    dmr_priority_score,
    gene_id = unlist(strsplit(associated_genes, ";", fixed = TRUE))
  ), by = seq_len(nrow(dmr_dt))]
  expanded <- expanded[gene_id != ""]
  gene_rank <- expanded[, .(
    n_dmrs = uniqueN(dmr_id),
    best_region = primary_region[which.max(fifelse(primary_region == "promoter", 3, fifelse(primary_region == "gene_body", 2, 1)))][1],
    max_abs_diff = max(abs_diff, na.rm = TRUE),
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    total_priority_score = sum(dmr_priority_score, na.rm = TRUE)
  ), by = .(context, contrast_id, gene_id)]
  setorder(gene_rank, contrast_id, context, -total_priority_score, -n_dmrs, -max_abs_diff, gene_id)
  gene_rank[, priority_rank := seq_len(.N), by = .(contrast_id, context)]
  as.data.frame(gene_rank)
}

make_contrast_dir <- function(contrast_id) {
  contrast_dir <- file.path(output_dir, contrast_id)
  dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)
  contrast_dir
}

all_annotated <- list()
all_gene_rank <- list()
summary_rows <- list()

for (spec in contrast_specs) {
  message(sprintf("Контраст: %s (%s)", spec$id, spec$method))
  contrast_dir <- make_contrast_dir(spec$id)
  contrast_dmr_count <- 0L
  for (ctx in c("CpG", "CHG")) {
    message(sprintf("  контекст %s", ctx))
    params <- context_params[[ctx]]
    dmr_table <- run_dmr_for_contrast(ctx, spec, params)
    write.table(
      dmr_table,
      file = file.path(contrast_dir, paste0("DMR_", ctx, "_raw.tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    annotated <- annotate_dmr_table(dmr_table)
    write.table(
      annotated,
      file = file.path(contrast_dir, paste0("DMR_", ctx, "_annotated.tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    gene_rank <- make_gene_priority(annotated)
    write.table(
      gene_rank,
      file = file.path(contrast_dir, paste0("DMR_", ctx, "_prioritized_genes.tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    if (nrow(annotated) > 0) {
      all_annotated[[paste(spec$id, ctx, sep = "|")]] <- annotated
    }
    if (nrow(gene_rank) > 0) {
      all_gene_rank[[paste(spec$id, ctx, sep = "|")]] <- gene_rank
    }
    contrast_dmr_count <- contrast_dmr_count + nrow(dmr_table)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      contrast_id = spec$id,
      contrast_title = spec$title,
      method = spec$method,
      context = ctx,
      group1_samples = paste(spec$group1, collapse = ";"),
      group2_samples = paste(spec$group2, collapse = ";"),
      n_dmrs = nrow(dmr_table),
      stringsAsFactors = FALSE
    )
  }
  message(sprintf("  найдено DMR: %d", contrast_dmr_count))
}

summary_df <- do.call(rbind, summary_rows)
write.table(
  summary_df,
  file = file.path(output_dir, "DMR_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

if (length(all_annotated) > 0) {
  all_annotated_df <- rbindlist(all_annotated, fill = TRUE, use.names = TRUE)
  write.table(
    as.data.frame(all_annotated_df),
    file = file.path(output_dir, "DMR_all_contrasts_annotated.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

all_gene_df <- if (length(all_gene_rank) > 0) rbindlist(all_gene_rank, fill = TRUE, use.names = TRUE) else data.table()
if (nrow(all_gene_df) > 0) {
  combined <- all_gene_df[, .(
    contexts = paste(sort(unique(context)), collapse = ";"),
    contrasts = paste(sort(unique(contrast_id)), collapse = ";"),
    total_n_dmrs = sum(n_dmrs, na.rm = TRUE),
    max_abs_diff = max(max_abs_diff, na.rm = TRUE),
    mean_abs_diff = mean(mean_abs_diff, na.rm = TRUE),
    integrated_priority_score = sum(total_priority_score, na.rm = TRUE)
  ), by = gene_id]
  setorder(combined, -integrated_priority_score, -total_n_dmrs, -max_abs_diff, gene_id)
  combined[, priority_rank := seq_len(.N)]
  write.table(
    as.data.frame(combined),
    file = file.path(output_dir, "DMR_prioritized_genes_combined.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

strategy_text <- c(
  "Стратегия поиска DMR (согласована с контрастами RNA-seq DE):",
  "",
  "Образцы WGBS (факторный дизайн 2x2, по 1 образцу на ячейку):",
  "  SRR10382875: DoubleMutant, NPN",
  "  SRR10382876: DoubleMutant, PN",
  "  SRR10382877: WT, NPN",
  "  SRR10382878: WT, PN",
  "",
  "Контрасты C1-C4 (аналоги DEG, method=pairwise, n=1 vs n=1):",
  "  C1: WT NPN vs WT PN (SRR10382877 vs SRR10382878)",
  "  C2: DoubleMutant NPN vs PN (SRR10382875 vs SRR10382876)",
  "  C3: DoubleMutant vs WT в PN (SRR10382876 vs SRR10382878)",
  "  C4: DoubleMutant vs WT в NPN (SRR10382875 vs SRR10382877)",
  "  Для C1-C4: site-level z-test пропорций, diff = group1 - group2, затем callDMR.",
  "  Ограничение: без биологических реплик p-value носят exploratory-характер.",
  "",
  "Контрасты P1-P2 (pooled, method=pooled, суммирование счётчиков 2 vs 2, z-test + callDMR):",
  "  P1: pooled генотип DoubleMutant vs WT ({75,76} vs {77,78})",
  "  P2: pooled среда NPN vs PN ({75,77} vs {76,78})",
  "  Ограничение: в pooled-сравнениях смешаны уровни второго фактора.",
  "",
  "C5 (интеракция genotype:condition): не выполняется — 4 образца, 4 ячейки дизайна, нет степеней свободы.",
  "",
  "Пороги callDMR (CpG и CHG): покрытие >= 3, p <= 0.05, |delta| >= 0.05, длина >= 30 п.н.,",
  "  минимум 2 сайта, merge distance = 200, минимум 30% значимых сайтов в регионе.",
  "",
  "Границы геномных областей:",
  "  promoter: -2000..+500 от TSS; gene_body: по GFF3; downstream_2kb: 0..2000 после 3'-конца; intergenic: без пересечений.",
  "",
  "Ассоциация DMR с генами: все пересекающиеся гены; primary_region по приоритету promoter > gene_body > downstream_2kb > intergenic.",
  "",
  "Приоритизация генов: score = region_score*1000 + |diff|*100 + log10(width+1)*10; суммирование по DMR внутри contrast_id и context."
)

writeLines(strategy_text, con = file.path(output_dir, "DMR_strategy_and_rules.txt"))
writeLines("DMR-анализ завершен. Результаты сохранены в папке WGBS/DMR_results/.")
