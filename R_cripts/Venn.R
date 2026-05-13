required_packages <- c("VennDiagram", "grid")
missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  stop(sprintf("Не установлены пакеты: %s", paste(missing_required, collapse = ", ")))
}

suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
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
output_dir <- file.path(project_root, "Venn_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

contrast_info <- data.frame(
  contrast_id = c(
    "C1_WT_NPN_vs_WT_PN",
    "C2_DoubleMutant_NPN_vs_DoubleMutant_PN",
    "C3_DoubleMutant_PN_vs_WT_PN",
    "C4_DoubleMutant_NPN_vs_WT_NPN",
    "C5_Interaction_effect"
  ),
  title = c(
    "Эффект среды в WT: NPN vs PN",
    "Эффект среды в DoubleMutant: NPN vs PN",
    "Эффект генотипа в PN: DoubleMutant vs WT",
    "Эффект генотипа в NPN: DoubleMutant vs WT",
    "Эффект взаимодействия genotype:condition"
  ),
  stringsAsFactors = FALSE
)

read_deg_set <- function(contrast_id) {
  path <- file.path(de_dir, paste0(contrast_id, "_DEG.tsv"))
  if (!file.exists(path)) {
    stop(sprintf("Не найден файл DEG: %s", path))
  }
  tab <- read.delim(path, sep = "\t", header = TRUE, check.names = FALSE)
  if (!("Geneid" %in% colnames(tab))) {
    stop(sprintf("В файле %s отсутствует колонка Geneid.", path))
  }
  unique(tab$Geneid[!is.na(tab$Geneid) & tab$Geneid != ""])
}

deg_sets <- setNames(lapply(contrast_info$contrast_id, read_deg_set), contrast_info$contrast_id)
set_sizes <- data.frame(
  contrast_id = names(deg_sets),
  n_deg = vapply(deg_sets, length, integer(1)),
  stringsAsFactors = FALSE
)
set_sizes <- merge(set_sizes, contrast_info, by = "contrast_id", all.x = TRUE, sort = FALSE)
set_sizes <- set_sizes[, c("contrast_id", "title", "n_deg")]
write.table(set_sizes, file.path(output_dir, "deg_set_sizes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

safe_label <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

analyze_pair <- function(id1, id2, name1, name2, out_prefix) {
  set1 <- deg_sets[[id1]]
  set2 <- deg_sets[[id2]]
  common <- intersect(set1, set2)
  only_1 <- setdiff(set1, set2)
  only_2 <- setdiff(set2, set1)
  union_set <- union(set1, set2)
  jaccard <- if (length(union_set) == 0) 0 else length(common) / length(union_set)

  write.table(data.frame(Geneid = common, stringsAsFactors = FALSE),
              file.path(output_dir, paste0(out_prefix, "_intersection.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(Geneid = only_1, stringsAsFactors = FALSE),
              file.path(output_dir, paste0(out_prefix, "_only_", safe_label(name1), ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(Geneid = only_2, stringsAsFactors = FALSE),
              file.path(output_dir, paste0(out_prefix, "_only_", safe_label(name2), ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  pair_stats <- data.frame(
    comparison = out_prefix,
    set1 = name1,
    set2 = name2,
    n_set1 = length(set1),
    n_set2 = length(set2),
    n_intersection = length(common),
    n_union = length(union_set),
    jaccard_index = round(jaccard, 4),
    stringsAsFactors = FALSE
  )

  venn_file <- file.path(output_dir, paste0(out_prefix, ".png"))
  png(filename = venn_file, width = 1800, height = 1600, res = 220)
  venn_grob <- draw.pairwise.venn(
    area1 = length(set1),
    area2 = length(set2),
    cross.area = length(common),
    category = c(name1, name2),
    fill = c("#4DAF4A", "#377EB8"),
    alpha = c(0.55, 0.55),
    cex = 1.5,
    cat.cex = 1.2,
    cat.pos = c(-20, 20),
    cat.dist = c(0.045, 0.045),
    scaled = FALSE
  )
  grid.newpage()
  grid.draw(venn_grob)
  dev.off()

  pair_stats
}

pair_stats_all <- list()
pair_stats_all[["nutrient_response_by_genotype"]] <- analyze_pair(
  id1 = "C1_WT_NPN_vs_WT_PN",
  id2 = "C2_DoubleMutant_NPN_vs_DoubleMutant_PN",
  name1 = "WT_NPN_vs_PN",
  name2 = "DoubleMutant_NPN_vs_PN",
  out_prefix = "venn_C1_WT_NPN_vs_WT_PN_C2_DoubleMutant_NPN_vs_DoubleMutant_PN"
)
pair_stats_all[["genotype_response_by_environment"]] <- analyze_pair(
  id1 = "C3_DoubleMutant_PN_vs_WT_PN",
  id2 = "C4_DoubleMutant_NPN_vs_WT_NPN",
  name1 = "DoubleMutant_vs_WT_in_PN",
  name2 = "DoubleMutant_vs_WT_in_NPN",
  out_prefix = "venn_C3_DoubleMutant_PN_vs_WT_PN_C4_DoubleMutant_NPN_vs_WT_NPN"
)
pair_stats_table <- do.call(rbind, pair_stats_all)
write.table(pair_stats_table, file.path(output_dir, "pairwise_overlap_stats.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

core_ids <- c(
  "C1_WT_NPN_vs_WT_PN",
  "C2_DoubleMutant_NPN_vs_DoubleMutant_PN",
  "C3_DoubleMutant_PN_vs_WT_PN",
  "C4_DoubleMutant_NPN_vs_WT_NPN"
)
core_sets <- deg_sets[core_ids]

common_all_core <- Reduce(intersect, core_sets)
write.table(data.frame(Geneid = common_all_core, stringsAsFactors = FALSE),
            file.path(output_dir, "common_genes_all_core_contrasts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

unique_by_core <- do.call(rbind, lapply(core_ids, function(id) {
  other <- core_sets[setdiff(core_ids, id)]
  only_this <- setdiff(core_sets[[id]], unique(unlist(other)))
  data.frame(
    contrast_id = id,
    Geneid = only_this,
    stringsAsFactors = FALSE
  )
}))
if (nrow(unique_by_core) > 0) {
  unique_by_core <- merge(unique_by_core, contrast_info, by = "contrast_id", all.x = TRUE, sort = FALSE)
  unique_by_core <- unique_by_core[, c("contrast_id", "title", "Geneid")]
}
write.table(unique_by_core, file.path(output_dir, "unique_genes_by_core_contrast.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

overlap_with_interaction <- do.call(rbind, lapply(core_ids, function(id) {
  common <- intersect(core_sets[[id]], deg_sets[["C5_Interaction_effect"]])
  data.frame(
    contrast_id = id,
    n_overlap_with_C5 = length(common),
    stringsAsFactors = FALSE
  )
}))
overlap_with_interaction <- merge(overlap_with_interaction, contrast_info, by = "contrast_id", all.x = TRUE, sort = FALSE)
overlap_with_interaction <- overlap_with_interaction[, c("contrast_id", "title", "n_overlap_with_C5")]
write.table(overlap_with_interaction, file.path(output_dir, "overlap_with_interaction_C5.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

quad_file <- file.path(output_dir, "venn_core_4contrasts.png")
png(filename = quad_file, width = 2100, height = 1900, res = 220)
venn_quad <- draw.quad.venn(
  area1 = length(core_sets[[1]]),
  area2 = length(core_sets[[2]]),
  area3 = length(core_sets[[3]]),
  area4 = length(core_sets[[4]]),
  n12 = length(intersect(core_sets[[1]], core_sets[[2]])),
  n13 = length(intersect(core_sets[[1]], core_sets[[3]])),
  n14 = length(intersect(core_sets[[1]], core_sets[[4]])),
  n23 = length(intersect(core_sets[[2]], core_sets[[3]])),
  n24 = length(intersect(core_sets[[2]], core_sets[[4]])),
  n34 = length(intersect(core_sets[[3]], core_sets[[4]])),
  n123 = length(Reduce(intersect, core_sets[c(1, 2, 3)])),
  n124 = length(Reduce(intersect, core_sets[c(1, 2, 4)])),
  n134 = length(Reduce(intersect, core_sets[c(1, 3, 4)])),
  n234 = length(Reduce(intersect, core_sets[c(2, 3, 4)])),
  n1234 = length(common_all_core),
  category = c("C1_WT_env", "C2_DM_env", "C3_genotype_PN", "C4_genotype_NPN"),
  fill = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
  alpha = c(0.5, 0.5, 0.5, 0.5),
  cex = 1.1,
  cat.cex = 1.0
)
grid.newpage()
grid.draw(venn_quad)
dev.off()

report_lines <- c(
  "Venn-анализ DEG завершен.",
  "Контекст варианта 4: дизайн 2x2 (генотип WT/DoubleMutant и среда PN/NPN).",
  "Биологический вопрос 1: какие DEG отвечают на смену среды в обоих генотипах и какие специфичны для WT или DoubleMutant (сравнение C1 и C2).",
  "Биологический вопрос 2: какие DEG отражают эффект генотипа в разных средах и какие специфичны для PN или NPN (сравнение C3 и C4).",
  sprintf("Всего контрастов с DEG: %d", nrow(set_sizes)),
  sprintf("Общих DEG для C1-C4: %d", length(common_all_core)),
  "Основные таблицы и изображения сохранены в папке Venn_results."
)
writeLines(report_lines, file.path(output_dir, "Venn_report.txt"))

writeLines("Venn-анализ завершен. Результаты сохранены в папке Venn_results.")
