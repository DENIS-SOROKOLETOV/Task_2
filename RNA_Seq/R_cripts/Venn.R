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
    "Эффект среды в ΔΔ: NPN vs PN",
    "Эффект генотипа в PN: ΔΔ vs WT",
    "Эффект генотипа в NPN: ΔΔ vs WT",
    "Эффект взаимодействия genotype:condition"
  ),
  stringsAsFactors = FALSE
)

deg_files <- list.files(de_dir, pattern = "_DEG\\.tsv$", recursive = TRUE, full.names = TRUE)
deg_lookup <- setNames(deg_files, sub("_DEG\\.tsv$", "", basename(deg_files)))

read_deg_set <- function(contrast_id) {
  path <- deg_lookup[[contrast_id]]
  if (is.null(path) || !nzchar(path)) {
    path <- NA_character_
  }
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

safe_label <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

analyze_pair <- function(id1, id2, file_name1, file_name2, label1, label2, out_prefix, plot_title) {
  set1 <- deg_sets[[id1]]
  set2 <- deg_sets[[id2]]
  common <- intersect(set1, set2)
  only_1 <- setdiff(set1, set2)
  only_2 <- setdiff(set2, set1)

  write.table(data.frame(Geneid = common, stringsAsFactors = FALSE),
              file.path(output_dir, paste0(out_prefix, "_intersection.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(Geneid = only_1, stringsAsFactors = FALSE),
              file.path(output_dir, paste0(out_prefix, "_only_", safe_label(file_name1), ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(Geneid = only_2, stringsAsFactors = FALSE),
              file.path(output_dir, paste0(out_prefix, "_only_", safe_label(file_name2), ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  venn_file <- file.path(output_dir, paste0(out_prefix, ".png"))
  png(filename = venn_file, width = 1800, height = 1600, res = 220)
  venn_grob <- draw.pairwise.venn(
    area1 = length(set1),
    area2 = length(set2),
    cross.area = length(common),
    category = c(label1, label2),
    fill = c("#4DAF4A", "#377EB8"),
    alpha = c(0.55, 0.55),
    cex = 1.5,
    cat.cex = 1.0,
    cat.pos = c(-45, 45),
    cat.dist = c(0.14, 0.14),
    scaled = FALSE
  )
  grid.newpage()
  grid.draw(venn_grob)
  grid.text(plot_title, y = unit(0.965, "npc"), gp = gpar(fontsize = 17, fontface = "bold"))
  dev.off()
}

analyze_pair(
  id1 = "C1_WT_NPN_vs_WT_PN",
  id2 = "C2_DoubleMutant_NPN_vs_DoubleMutant_PN",
  file_name1 = "WT_NPN_vs_PN",
  file_name2 = "DoubleMutant_NPN_vs_PN",
  label1 = "WT NPN vs WT PN",
  label2 = "ΔΔ NPN vs ΔΔ PN",
  out_prefix = "venn_C1_WT_NPN_vs_WT_PN_C2_DoubleMutant_NPN_vs_DoubleMutant_PN",
  plot_title = "Ответ на среду: WT и ΔΔ (NPN vs PN)"
)
analyze_pair(
  id1 = "C3_DoubleMutant_PN_vs_WT_PN",
  id2 = "C4_DoubleMutant_NPN_vs_WT_NPN",
  file_name1 = "DoubleMutant_vs_WT_in_PN",
  file_name2 = "DoubleMutant_vs_WT_in_NPN",
  label1 = "ΔΔ PN vs WT PN",
  label2 = "ΔΔ NPN vs WT NPN",
  out_prefix = "venn_C3_DoubleMutant_PN_vs_WT_PN_C4_DoubleMutant_NPN_vs_WT_NPN",
  plot_title = "Эффект генотипа: ΔΔ vs WT (PN и NPN)"
)

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
  category = c("WT NPN vs WT PN", "ΔΔ NPN vs ΔΔ PN", "ΔΔ PN vs WT PN", "ΔΔ NPN vs WT NPN"),
  fill = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
  alpha = c(0.5, 0.5, 0.5, 0.5),
  cex = 1.1,
  cat.cex = 1.0
)
grid.newpage()
grid.draw(venn_quad)
grid.text("Общий и уникальный вклад контрастов", y = unit(0.97, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))
dev.off()

writeLines("Venn-анализ завершен. Результаты сохранены в папке Venn_results.")
