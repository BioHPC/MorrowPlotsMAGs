#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("BASE", "/lab/biohpc/Metagenomics_Huang/new_data_work")

OUTDIR <- file.path(BASE, "fig2_like_from_checkm_coassembly")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OUTPDF <- file.path(OUTDIR, "Fig2_PanelB_C.pdf")

TREE_META <- file.path(BASE, "tree_all_groups", "tree_meta.tsv")

tm <- fread(TREE_META)

tm <- tm %>%
  mutate(
    group  = as.character(group),
    phylum = ifelse(is.na(phylum) | phylum == "" | phylum == "NA", "Unclassified", phylum)
  )

grp_levels <- tm %>%
  distinct(group) %>%
  mutate(
    ord = case_when(
      str_starts(group, "CCC") ~ 1L,
      str_starts(group, "COA") ~ 2L,
      str_starts(group, "CS")  ~ 3L,
      str_starts(group, "sod_grass_") ~ 4L,
      TRUE ~ 9L
    )
  ) %>%
  arrange(ord, group) %>%
  pull(group)

tm$group <- factor(tm$group, levels = grp_levels)

counts_phylum <- tm %>%
  group_by(group, phylum) %>%
  summarise(n_MAG = n(), .groups = "drop")

pal_phylum <- c(
  "Acidobacteriota" = "#4DBBD5",
  "Actinomycetota"  = "#00A087",
  "Chloroflexota"   = "#3C5488",
  "Desulfobacterota_B" = "#F39B7F",
  "Eisenbacteria"   = "#8491B4",
  "Gemmatimonadota" = "#91D1C2",
  "Methylomirabilota" = "#DC0000",
  "Nitrospirota"    = "#7E6148",
  "Pseudomonadota"  = "#E64B35",
  "Thermoproteota"  = "#B09C85",
  "Verrucomicrobiota" = "#A1CAF1",
  "Unclassified"    = "#999999"
)

p_phylum <- ggplot(counts_phylum, aes(x = group, y = n_MAG, fill = phylum)) +
  geom_col(width = 0.9, colour = "grey20", linewidth = 0.1) +
  scale_fill_manual(values = pal_phylum) +
  labs(
    title = "Fig2 Panel B",
    x = "Co-assembly group",
    y = "Number of MAGs",
    fill = "Phylum"
  ) +
  theme_light(12) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "top"
  )

qa_files <- list.files(BASE, pattern = "qa_summary.tsv", recursive = TRUE, full.names = TRUE)
qa_files <- qa_files[grepl("/checkm/", qa_files)]

read_checkm <- function(fp){
  lines <- readLines(fp)
  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[!grepl("^-{3,}$", lines)]
  header_idx <- grep("Bin\\s+Id", lines)[1]
  header <- strsplit(trimws(lines[header_idx]), "\\s{2,}")[[1]]
  dat <- gsub("\\s{2,}", "\t", lines[(header_idx+1):length(lines)])
  fread(text = dat, sep = "\t", header = FALSE, col.names = header)
}

all_list <- list()

for (fp in qa_files) {
  group <- sub(paste0("^", BASE, "/"), "", fp)
  group <- sub("/checkm/qa_summary.tsv$", "", group)

  qa <- tryCatch(read_checkm(fp), error = function(e) NULL)
  if (is.null(qa)) next

  df <- qa %>%
    transmute(
      group = group,
      gc = as.numeric(GC),
      length_mbp = as.numeric(`Genome size (bp)`) / 1e6
    ) %>%
    filter(!is.na(gc), !is.na(length_mbp))

  all_list[[group]] <- df
}

mag <- bind_rows(all_list)
mag$group <- factor(mag$group, levels = grp_levels)

p_gc <- ggplot(mag, aes(x = group, y = gc)) +
  geom_boxplot(fill = "grey85", colour = "black",
               width = 0.7, outlier.size = 0.6, linewidth = 0.35) +
  labs(y = "GC content (%)", x = NULL, title = "Fig2 Panel C") +
  theme_light(12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )

p_len <- ggplot(mag, aes(x = group, y = length_mbp)) +
  geom_boxplot(fill = "grey85", colour = "black",
               width = 0.7, outlier.size = 0.6, linewidth = 0.35) +
  labs(y = "Length (Mbp)", x = "Co-assembly group") +
  theme_light(12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )

fig <- p_phylum / p_gc / p_len +
  plot_layout(heights = c(1.2, 1, 1))

ggsave(OUTPDF, fig, width = 14, height = 12, device = "pdf", bg = "white")