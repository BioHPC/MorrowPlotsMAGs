#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(tibble)
})

BASE   <- Sys.getenv("BASE", "/lab/biohpc/Metagenomics_Huang/new_data_work")
OUTSUB <- Sys.getenv("OUTSUB", "kraken2_reads_group")

TREE_META <- file.path(BASE, "tree_all_groups", "tree_meta.tsv")
OUTDIR <- file.path(BASE, "fig_kraken_vs_mags_phylum_rel_abund")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OUTPDF <- file.path(OUTDIR, "Kraken_vs_MAGs_Phylum_relative_abundance.pdf")

TOP_N <- as.integer(Sys.getenv("TOP_N", "18"))

OTHER_LABEL <- "Other (rare phyla)"
SPECIAL <- c("Unclassified", "Classified above phylum/other ranks")

all_reports <- sort(list.files(
  path = BASE,
  pattern = "_group_reads\\.k2\\.report$",
  recursive = TRUE,
  full.names = TRUE
))
reports <- all_reports[grepl(paste0("/", OUTSUB, "/"), all_reports)]

grp_from_report <- function(p) {
  bn <- basename(p)
  sub("_group_reads\\.k2\\.report$", "", bn)
}
grp_kraken <- vapply(reports, grp_from_report, character(1))

tm <- fread(TREE_META, sep = "\t", header = TRUE, data.table = FALSE)

tm <- tm %>%
  mutate(
    group  = as.character(group),
    phylum = as.character(phylum),
    phylum = ifelse(is.na(phylum) | phylum == "" | phylum == "NA", "Unclassified", phylum)
  )

all_groups <- union(tm$group, grp_kraken)

grp_levels <- tibble(group = all_groups) %>%
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

mag_phylum <- tm %>%
  group_by(group, phylum) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(rel_abund = n / sum(n)) %>%
  ungroup() %>%
  transmute(
    group = factor(group, levels = grp_levels),
    source = "MAGs",
    phylum,
    rel_abund
  )

read_kraken_from_report <- function(report_path, grp) {
  dt <- fread(report_path, sep = "\t", header = FALSE, fill = TRUE, quote = "")
  setnames(dt, paste0("V", seq_len(ncol(dt))))

  x <- dt %>%
    transmute(
      clade     = suppressWarnings(as.numeric(V2)),
      rank_code = str_trim(as.character(V4)),
      name      = str_trim(as.character(V6))
    )

  total <- x %>%
    filter(rank_code == "R") %>%
    summarise(t = max(clade, na.rm = TRUE)) %>%
    pull(t)

  uncl_reads <- x %>%
    filter(rank_code == "U" | tolower(name) == "unclassified") %>%
    summarise(u = sum(clade, na.rm = TRUE)) %>%
    pull(u)

  ph <- x %>%
    filter(rank_code == "P") %>%
    transmute(
      group = grp,
      source = "Kraken",
      phylum = name,
      rel_abund = clade / total
    )

  sumP <- sum(ph$rel_abund, na.rm = TRUE)
  rel_uncl <- uncl_reads / total
  other_ranks <- max(0, 1 - rel_uncl - sumP)

  bind_rows(
    ph,
    tibble(group = grp, source = "Kraken", phylum = "Unclassified", rel_abund = rel_uncl),
    tibble(group = grp, source = "Kraken", phylum = "Classified above phylum/other ranks", rel_abund = other_ranks)
  ) %>%
    filter(is.finite(rel_abund), rel_abund > 0)
}

kraken_phylum <- bind_rows(lapply(seq_along(reports), function(i) {
  read_kraken_from_report(reports[i], grp_kraken[i])
}))

kraken_phylum <- kraken_phylum %>%
  mutate(group = factor(group, levels = grp_levels))

dat0 <- bind_rows(mag_phylum, kraken_phylum) %>%
  mutate(source = factor(source, levels = c("Kraken", "MAGs")))

rank_tbl <- dat0 %>%
  filter(!(phylum %in% SPECIAL)) %>%
  group_by(phylum) %>%
  summarise(mean_abund = mean(rel_abund, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund))

keep_phyla <- head(rank_tbl$phylum, TOP_N)

dat <- dat0 %>%
  mutate(
    phylum2 = case_when(
      phylum %in% SPECIAL ~ phylum,
      phylum %in% keep_phyla ~ phylum,
      TRUE ~ OTHER_LABEL
    )
  ) %>%
  group_by(group, source, phylum2) %>%
  summarise(rel_abund = sum(rel_abund, na.rm = TRUE), .groups = "drop")

mag_phyla_levels <- sort(unique(mag_phylum$phylum))
pal_mag <- setNames(hcl.colors(length(mag_phyla_levels), palette = "Dark 3"), mag_phyla_levels)

new_phyla <- setdiff(keep_phyla, names(pal_mag))
pal_new <- setNames(grDevices::hcl(seq(15, 375, length.out = length(new_phyla) + 1)[-1], c = 110, l = 55), new_phyla)

pal_misc <- c(
  OTHER_LABEL = "#BDBDBD",
  "Classified above phylum/other ranks" = "#E0E0E0",
  "Unclassified" = "#4D4D4D"
)

pal_all <- c(pal_mag, pal_new, pal_misc)
pal_all <- pal_all[!duplicated(names(pal_all))]

levels_final <- c(
  intersect(names(pal_mag), keep_phyla),
  setdiff(keep_phyla, names(pal_mag)),
  OTHER_LABEL,
  "Classified above phylum/other ranks",
  "Unclassified"
)
dat$phylum2 <- factor(dat$phylum2, levels = levels_final)

dat <- dat %>%
  mutate(
    group_chr  = as.character(group),
    group_disp = str_replace_all(group_chr, "^sod_grass_", "sg_"),
    xcat = interaction(group_disp, source, sep = " | ", lex.order = TRUE)
  ) %>%
  group_by(xcat, phylum2) %>%
  summarise(rel_abund = sum(rel_abund, na.rm = TRUE), .groups = "drop") %>%
  group_by(xcat) %>%
  mutate(rel_abund = rel_abund / sum(rel_abund, na.rm = TRUE)) %>%
  ungroup()

x_levels <- levels(dat$xcat)
x_labels <- setNames(
  sapply(x_levels, function(z) {
    parts <- str_split_fixed(z, " \\| ", 2)
    paste0(parts[1], "\n", parts[2])
  }),
  x_levels
)

p <- ggplot(dat, aes(x = xcat, y = rel_abund, fill = phylum2)) +
  geom_col(width = 0.9, colour = "grey20", linewidth = 0.1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_x_discrete(labels = x_labels) +
  scale_fill_manual(values = pal_all, drop = FALSE) +
  labs(x = "Co-assembly group", y = "Relative abundance (phylum)", fill = "Phylum") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    legend.key.height = unit(5, "pt")
  )

ggsave(OUTPDF, p, width = 18, height = 7, device = "pdf", bg = "white", limitsize = FALSE)