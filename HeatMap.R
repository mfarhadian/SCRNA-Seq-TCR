

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(ggplot2)
  library(stringr)
  library(patchwork)
})

## ------------------------- INPUTS -------------------------------------------
xlsx_in <- "/data_storage/Top30_UPandDOWN.xlsx"
sheet_in <- "Union_all_rows"

desired_contrast_order <- c(
  "Motives vs. All cells",
  "Motives vs. All CD4 cells",
  "Motives vs. CD4 (TCM,TEM) + CSF",
  "Motives vs. CD4 (TCM,TEM) + CSF + CTR",
  "Motives vs. CD4 (TCM,TEM) + CSF + MS"
)
sort_contrast <- "Motives vs. CD4 (TCM,TEM) + CSF + MS"

outdir <- "/data_storage/New_HeatMap_without_star_signs_increase_title_font/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

png_out <- file.path(outdir, "Heatmap_Up_vs_Down_top30.png")
pdf_out <- file.path(outdir, "Heatmap_Up_vs_Down_top30.pdf")

## ------------------------- READ + PREP --------------------------------------
df <- read.xlsx(xlsx_in, sheet = sheet_in)
fc_col <- grep("^avg_log2FC", colnames(df), value = TRUE)[1]

df <- df %>%
  transmute(
    gene     = as.character(gene),
    contrast = as.character(contrast),
    log2FC   = as.numeric(.data[[fc_col]])
  )

df$contrast <- factor(df$contrast, levels = desired_contrast_order)

df <- df %>%
  mutate(
    contrast_id = factor(
      match(contrast, desired_contrast_order),
      levels = seq_along(desired_contrast_order),
      labels = as.character(seq_along(desired_contrast_order))
    )
  )

## ------------------------- SPLIT --------------------------------------------
df_up   <- df %>% filter(!is.na(log2FC) & log2FC >  0)
df_down <- df %>% filter(!is.na(log2FC) & log2FC <  0)

## ------------------------- ORDER ROWS ---------------------------------------
get_order_levels <- function(d, ascending = FALSE) {
  sc <- d %>%
    filter(contrast == sort_contrast) %>%
    group_by(gene) %>%
    summarise(score = mean(log2FC, na.rm = TRUE), .groups = "drop")
  if (ascending) {
    sc <- arrange(sc, score, gene)
  } else {
    sc <- arrange(sc, desc(score), gene)
  }
  sc$gene
}

up_levels   <- get_order_levels(df_up,   ascending = FALSE)
down_levels <- get_order_levels(df_down, ascending = TRUE)

df_up$gene   <- factor(df_up$gene,   levels = up_levels)
df_down$gene <- factor(df_down$gene, levels = down_levels)

## ------------------------- COMPLETE GRID (NA TILES) -------------------------
contrast_ids <- levels(df$contrast_id)

complete_grid <- function(genes, contrast_ids) {
  expand.grid(
    gene = genes,
    contrast_id = contrast_ids,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

na_tiles <- function(grid, have) {
  have_keys <- paste(have$gene, have$contrast_id)
  grid %>%
    mutate(key = paste(gene, contrast_id)) %>%
    filter(!(key %in% have_keys)) %>%
    transmute(gene, contrast_id, log2FC = NA_real_)
}

df_up_na   <- na_tiles(complete_grid(up_levels,   contrast_ids), df_up)
df_down_na <- na_tiles(complete_grid(down_levels, contrast_ids), df_down)

max_up   <- max(df_up$log2FC,   na.rm = TRUE)
min_down <- min(df_down$log2FC, na.rm = TRUE)

## ------------------------- COLORS -------------------------------------------
col_blue <- "#2166AC"
col_mid  <- "#F7F7F7"
col_red  <- "#B2182B"

## ------------------------- UPREGULATED --------------------------------------
p_up <- ggplot() +
  geom_tile(
    data = bind_rows(df_up %>% select(gene, contrast_id, log2FC), df_up_na),
    aes(x = contrast_id, y = gene, fill = log2FC)
  ) +
  scale_fill_gradient(
    low = col_mid, high = col_red,
    limits = c(0, max_up), na.value = "#EEEEEE",
    name = "Log2FC (>0)"
  ) +
  scale_x_discrete(limits = contrast_ids) +
  scale_y_discrete(position = "right", limits = rev(up_levels)) +
  ggtitle("Top 30 Upregulated Genes") +
  labs(x = "Comparison group IDs", y = NULL) +
  coord_cartesian(expand = FALSE, clip = "off") +
  theme_minimal(base_size = 28) +
  theme(
    plot.title        = element_text(face = "bold", size = 28, hjust = 0.5),
    axis.title.x      = element_text(size = 28, face = "bold"),
    axis.text.x       = element_text(face = "bold", size = 24),
    axis.text.y.right = element_text(size = 24, face = "bold"),
    axis.text.y.left  = element_blank(),
    axis.ticks.y.left = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid        = element_blank(),
    legend.title      = element_text(face = "bold", size = 17),
    legend.text       = element_text(size = 15),
    plot.margin       = margin(t = 16, r = 30, b = 30, l = 80)
  )

## ------------------------- DOWNREGULATED ------------------------------------
p_down <- ggplot() +
  geom_tile(
    data = bind_rows(df_down %>% select(gene, contrast_id, log2FC), df_down_na),
    aes(x = contrast_id, y = gene, fill = log2FC)
  ) +
  scale_fill_gradient(
    low = col_blue, high = col_mid,
    limits = c(min_down, 0), na.value = "#EEEEEE",
    name = "Log2FC (<0)"
  ) +
  scale_x_discrete(limits = contrast_ids) +
  scale_y_discrete(position = "right", limits = rev(down_levels)) +
  ggtitle("Top 30 Downregulated Genes") +
  labs(x = "Comparison group IDs", y = NULL) +
  coord_cartesian(expand = FALSE, clip = "off") +
  theme_minimal(base_size = 28) +
  theme(
    plot.title        = element_text(face = "bold", size = 28, hjust = 0.5),
    axis.title.x      = element_text(size = 28, face = "bold"),
    axis.text.x       = element_text(face = "bold", size = 20),
    axis.text.y.right = element_text(size = 21, face = "bold"),
    axis.text.y.left  = element_blank(),
    axis.ticks.y.left = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid        = element_blank(),
    legend.title      = element_text(face = "bold", size = 17),
    legend.text       = element_text(size = 15),
    plot.margin       = margin(t = 16, r = 30, b = 30, l = 80)
  )


## ------------------------- COMBINE & SAVE -----------------------------------
combined <- p_up + p_down + plot_layout(ncol = 2)

ggsave(png_out, combined, width = 23, height = 19, dpi = 1100, type = "cairo-png")
ggsave(pdf_out, combined, width = 23, height = 19, device = cairo_pdf)

message("PURE heatmaps saved:\n", png_out, "\n", pdf_out)
