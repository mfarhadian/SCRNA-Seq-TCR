
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(grid) 
})

#==============================================================================
# INPUT DATA
#==============================================================================
df <- tibble::tribble(
  ~Motife,   ~Blood, ~CSF,
  "Motif 1",     1,   63,
  "Motif 2",     3,   38,
  "Motif 3",     1,   17
)

#==============================================================================
# SETTINGS
#==============================================================================
blood_color <- "#F8766D"
csf_color   <- "cornflowerblue"

out_dir <- "~/data_storage/Motif_Donut/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

FONT_FAMILY    <- "Arial"
TEXT_SIZE_PT   <- 35
GEOM_TEXT_SIZE <- TEXT_SIZE_PT / ggplot2::.pt

# Donut geometry (UNCHANGED)
ring_x    <- 2
bar_width <- 0.95
xlim_min  <- 0.9
xlim_max  <- 2.6

TITLE_VJUST  <- -0.9
STACK_GAP_PT <- 2

#==============================================================================
# CSF LABEL CONTROL
#==============================================================================
CSF_LABEL_ANGLE_FRAC <- 0.50

#==============================================================================
# LEGEND SIZE CONTROLS (NEW – legend ONLY)
#==============================================================================
LEGEND_TEXT_SIZE_PT <- 48
LEGEND_KEY_W        <- unit(1.9, "cm")
LEGEND_KEY_H        <- unit(0.85, "cm")
LEGEND_SPACING_Y    <- unit(0.45, "cm")

#==============================================================================
# DONUT FUNCTION
#==============================================================================
make_donut <- function(motif_name, blood_n, csf_n, show_legend = FALSE){
  
  dat <- tibble(
    Compartment = c("CSF","Blood"),
    Count       = c(csf_n, blood_n)
  ) %>%
    mutate(Compartment = factor(Compartment, levels = c("CSF","Blood"))) %>%
    arrange(desc(Compartment)) %>%
    mutate(ypos = cumsum(Count) - Count/2)
  
  total <- sum(dat$Count)
  
  label_dat <- dat %>%
    mutate(
      y_label = ifelse(
        Compartment == "CSF",
        total * CSF_LABEL_ANGLE_FRAC,
        ypos
      )
    )
  
  ggplot(dat, aes(x = ring_x, y = Count, fill = Compartment)) +
    geom_col(width = bar_width, color = "white", linewidth = 1.2) +
    coord_polar(theta = "y") +
    xlim(xlim_min, xlim_max) +
    
    geom_text(
      data = label_dat %>% filter(Count > 0),
      aes(x = ring_x, y = y_label, label = Count),
      inherit.aes = FALSE,
      size = GEOM_TEXT_SIZE,
      family = FONT_FAMILY,
      fontface = "bold",
      color = "black"
    ) +
    
    scale_fill_manual(
      values = c("CSF" = csf_color, "Blood" = blood_color),
      limits = c("CSF","Blood"),
      guide = guide_legend(
        ncol = 1,
        byrow = TRUE,
        keywidth  = LEGEND_KEY_W,
        keyheight = LEGEND_KEY_H
      )
    ) +
    
    ggtitle(motif_name) +
    theme_void() +
    theme(
      text = element_text(family = FONT_FAMILY, face = "bold", size = TEXT_SIZE_PT),
      plot.title = element_text(hjust = 0.5, vjust = TITLE_VJUST),
      plot.margin = margin(t = 0, r = 0, b = STACK_GAP_PT, l = 0, unit = "pt"),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_blank(),
      legend.text  = element_text(size = LEGEND_TEXT_SIZE_PT),
      legend.key.width  = LEGEND_KEY_W,
      legend.key.height = LEGEND_KEY_H,
      legend.spacing.y  = LEGEND_SPACING_Y
    )
}

#==============================================================================
# BUILD PLOTS
#==============================================================================
plots <- vector("list", nrow(df))
for(i in seq_len(nrow(df))){
  plots[[i]] <- make_donut(
    motif_name  = df$Motife[i],
    blood_n     = df$Blood[i],
    csf_n       = df$CSF[i],
    show_legend = (i == 1)
  )
}

#==============================================================================
# COMBINE + FORCE LEGEND BOTTOM
#==============================================================================
final_plot <-
  wrap_plots(plots, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.text = element_text(size = LEGEND_TEXT_SIZE_PT),
    legend.key.width  = LEGEND_KEY_W,
    legend.key.height = LEGEND_KEY_H,
    legend.spacing.y  = LEGEND_SPACING_Y
  )

#==============================================================================
# SAVE
#==============================================================================
ggsave(
  filename = file.path(out_dir, "Motif_donuts_FINAL_legendONLY_enlarged.png"),
  plot = final_plot,
  width = 6.5,
  height = 18,
  dpi = 500,
  bg = "white"
)

