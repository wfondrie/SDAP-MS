# ggplot theme
coolTheme <- theme_bw() +
    theme(axis.line = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          text = element_text(size = 8, color = "black"),
          panel.border = element_rect(color = "black", size = 1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 1),
          strip.background = element_rect(color = "black"))

# function to get default ggplot colors
ggColors <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

# save function
savePlot <- function(p, w, h) {
  ggsave(plot = p, filename = paste0("results/", substitute(p), ".pdf"), useDingbats = F,
         width = w, height = h, unit = "in")
}