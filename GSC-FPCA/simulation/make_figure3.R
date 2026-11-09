library(data.table); library(ggplot2)
s <- fread("summary2.csv")
meth <- c(gsc1 = "GSC-FPCA", gsc0 = "GSC-FPCA, prior mean 0", did = "DiD", lme2 = "LME-2", mc = "MC")
x <- s[scen %in% c("none", "select") & size == "m" & sp == "B" & method %in% names(meth) & shape %in% c("any", "linear") & bin %in% sprintf("%.1f", (1:7) / 10)]
x[, t := as.numeric(bin)]; x[, method := factor(meth[method], levels = meth)]
x[, panel := factor(ifelse(scen == "none", "No selection", "Selection on scores"), levels = c("No selection", "Selection on scores"))]
cols <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"), meth)
shp <- setNames(c(16, 17, 15, 18, 4), meth)
p <- ggplot(x, aes(t, bias, color = method, shape = method, group = method)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = .4) +
  geom_line(linewidth = .7) + geom_point(size = 2.2) +
  facet_wrap(~panel, nrow = 1) +
  scale_color_manual(values = cols, name = NULL) + scale_shape_manual(values = shp, name = NULL) +
  scale_x_continuous(breaks = (1:7) / 10) +
  labs(x = "Scaled time since the last pre-exposure observation", y = "Mean error of the estimated ATT") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey92", linewidth = .3), strip.text = element_text(face = "bold"), axis.title = element_text(color = "grey25"), axis.text = element_text(color = "grey35"))
ggsave("Figure3.png", p, width = 8, height = 3.8, dpi = 300, bg = "white")
