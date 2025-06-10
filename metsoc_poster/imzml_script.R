library(tidyverse)
library(RaMS)
arrow::open_dataset("imzml_data/pqds") %>%
  filter(mz%between%pmppm(421.3111, 10) | mz%between%pmppm(281.2484, 10)) %>%
  mutate(mz=paste("m/z =", round(mz, 2))) %>%
  collect() %>%
  ggplot(aes(x, y, fill=int)) + geom_raster() + facet_grid(mz~filename)


out_gp <- arrow::open_dataset("imzml_data/pqds") %>%
  filter(mz%between%pmppm(421.3111, 10) | mz%between%pmppm(281.2484, 10)) %>%
  mutate(mz=paste("m/z =", round(mz, 2))) %>%
  mutate(filename=str_remove(filename, "-centroid")) %>%
  mutate(filename=str_replace(filename, " S", "\n S")) %>%
  collect() %>%
  ggplot() +
  geom_raster(aes(x, y, fill=log10(int))) +
  ggh4x::facet_grid2(mz~filename, scales="free", independent="y") +
  scale_fill_gradientn(colors=viridis::viridis(10, direction = -1), 
                       breaks=c(1:4), labels=10^(1:4), name="Intensity") +
  scale_x_continuous(expand = expansion(), name = "DESI IMS X-coordinate") +
  scale_y_continuous(expand = expansion(), name = "DESI IMS Y-coordinate") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 0.5),
        legend.justification = c(1, 0.5),
        legend.background = element_rect(fill="#FFFFFF", color="black"))
ggsave("imzml_arrow_R_fig.png", plot = out_gp, device = "png", width = 6, height = 5, units = "in", dpi = 600)
