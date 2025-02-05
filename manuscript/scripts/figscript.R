
library(tidyverse)

method_names <- c("pyteomics", "pymzml", "pyopenms", "pyopenms_2DPeak",
                  "mzMLb", "mzDB", "MZA", "mz5", "MZTree", "mzMD", 
                  "SQLite", "DuckDB", "Parquet")

all_timings <- read_csv("data/singlefile_times.csv") %>% 
  mutate(method=factor(method, levels=method_names)) %>%
  mutate(query=factor(query, levels=c("ms1_spec", "chrom", "rtrange",
                                      "ms2_spec", "premz", "fragmz"),
                      labels=c("MS1 scan", "Chromatogram", "RT Range",
                               "MS2 scan", "Precursor search",
                               "Fragment search")))
all_sizes <- read_csv("data/file_sizes.csv") %>%
  mutate(method=factor(method, levels=method_names))

time_boxplot <- all_timings %>%
  ggplot() +
  geom_boxplot(aes(x=method, y=time, color=method)) +
  facet_wrap(~query) +
  scale_y_log10() +
  labs(x=NULL, y="Time (seconds)", color=NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        legend.position = "none")
size_barplot <- all_sizes %>%
  filter(!method%in%c("pyteomics", "pyopenms", "pyopenms_2DPeak")) %>%
  mutate(method=as.character(method)) %>%
  mutate(method=ifelse(method=="pymzml", "mzml", method)) %>%
  mutate(method=factor(method, levels=c("mzml", method_names))) %>%
  mutate(query="File size") %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_col(aes(x=method, y=file_size/1e6, fill=method), color="black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~query) +
  labs(x=NULL, y="File size (MB)", fill="File type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), legend.position = "none")
library(patchwork)
time_boxplot + size_barplot + plot_layout(widths = c(3.5, 1))

ggsave("figures/singlefile_fig.png", width = 6.5, height=5)

