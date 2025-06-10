
library(tidyverse)
library(ggtext)
library(patchwork)

method_names <- c("pyteomics", "pymzml", "pyopenms", "2DPeak",
                  "mzMLb", "mzDB", "MZA", "mz5", 
                  "MZTree", "mzMD", 
                  "SQLite", "DuckDB", "Parquet")
method_colors <- c("#393B79", "#5254A3", "#6B6ECF", "#9C9EDE", 
                   "#0F8299", "#637939", "#E6550D", "#843C39", 
                   "#7B4173", "#CE6DBD", 
                   "#8C6D31", "#E7BA52", "#efddb9")
backup_colors <- c("#393B79", "#5254A3", "#6B6ECF", "#9C9EDE", 
                   "#0F8299", "#637939", "#E6550D", "#843C39", 
                   "#7B4173", "#CE6DBD", 
                   "#8C6D31", "#E7BA52", "grey20")

read_csv("timing_data/singlefile_times.csv") %>% 
  mutate(method=str_remove(method, "pyopenms_")) %>%
  mutate(method=factor(method, levels=method_names)) %>%
  mutate(query=factor(query, levels=c("ms1_spec", "chrom", "rtrange",
                                      "ms2_spec", "premz", "fragmz"),
                      labels=c("A. Full scan", "B. Chromatogram", "C. RT Range",
                               "D. MS/MS scan", "E. Precursor search",
                               "F. Fragment search"))) %>%
  ggplot() +
  geom_boxplot(aes(x=method, y=time, color=method)) +
  facet_wrap(~query) +
  scale_y_log10() +
  scale_color_manual(breaks = method_names, values=method_colors) +
  labs(x=NULL, y="Time (seconds)", color=NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        legend.position = "none",
        strip.text.x = element_text(hjust = 0),
        text=element_text(size=24))
ggsave("timing_boxplots.png", width = 9.5, height = 6, dpi=600)

read_csv("timing_data/file_sizes.csv") %>%
  mutate(method=str_remove(method, "pyopenms_")) %>%
  mutate(method=factor(method, levels=method_names)) %>%
  filter(!method%in%c("pyteomics", "pyopenms", "2DPeak")) %>%
  mutate(method=as.character(method)) %>%
  mutate(method=ifelse(method=="pymzml", "mzML", method)) %>%
  mutate(method=factor(method, levels=c("mzML", method_names))) %>%
  mutate(lab_vjust=ifelse(method%in%c("SQLite", "mzMD", "Parquet"), 1.5, -0.5)) %>%
  mutate(lab_vjust=ifelse(method%in%c("Parquet"), -2.5, lab_vjust)) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_col(aes(x=method, y=file_size/1e6, fill=method), color="black") +
  geom_label(aes(x=method, y=file_size/1e6, vjust=lab_vjust, label=method, color=method),
             size = 18/.pt, label.size=1.5) +
  scale_x_discrete(expand=expansion(mult = 0.08)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(breaks = c("mzML", method_names), values=c("#393B79", method_colors)) +
  scale_color_manual(breaks = c("mzML", method_names), values=c("#393B79", backup_colors)) +
  labs(x=NULL, y="File size (MB)", fill="File type") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(hjust = 0),
        text=element_text(size=24))
ggsave("sizing_barplot.png", width = 9.5, height = 3, dpi=600)
