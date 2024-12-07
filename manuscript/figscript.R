

library(tidyverse)


all_timings <- c("chrom", "spec", "rtrange") %>%
  map(function(metric){
  read_csv(paste0(metric, "_timing_data.csv")) %>%
    pivot_longer(everything()) %>%
    mutate(metric=metric)
}) %>%
  bind_rows() %>%
  mutate(category=case_when(
    str_detect(name, "pyteomics|pyopenms|pymzml")~"mzML",
    str_detect(name, "SQLite|DuckDB")~"Database",
    TRUE~name
  )) %>%
  mutate(name=ifelse(name=="MZA (package)", "mzapy", name)) %>%
  mutate(name=fct_inorder(name)) %>%
  mutate(category=fct_inorder(category)) %>%
  mutate(metric=factor(metric, levels=c("spec", "chrom", "rtrange"),
                       labels=c("Spectrum extraction",
                                "Chromatogram extraction",
                                "Retention time range")))

all_timings %>%
  filter(category=="mzML") %>%
  filter(!str_detect(name, "idx")) %>%
  filter(metric=="Spectrum extraction") %>%
  ggplot(aes(x=name, y=value)) +
  geom_boxplot(color="black") +
  geom_jitter(width = 0.1, alpha=0.5) +
  geom_hline(yintercept = 0) +
  facet_wrap(~metric, ncol=1, scales = "free_y") +
  theme_bw() +
  theme(text=element_text(size=24)) +
  labs(x=NULL, y="Query time (seconds)")

all_timings %>%
  filter(category=="mzML") %>%
  filter(!str_detect(name, "idx")) %>%
  mutate(name=factor(name, levels=c(
    "pyteomics", "pyopenms", "pyopenms_2d", "pymzml"
  ), labels= c(
    "pyteomics", "pyopenms", "pyopenms\n(2D peak method)", "pymzml"
  ))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_boxplot(color="black", outliers = FALSE) +
  geom_jitter(width = 0.1, alpha=0.5) +
  geom_hline(yintercept = 0) +
  facet_wrap(~metric, ncol=1, scales = "free_y") +
  theme_bw() +
  theme(text=element_text(size=24)) +
  labs(x=NULL, y="Query time (seconds)")

all_timings %>%
  filter(category=="mzML") %>%
  mutate(mzml_type=ifelse(str_detect(name, "idx"), "indexed", "default")) %>%
  mutate(name=str_remove(name, "\r\n\\(idx mzML\\)")) %>%
  mutate(name=factor(name, levels=c(
    "pyteomics", "pyopenms", "pyopenms_2d", "pymzml"
  ), labels= c(
    "pyteomics", "pyopenms", "pyopenms\n(2D peak method)", "pymzml"
  ))) %>%
  ggplot(aes(x=name, y=value, color=mzml_type)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(alpha=0.5, position=position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0) +
  facet_wrap(~metric, ncol=1, scales = "free_y") +
  theme_bw() +
  theme(text=element_text(size=24), legend.position = "none") +
  labs(x=NULL, y="Query time (seconds)") +
  scale_color_manual(breaks=c("default", "indexed"), values = c("black", "#a41118"))


all_timings %>%
  filter(category=="mzML") %>%
  filter(metric=="Spectrum extraction") %>%
  mutate(mzml_type=ifelse(str_detect(name, "idx"), "indexed", "default")) %>%
  mutate(name=str_remove(name, "\r\n\\(idx mzML\\)")) %>%
  mutate(name=factor(name, levels=c(
    "pyteomics", "pyopenms", "pyopenms_2d", "pymzml"
  ), labels= c(
    "pyteomics", "pyopenms", "pyopenms\n(2D peak method)", "pymzml"
  ))) %>%
  ggplot(aes(x=name, y=value, color=mzml_type)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(alpha=0.5, position=position_dodge(width = 0.75)) +
  facet_wrap(~metric, ncol=1, scales = "free_y") +
  theme_bw() +
  theme(text=element_text(size=24), legend.position = "none") +
  labs(x=NULL, y="Query time (seconds)") +
  scale_color_manual(breaks=c("default", "indexed"), values = c("black", "#a41118")) +
  scale_y_log10()



mzml_baselines <- all_timings %>%
  filter(category=="mzML") %>%
  group_by(name, metric) %>%
  summarise(mean_val=mean(value)) %>%
  ungroup() %>%
  filter(mean_val==min(mean_val), .by=metric)



all_timings %>%
  filter(!category%in%c("mzML", "Database")) %>%
  ggplot(aes(x=name, y=value)) +
  geom_hline(aes(yintercept=mean_val), data=mzml_baselines, color="#a41118", lwd=1.2) +
  geom_label(aes(y=mean_val), x=0, label = "Best mzML", data=mzml_baselines, 
             hjust=0.5, color="#a41118", label.r = unit(0, "in")) +
  geom_boxplot(outliers = FALSE, color="white") +
  geom_jitter(alpha=0.5, position=position_dodge(width = 0.75), color="white") +
  facet_wrap(~metric, ncol=1, scales = "free_y") +
  theme_bw() +
  theme(text=element_text(size=24), legend.position = "none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(x=NULL, y="Query time (seconds)") +
  scale_y_log10() +
  scale_x_discrete(expand = expansion(add=c(1.7, 0)))


all_timings %>%
  filter(!category%in%c("mzML", "Database")) %>%
  ggplot(aes(x=name, y=value)) +
  geom_hline(aes(yintercept=mean_val), data=mzml_baselines, color="#a41118", lwd=1.2) +
  geom_label(aes(y=mean_val), x=0, label = "Best mzML", data=mzml_baselines, 
             hjust=0.5, color="#a41118", label.r = unit(0, "in")) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(alpha=0.5, position=position_dodge(width = 0.75)) +
  facet_wrap(~metric, ncol=1, scales = "free_y") +
  theme_bw() +
  theme(text=element_text(size=24), legend.position = "none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(x=NULL, y="Query time (seconds)") +
  scale_y_log10() +
  scale_x_discrete(expand = expansion(add=c(1.7, 0)))



nondb_baselines <- all_timings %>%
  filter(category!="Database") %>%
  group_by(name, metric) %>%
  summarise(mean_val=mean(value)) %>%
  ungroup() %>%
  filter(mean_val==min(mean_val), .by=metric)



all_timings %>%
  filter(category=="Database") %>%
  ggplot(aes(x=name, y=value)) +
  geom_boxplot(outliers = FALSE, color="white") +
  geom_jitter(alpha=0.5, position=position_dodge(width = 0.75), color="white") +
  geom_hline(aes(yintercept=mean_val), data=nondb_baselines, color="#028e34", lwd=1.2) +
  geom_label(aes(y=mean_val), x=0, label = "Best non-DB", data=nondb_baselines, 
             hjust=0.5, color="#028e34", label.r = unit(0, "in"), vjust=c(0, 1, 1)) +
  facet_wrap(~metric, ncol=1, scales = "free_y") +
  theme_bw() +
  theme(text=element_text(size=24), legend.position = "none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(x=NULL, y="Query time (seconds)") +
  scale_y_log10() +
  scale_x_discrete(expand = expansion(add=c(1.7, 0)))
all_timings %>%
  filter(category=="Database") %>%
  ggplot(aes(x=name, y=value)) +
  geom_hline(aes(yintercept=mean_val), data=nondb_baselines, color="#028e34", lwd=1.2) +
  geom_label(aes(y=mean_val), x=0, label = "Best non-DB", data=nondb_baselines, 
             hjust=0.5, color="#028e34", label.r = unit(0, "in"), vjust=c(0, 1, 1)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(alpha=0.5, position=position_dodge(width = 0.75)) +
  facet_wrap(~metric, ncol=1, scales = "free_y") +
  theme_bw() +
  theme(text=element_text(size=24), legend.position = "none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(x=NULL, y="Query time (seconds)") +
  scale_y_log10() +
  scale_x_discrete(expand = expansion(add=c(1.7, 0)))



all_timings %>%
  mutate(name=str_remove(name, "(\r)?\n.*")) %>%
  arrange(mean(log10(value)), .by=name) %>%
  mutate(name=fct_inorder(name)) %>%
  ggplot() +
  geom_boxplot(aes(x=name, y=value, color=category)) +
  facet_wrap(~metric, ncol=1, scales = "free_y") +
  scale_y_log10() +
  labs(x=NULL, y="Query time (seconds)") +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

file_sizes <- sapply(list.files("../demo_data", full.names = TRUE), file.size)

plot_cats <- all_timings %>%
  distinct(category, name) %>%
  mutate(category=as.character(category)) %>%
  mutate(name=as.character(name)) %>%
  mutate(category=ifelse(category=="Database", name, category)) %>% 
  mutate(category=str_remove(category, "\r\n.*"))
as.data.frame(file_sizes) %>%
  rownames_to_column("filename") %>%
  mutate(filename=str_remove(filename, "../demo_data/180205_Poo_TruePoo_Full1")) %>%
  filter(!filename%in%c(".json", ".raw.mzDB", ".mzTree", "../demo_data/README.md", "_idx.mzML", ".aird")) %>%
  mutate(category=c("DuckDB", "mz5 (*)", "MZA", "mzDB", "mzMD (+)", "mzML", "mzMLb", "MZTree (+)", "SQLite")) %>%
  left_join(plot_cats) %>%
  select(name, file_sizes) %>%
  mutate(file_sizes=file_sizes/1e6) %>%
  left_join(all_timings) %>%
  ggplot() +
  geom_point(aes(x=file_sizes, y=value, fill=category), size=4, color="black", pch=21) +
  facet_wrap(~metric, ncol=1, scales = "free_y") +
  scale_y_log10() +
  labs(x="File size on disk (MB)", y="Query time (seconds)", fill="File type") +
  lims(x=c(0, NA)) +
  theme_bw() +
  theme(text=element_text(size=24))
