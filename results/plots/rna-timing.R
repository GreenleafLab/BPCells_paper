library(tidyverse)
library(patchwork)

################################################################################
# Data loading 
################################################################################

repo_root <- "/home/bparks/Sync/BPCells_paper/final_upload_version/"

results_dir <- file.path(repo_root, "/results/data_tables")
plot_dir <- file.path(repo_root, "/results/plots/rna-timing")
dir.create(plot_dir, showWarnings = FALSE)

rna_datasets <- read_tsv(file.path(results_dir, "datasets/rna.tsv"))
pca_benchmark_raw <- read_tsv(file.path(results_dir, "rna-timing/pca-benchmark/performance.tsv"),
                          col_types = "cccicdddd")

pca_benchmark <- pca_benchmark_raw |>
  mutate(tool_split=if_else(tool=="scanpy" & dataset=="4m_fetal", "scanpy_split", tool))

pca_benchmark_file_sizes <- read_tsv(file.path(results_dir, "rna-timing/pca-benchmark/staged_file_size.tsv"), col_types="cdccdd")

pca_accuracy <- read_tsv(file.path(results_dir, "rna-timing/pca-benchmark/accuracy.tsv.gz"),
                          col_types = "cccciidc")

matrix_transpose <- read_tsv(file.path(results_dir, "rna-timing/matrix-transpose.tsv"))

marker_performance <- read_tsv(file.path(results_dir, "rna-timing/marker-genes/performance.tsv"))
marker_accuracy <- read_tsv(file.path(results_dir, "rna-timing/marker-genes/accuracy.tsv"))

cell_metadata_size <- tribble(
  ~dataset, ~cell_metadata_bytes,
  "11m_jax", 477138247,
  "130k_thymus_atlas", 3734668,
  "1m_neurons", 26375993,
  "20k_pbmc", 452903,
  "22m_pansci", 1293335242,
  "2m_perturbseq", 40973855,
  "480k_tabula_sapiens", 21757014,
  "4m_fetal", 155771190,
  "500k_drugscreen", 10717711,    
)

################################################################################
# Plot styling helpers
################################################################################

set1 <- RColorBrewer::brewer.pal(4, name="Set1")
paired <- RColorBrewer::brewer.pal(9, "Paired")
tool_colors <- c("BPCells"=set1[1], "ArchR"=set1[2], "10x"=set1[3], 
                 "dgCMatrix"=set1[2], "SciPy"=set1[3], "Presto"=set1[2],
                 "Scanpy"=set1[3], "Seurat"=set1[2], "DelayedArray"=set1[4],
                 "Scanpy + Dask"=set1[3], "BPCells randomized"=paired[5],
                 "Presto"=set1[2], "Scanpy Uncorrected"=paired[3])



normalize_name <- c(
  bpcells_stage_float="BPCells", seurat="Seurat", delayedarray="DelayedArray", scanpy="Scanpy",
  bpcells_randomized="BPCells randomized", scanpy_dask = "Scanpy + Dask", bpcells="BPCells",
  scipy="SciPy", dgCMatrix="dgCMatrix", presto="Presto"
)

cell_axis <- rna_datasets |>
  filter(dataset != "1m_neurons") %>%
  mutate(label=str_split_fixed(dataset, "_", 2)[,1]) %>%
  mutate(label=case_when(dataset=="4m_fetal" ~ "1m\n4m", dataset=="500k_drugscreen" ~ "", TRUE ~ label)) %>%
  mutate(label=str_replace_all(label, "m", "M")) %>%
  arrange(nonzero_entries) %>%
  {sec_axis(transform="identity", breaks=.$nonzero_entries, labels=.$label, guide=guide_axis(angle=90))}

# Slight concern of B being confused as 10^12 for European readers, but I find B much clearer than G for US audiences 
label_si_custom <- function(x) {
  ret <- scales::label_number(scale_cut=scales::cut_short_scale())(x)
  return(ret)
}
  
style_pca <- list(
  stat_summary(fun="mean", geom="line", key_glyph=draw_key_rect),
  stat_summary(fun.min=min, fun.max=max, fun="mean", geom="pointrange", size=0.25),
  scale_x_continuous(transform="log10", guide="axis_logticks", sec.axis=cell_axis, labels = label_si_custom),
  scale_y_continuous(transform="log10", guide="axis_logticks", labels = label_si_custom),
  scale_color_manual(values=tool_colors),
  scale_fill_manual(values=tool_colors),
  guides(color="none", fill="none"),
  labs(x="Nonzero entries", fill="Tool"),
  theme_classic(base_size = 14, base_family="Arial"),
  theme(panel.border = element_rect(color="black", linewidth = rel(1), fill=NA))
)

redundant_cell_axis <- rna_datasets |>
  mutate(label=str_split_fixed(dataset, "_", 2)[,1]) %>%
  mutate(label=case_when(dataset=="480k_tabula_sapiens" ~ "480k\n500k", dataset=="500k_drugscreen" ~ "", TRUE ~ label)) %>%
  mutate(label=str_replace_all(label, "m", "M")) %>%
  arrange(cells) %>%
  {sec_axis(transform="identity", breaks=.$cells, labels=.$label, guide=guide_axis(angle=90))}


style_pca_by_cell <- list(
  stat_summary(fun="mean", geom="line", key_glyph=draw_key_rect),
  stat_summary(fun.min=min, fun.max=max, fun="mean", geom="pointrange", size=0.25),
  scale_x_continuous(transform="log10", guide="axis_logticks", sec.axis=redundant_cell_axis, labels = label_si_custom),
  scale_y_continuous(transform="log10", guide="axis_logticks", labels=label_si_custom),
  scale_color_manual(values=tool_colors),
  scale_fill_manual(values=tool_colors),
  guides(color="none", fill="none"),
  labs(x="Cells", fill="Tool"),
  theme_classic(base_size = 14, base_family="Arial"),
  theme(panel.border = element_rect(color="black", linewidth = rel(1), fill=NA))
)

################################################################################
# Standard PCA memory + timing
################################################################################

# Fig 1b
pca_memory_combined_by_cell <- pca_benchmark |>
  filter(threads == 1, tool %in% c("bpcells_stage_float", "delayedarray", "scanpy", "seurat")) |>
  group_by(tool, tool_split, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  inner_join(rna_datasets, by="dataset") |>
  ggplot(aes(cells, max_rss, color=normalize_name[tool], fill=normalize_name[tool], group=tool_split)) +
  style_pca_by_cell +
  labs(y="Memory (GB)", subtitle="Normalize & PCA")
ggsave(plot=pca_memory_combined_by_cell, file.path(plot_dir, "pca_memory_combined_by_cell.svg"), width=3.5, height=3.75)


# Fig 1c
pca_time_combined <- pca_benchmark |>
  filter(threads == 1, tool %in% c("bpcells_stage_float", "delayedarray", "scanpy", "seurat")) |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  inner_join(rna_datasets, by="dataset") |>
  ggplot(aes(nonzero_entries, time_cpu, color=normalize_name[tool], fill=normalize_name[tool])) +
  style_pca +
  labs(y="CPU Time (seconds)", subtitle="Normalize & PCA")
ggsave(plot=pca_time_combined, file.path(plot_dir, "pca_time_combined.svg"), width=3.5, height=3.75)


################################################################################
# Randomized PCA
################################################################################

# Fig S1a
rpca_memory_combined_by_cell <- pca_benchmark |>
  filter(threads == 1, tool %in% c("bpcells_randomized", "scanpy_dask", "bpcells_stage_float")) |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  inner_join(rna_datasets, by="dataset") |>
  ggplot(aes(cells, max_rss, color=normalize_name[tool], group=tool, fill=normalize_name[tool])) +
  style_pca_by_cell +
  labs(y="Memory (GB)", subtitle="Normalize & rPCA (iter=0)")
ggsave(plot=rpca_memory_combined_by_cell, file.path(plot_dir, "rpca_memory_combined_by_cell.svg"), width=3.5, height=3.75)

# Fig S1b
rpca_time_combined <- pca_benchmark |>
  filter(threads == 1, tool %in% c("bpcells_randomized", "scanpy_dask", "bpcells_stage_float")) |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_elapsed), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  inner_join(rna_datasets, by="dataset") |>
  ggplot(aes(nonzero_entries, time_cpu, color=normalize_name[tool], group=tool, fill=normalize_name[tool])) +
  style_pca +
  labs(y="Elapsed time (seconds)", subtitle="Normalize & rPCA (iter=0)")
ggsave(plot=rpca_time_combined, file.path(plot_dir, "rpca_time_combined.svg"), width=3.5, height=3.75)

# Fig S1c
pca_accuracy_plot <- pca_accuracy |>
  filter(dataset=="1m_neurons", axis=="cell", PC_ref==PC_alt, replicate %in% c("rep1", "rep2", "rep3")) |>
  mutate(
    method=if_else(str_detect(alt, "randomized|dask"), "randomized", "exact"),
    alt = case_match(alt, "bpcells_randomized" ~ "BPCells randomized", "bpcells_stage_float" ~ "BPCells exact", "scanpy_dask" ~ "Dask + Scanpy")) |>
  ggplot(aes(PC_ref, abs(pearson), color=alt, group=paste0(replicate, alt))) +
  geom_point() + 
  geom_line() +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(1,NA), breaks=c(1,10,20,30,40,50)) +
  scale_color_manual(values=c("Dask + Scanpy"=set1[3], "BPCells randomized"=RColorBrewer::brewer.pal(5, "Paired")[5], "BPCells exact"=set1[1])) + 
  labs(x="PC", y="Pearson r (abs. value)", subtitle="rPCA accuracy (iter = 0)\n", color="method") +
  guides(color="none") +
  theme_classic(base_size = 14, base_family="Arial")

ggsave(plot=pca_accuracy_plot, file.path(plot_dir, "pca_accuracy_plot.svg"), width=3.5, height=3.75)


################################################################################
# Multi-threaded speedup
################################################################################
# Fig 1d
pca_multithread_speedup <- pca_benchmark |>
  filter(tool == "bpcells_stage_float", step=="pca") |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_elapsed=sum(time_elapsed), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  group_by(tool, dataset, replicate) |>
  mutate(
    relative_time=time_elapsed[threads==1]/time_elapsed,
    threads=factor(threads, levels=c(1,2,4,8,16)),
    .groups="drop"
  ) |>
  inner_join(rna_datasets, by="dataset") |>
  ggplot(aes(nonzero_entries, relative_time, color=threads)) +
  style_pca +
  scale_color_manual(values=RColorBrewer::brewer.pal(9, "Reds")[9:5]) +
  scale_y_continuous(breaks=scales::pretty_breaks()) +
  guides(color="legend") +
  theme(legend.position="inside", legend.position.inside = c(0.01, 1), 
      legend.direction="horizontal", legend.text.position="bottom", legend.title.position="top", 
      legend.key.spacing.x=unit(0, "pt"), legend.justification = c(0, 1.05)) +
  theme(panel.border = element_rect(color="black", linewidth = rel(1), fill=NA)) +
  labs(x="Nonzero entries", y="Relative speedup", color="Threads", subtitle="PCA (multi-threaded BPCells)")

ggsave(plot=pca_multithread_speedup, file.path(plot_dir, "pca_multithread_speedup.svg"), width=3.5, height=3.75)

################################################################################
# Matrix transpose
################################################################################

# Fig S3b
transpose_time <- matrix_transpose |>
  inner_join(rna_datasets, by="dataset") |>
  ggplot(aes(nonzero_entries, time_cpu, color=normalize_name[tool], fill=normalize_name[tool])) +
  style_pca +
  labs(y="CPU Time (seconds)")

# Fig S3b
transpose_memory <- matrix_transpose |>
  inner_join(rna_datasets, by="dataset") |>
  ggplot(aes(nonzero_entries, max_rss/1e9, color=normalize_name[tool], fill=normalize_name[tool])) +
  style_pca +
  scale_y_continuous(transform="log10", guide="axis_logticks", labels = label_si_custom, limits=c(1, NA), expand = expansion(c(0, 0.05))) +
  labs(y="Memory (GB)")

transpose_combined <- transpose_memory + transpose_time + patchwork::plot_annotation(title="Storage order transpose", theme=theme(plot.title=element_text(size=14, family="Arial", hjust=0.5)))
ggsave(plot=transpose_combined, file.path(plot_dir, "transpose_combined.svg"), width=7, height=3.75)


################################################################################
# Marker genes
################################################################################

# Fig S1d
marker_time <- marker_performance |>
  inner_join(rna_datasets, by="dataset") |>
  mutate(tool = case_match(tool, "bpcells"~"BPCells","presto"~"Presto","scanpy"~"Scanpy Uncorrected","scanpy_tiecorrect"~"Scanpy")) |>
  ggplot(aes(nonzero_entries, time_cpu, color=tool, fill=tool)) +
  style_pca +
  labs(y="CPU Time (seconds)")

# Fig S1d
marker_memory_by_cell <- marker_performance |>
  inner_join(rna_datasets, by="dataset") |>
  mutate(tool = case_match(tool, "bpcells"~"BPCells","presto"~"Presto","scanpy"~"Scanpy Uncorrected","scanpy_tiecorrect"~"Scanpy")) |>
  ggplot(aes(cells, max_rss/1e9, color=tool, fill=tool)) +
  style_pca_by_cell +
  labs(y="Memory (GB)")

marker_combined_by_cell <- marker_memory_by_cell + marker_time + patchwork::plot_annotation(title="Marker genes (Wilcoxon test)", theme=theme(plot.title=element_text(size=14, family="Arial", hjust=0.5)))
ggsave(plot=marker_combined_by_cell, file.path(plot_dir, "marker_combined_by_cell.svg"), width=7, height=3.75)


################################################################################
# Comparison labels for figures
################################################################################

# Fig 1b Scanpy 68x more memory, 1c DelayedArray 41x more time on 1M cell dataset
pca_benchmark |>
  filter(threads == 1, tool %in% c("bpcells_stage_float", "delayedarray", "scanpy", "seurat")) |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  summarize(min_time=min(time_cpu), max_time=max(time_cpu), time_cpu=mean(time_cpu), 
            min_mem=min(max_rss), max_mem=max(max_rss), max_rss=mean(max_rss), .groups="drop") |>
  inner_join(rna_datasets, by="dataset") |>
  select(tool, dataset, max_rss, time_cpu) |>
  group_by(dataset) |>
  mutate(relative_mem = max_rss/max_rss[tool=="bpcells_stage_float"], relative_time_cpu = time_cpu/time_cpu[tool=="bpcells_stage_float"]) |>
  filter(dataset == "1m_neurons")

# Fig S1a Scanpy 4x more time than randomized, 24x more memory
pca_benchmark |>
  filter(threads == 1, tool %in% c("bpcells_randomized", "bpcells_stage_float", "scanpy_dask")) |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_elapsed), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  summarize(time_cpu=mean(time_cpu), max_rss=mean(max_rss), .groups="drop") |>
  filter(dataset == "1m_neurons") |>
  group_by(dataset) |>
  mutate(improved_time=time_cpu[tool=="scanpy_dask"] / time_cpu, improved_memory=max_rss[tool=="scanpy_dask"]/max_rss)


################################################################################
# Specific results for text
################################################################################

# Approximate speedup of marker gene calculations is ~2x relative to Presto
marker_performance |>
  group_by(tool, dataset) |>
  summarize(time_cpu=mean(time_cpu), max_rss=mean(max_rss), .groups="drop") |>
  group_by(dataset) |>
  mutate(relative_time=time_cpu/time_cpu[tool=="bpcells"], relative_memory=max_rss/max_rss[tool=="bpcells"]) |>
  arrange(dataset) |>
  group_by(tool) |>
  summarize(geo_mean_relative_time = exp(mean(log(relative_time))))

# PCA requiring ~250-300 passes over the dataset
# Each "op" is a multiply by X and t(X), so two passes. But the final 
# matrix calculation of all 50 PCs takes just one pass but is counted as 50 ops
pca_benchmark |> filter(!is.na(n_ops), str_detect(tool, "bpcells")) |>
  distinct(dataset, n_ops) |> mutate(n_passes = (n_ops-50)*2 + 1)

# BPCells can analyze a 22m cell dataset in 24 minutes
pca_benchmark |>
  filter(threads == 1, tool == "bpcells_stage_float", dataset=="22m_pansci") |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)/1e9, time_elapsed=sum(time_elapsed), .groups="drop_last") |>
  summarize(min_time=min(time_cpu), max_time=max(time_cpu), time_cpu=mean(time_cpu), time_elapsed=mean(time_elapsed),
            min_mem=min(max_rss), max_mem=max(max_rss), max_rss=mean(max_rss), .groups="drop")

# Relative BPCells vs Scanpy 17% geo-mean relative time_cpu
pca_benchmark |>
  filter(threads == 1, tool %in% c("bpcells_stage_float", "scanpy")) |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  summarize(min_time=min(time_cpu), max_time=max(time_cpu), time_cpu=mean(time_cpu), 
            min_mem=min(max_rss), max_mem=max(max_rss), max_rss=mean(max_rss), .groups="drop") |>
  inner_join(rna_datasets, by="dataset") |>
  filter(dataset %in% dataset[tool=="scanpy"]) |>
  select(dataset, tool, time_cpu) |>
  pivot_wider(names_from="tool", values_from="time_cpu") |>
  mutate(relative = bpcells_stage_float/scanpy) |>
  summarize(avg_relative = exp(mean(log(bpcells_stage_float/scanpy))))

# Scanpy memory usage compared to the minimal required on 1m neurons (~5x greater)
pca_benchmark |>
  filter(tool=="scanpy", dataset=="1m_neurons") |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  summarize(min_time=min(time_cpu), max_time=max(time_cpu), time_cpu=mean(time_cpu), 
            min_mem=min(max_rss), max_mem=max(max_rss), max_rss=mean(max_rss), .groups="drop") |>
  inner_join(rna_datasets, by="dataset") |>
  summarize(memory_used=max_rss, memory_required=nonzero_entries*8/1e9, ratio=memory_used/memory_required)
  
# BPCells 22m cells uses about half the memory of scanpy / seurat on 480k cells
pca_benchmark |>
  filter(threads == 1, tool %in% c("bpcells_stage_float", "seurat", "scanpy"), dataset %in% c("22m_pansci", "480k_tabula_sapiens")) |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  summarize(min_time=min(time_cpu), max_time=max(time_cpu), time_cpu=mean(time_cpu), 
            min_mem=min(max_rss), max_mem=max(max_rss), max_rss=mean(max_rss), .groups="drop") |>
  inner_join(rna_datasets, by="dataset")
  
# BPCells performance uplift with 16 threads up to 12x faster
pca_benchmark |>
  filter(tool == "bpcells_stage_float", step=="pca") |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_elapsed=sum(time_elapsed), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  group_by(tool, dataset, replicate) |>
  mutate(
    relative_time=time_elapsed[threads==1]/time_elapsed,
    threads=factor(threads, levels=c(1,2,4,8,16))
  ) |> arrange(desc(relative_time))

# BPCells is fastest tool on all datasets when using 2 threads
pca_benchmark |>
  filter(threads <= 2, tool %in% c("bpcells_stage_float", "delayedarray", "scanpy", "seurat")) |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)/1e9, time_elapsed=sum(time_elapsed), .groups="drop_last") |>
  summarize(min_time=min(time_cpu), max_time=max(time_cpu), time_cpu=mean(time_cpu), time_elapsed=mean(time_elapsed),
            min_mem=min(max_rss), max_mem=max(max_rss), max_rss=mean(max_rss), .groups="drop") |>
  inner_join(rna_datasets, by="dataset") |>
  select(tool, threads, dataset, max_rss, time_cpu, time_elapsed) |>
  group_by(dataset) |>
  mutate(is_fastest = time_elapsed == min(time_elapsed)) |>
  group_by(tool) |>
  summarize(total_fastest = sum(is_fastest))

# Randomized scanpy 24x more memory, 3.9x slower than BPCells randomized,
# Randomizes scanpy 26x more memory, 1.8x slower than BPCells exact
pca_benchmark |>
  filter(threads == 1, tool %in% c("bpcells_randomized", "bpcells_stage_float", "scanpy_dask")) |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_elapsed), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  summarize(time_cpu=mean(time_cpu), max_rss=mean(max_rss), .groups="drop") |>
  filter(dataset == "1m_neurons") |>
  group_by(dataset) |>
  mutate(improved_time=time_cpu[tool=="scanpy_dask"] / time_cpu, improved_memory=max_rss[tool=="scanpy_dask"]/max_rss)

# Marker genes: 154x and 83x more memory, 2.3x and 40x more time 
marker_performance |>
  group_by(tool, dataset) |>
  summarize(time_cpu=mean(time_cpu), max_rss=mean(max_rss), .groups="drop") |>
  group_by(dataset) |>
  mutate(relative_time=time_cpu/time_cpu[tool=="bpcells"], relative_memory=max_rss/max_rss[tool=="bpcells"]) |>
  filter(dataset=="500k_drugscreen")

# Single-threaded PCA can use over 3 GB/s of uncompressed matrix data 
pcs <- 50
pca_benchmark |>
  filter(tool == "bpcells_stage_float", step == "pca", threads==1) |>
  inner_join(pca_benchmark_file_sizes, by=c("dataset", "tool", "threads", "replicate")) |>
  inner_join(cell_metadata_size, by="dataset") |>
  mutate(
    data_passes = (2*(n_ops - pcs) + 1),
    estimated_read_gbps = (file_bytes-cell_metadata_bytes) * data_passes / time_elapsed / 1e9,
    in_memory_gbps = subset_nonzeros * 8 * data_passes / time_elapsed / 1e9
  ) |>
  group_by(tool, dataset) |>
  summarize(estimated_read_gbps = mean(estimated_read_gbps), in_memory_gbps = mean(in_memory_gbps), file_bytes=unique(file_bytes-cell_metadata_bytes))

# Matrix transpose for BPCells uses at worst 1.4 GB of RAM
matrix_transpose |>
  filter(tool == "bpcells") |>
  summarize(max_rss = max(max_rss)/1e9) |>
  pull(max_rss)

max(rna_datasets$nonzero_entries) * 8 / 1e9

# Matrix transpose is 39-52% slower than scipy / R Matrix library
matrix_transpose |>
  group_by(tool, dataset) |>
  summarize(time_cpu = mean(time_cpu), max_rss=mean(max_rss)) |>
  group_by(dataset) |>
  mutate(relative_time_scipy = time_cpu/ifelse("scipy" %in% tool, time_cpu[tool=="scipy"], NA),
         relative_time_R = time_cpu/ifelse("dgCMatrix" %in% tool, time_cpu[tool=="dgCMatrix"], NA)) |>
  inner_join(rna_datasets, by="dataset") |>
  ungroup() |>
  filter(tool=="bpcells") |>
  summarize(relative_time_scipy = exp(mean(log(relative_time_scipy), na.rm=TRUE)),
            relative_time_R = exp(mean(log(relative_time_R), na.rm=TRUE)))

# Uses 147x less memory on 2M cell dataset
matrix_transpose |>
  group_by(tool, dataset) |>
  summarize(time_cpu = mean(time_cpu), max_rss=mean(max_rss)) |>
  group_by(dataset) |>
  mutate(relative_mem_bpcells = max_rss/max_rss[tool=="bpcells"])
  
# Estimated performance for 44M cell dataset with Scanpy
# - On datasets >=1M cells, Scanpy averages 5.5x excess memory usage
scanpy_excess_memory <- pca_benchmark |>
  filter(tool=="scanpy") |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  summarize(min_time=min(time_cpu), max_time=max(time_cpu), time_cpu=mean(time_cpu), 
            min_mem=min(max_rss), max_mem=max(max_rss), max_rss=mean(max_rss), .groups="drop") |>
  inner_join(rna_datasets, by="dataset") |>
  filter(cells > 1e6) |>
  mutate(memory_used=max_rss, memory_required=nonzero_entries*8/1e9, ratio=memory_used/memory_required) |>
  pull(ratio) |>
  mean()
# - Census has 95.6B non-zero entries
cellxgene_nonzero_entries <- read_tsv(file.path(results_dir, "cellxgene/02_matrix_slicing.tsv")) |>
  pull(entries_loaded) |>
  max()
# - Scanpy would estimate memory usage of 4.2TB
cellxgene_nonzero_entries*8*scanpy_excess_memory

# DelayedArray estimated performance for 44M cell dataset:
# - Assume the 41x relative performance gap holds
# - This bakes in an optimistic assumption that DelayedArray will parallelize
#   as well as BPCells despite having to use a slower file format
relative_delayedarray_speed <- pca_benchmark |>
  filter(threads == 1, tool %in% c("bpcells_stage_float", "delayedarray")) |>
  group_by(tool, dataset, threads, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)/1e9, .groups="drop_last") |>
  summarize(min_time=min(time_cpu), max_time=max(time_cpu), time_cpu=mean(time_cpu), 
            min_mem=min(max_rss), max_mem=max(max_rss), max_rss=mean(max_rss), .groups="drop") |>
  inner_join(rna_datasets, by="dataset") |>
  select(tool, dataset, max_rss, time_cpu) |>
  group_by(dataset) |>
  mutate(relative_mem = max_rss/max_rss[tool=="bpcells_stage_float"], relative_time_cpu = time_cpu/time_cpu[tool=="bpcells_stage_float"]) |>
  filter(tool=="delayedarray") |>
  pull(relative_time_cpu) |>
  log() |> mean() |> exp()
# Estimate 1.5-10 days for DelayedArray
bind_rows(
    laptop = read_tsv(file.path(results_dir, "cellxgene/04_pca-laptop/timing.tsv")),
    server = read_tsv(file.path(results_dir, "cellxgene/04_pca/timing.tsv")),
    .id="computer"
  ) |> group_by(computer, input, method) |> 
  summarize(time_elapsed = mean(time_elapsed)) |>
  filter(method=="exact", input=="compress") |>
  mutate(days_delayedarray_est=relative_delayedarray_speed*time_elapsed / (24*60*60))

# Supplemental discussion: Average of 8% of nonzero entries in the variable genes PCA subset
pca_benchmark_file_sizes |>
  inner_join(rna_datasets, by="dataset") |>
  distinct(dataset, subset_nonzeros, nonzero_entries) |>
  mutate(ratio = subset_nonzeros/nonzero_entries) |>
  pull(ratio) |> log() |> mean() |> exp()

# View GB/s disk bandwidth for different thread counts:
# - All of the stage_float results for 1 thread are below 2GB/s, and 2 threads below 4GB/s
# - Over 20GB/s with 16 threads in one instance
pcs <- 50
pca_benchmark |>
  filter(tool %in% c("bpcells_stage_float", "bpcells_stage_int"), step == "pca") |>
  inner_join(pca_benchmark_file_sizes, by=c("dataset", "tool", "threads", "replicate")) |>
  inner_join(cell_metadata_size, by="dataset") |>
  mutate(
    data_passes = (2*(n_ops - pcs) + 1),
    estimated_read_gbps = (file_bytes-cell_metadata_bytes) * data_passes / time_elapsed / 1e9,
    in_memory_gbps = subset_nonzeros * 8 * data_passes / time_elapsed / 1e9
  ) |>
  group_by(tool, threads, dataset) |>
  summarize(estimated_read_gbps = mean(estimated_read_gbps), in_memory_gbps = mean(in_memory_gbps), file_bytes=unique(file_bytes-cell_metadata_bytes)) |>
  pivot_wider(names_from="threads", values_from="estimated_read_gbps", id_cols=c("dataset", "tool"))


# Largest intermediate matrix size was just 7.5GB in size, and in all but 2 cases the intermediate
# matrix size was less than the maximum memory usage for the rest of the process
pca_benchmark |> filter(tool=="bpcells_stage_float", threads==1, step=="pca") |>
  group_by(tool, dataset) |>
  summarize(max_rss=mean(max_rss)) |>
  inner_join(distinct(pca_benchmark_file_sizes, tool, dataset, file_bytes), by=c("tool", "dataset")) |>
  inner_join(cell_metadata_size) |>
  mutate(staged_file_bytes = file_bytes-cell_metadata_bytes, ratio = (file_bytes-cell_metadata_bytes)/max_rss)

# Intermediate matrix got geo-mean 1.6x compression ratio
pca_benchmark_file_sizes |>
  filter(tool=="bpcells_stage_float") |>
  distinct(dataset, file_bytes, subset_nonzeros) |>
  inner_join(cell_metadata_size, by="dataset") |>
  mutate(compression_ratio = 8*subset_nonzeros / (file_bytes-cell_metadata_bytes)) |>
  pull(compression_ratio) |> log() |> mean() |> exp()

# Marker genes with cost of storage order transpose accounted for:
# - On 500k dataset, BPCells 128 with markers only vs 295 for
marker_performance |>
  filter(dataset=="500k_drugscreen", tool %in% c("bpcells", "presto")) |>
  mutate(matrix_format = if_else(tool=="bpcells", "bpcells", "dgCMatrix")) |>
  inner_join(matrix_transpose, by=c(matrix_format="tool", "dataset", "replicate"), suffix = c(".markers", ".transpose")) |>
  group_by(tool, dataset) |>
  summarize(time_cpu.markers=mean(time_cpu.markers), time_cpu.transpose=mean(time_cpu.transpose)) |>
  mutate(time_cpu.markers = if_else(tool=="bpcells", time_cpu.markers, time_cpu.markers-time_cpu.transpose)) |>
  mutate(time_cpu.total = time_cpu.markers + time_cpu.transpose) |>
  ungroup() |>
  mutate(relative_time.total = time_cpu.total/min(time_cpu.total))

# Confirm that all the methods are calculating the same values
pca_accuracy |>
  filter(dataset=="130k_thymus_atlas", axis=="cell", PC_ref==PC_alt) |>
  group_by(ref, alt) |>
  summarize(min_pearson = sprintf("%0.14f", min(abs(pearson))))


# Latex-formatted dataset summary
rna_citation <- tibble::tribble(
  ~dataset, ~citation, ~link,
  "20k_pbmc", "(\\href{https://www.10xgenomics.com/resources/datasets/20-k-human-pbm-cs-3-ht-v-3-1-chromium-x-3-1-high-6-1-0}{10x Genomics})", "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
  "130k_thymus_atlas",  "\\cite{park2020a}", "https://zenodo.org/record/5500511/files/HTA07.A01.v02.entire_data_raw_count.h5ad?download=1",
  "480k_tabula_sapiens", "\\cite{thetabulasapiensconsortium2022}", "https://figshare.com/ndownloader/files/34702114",
  "500k_drugscreen", "(\\href{https://www.10xgenomics.com/resources/datasets/full-chip-mixture-of-drug-treated-h1975-and-a549-cells-3-ht-v3-1-3-1-high}{10x Genomics})", "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/H1975_A549_DrugScreen_3p_HT_nextgem/H1975_A549_DrugScreen_3p_HT_nextgem_count_filtered_feature_bc_matrix.h5",
  "1m_neurons", "(\\href{https://www.10xgenomics.com/resources/datasets/1-3-million-brain-cells-from-e-18-mice-2-standard-1-3-0}{10x Genomics})", "https://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5",
  "2m_perturbseq", "\\cite{replogle2022}", "https://plus.figshare.com/ndownloader/files/35775507",
  "4m_fetal", "\\cite{cao2020}", "https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/",
  "11m_jax", "\\cite{qiu2024}", "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax/download/adata/",
  "22m_pansci", "\\cite{zhang2024a}", "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE247719&format=file&file=GSE247719%5F20240213%5FPanSci%5Fall%5Fcells%5Fadata%2Eh5ad%2Egz"
)

rna_datasets |> 
  arrange(cells) |>
  inner_join(rna_citation, by="dataset") |>
  mutate(
    fraction_nonzero = scales::percent(fraction_nonzero, accuracy=0.1) |> str_replace("%", "\\\\%"),
    nonzero_entries = scales::label_number(accuracy=0.1, scale_cut =scales::cut_short_scale())(nonzero_entries),
    cells = scales::label_comma()(cells),
    genes = scales::label_comma()(genes),
    median_reads = scales::label_comma()(median_reads),
    median_genes = scales::label_comma()(median_genes),
    dataset=str_c(str_replace_all(dataset, "_", " "), " ", citation),
    link=str_c("\\href{", link, "}{\\LinkSymbol}")
    ) |>
  select(dataset, link, cells, genes, median_reads, median_genes, nonzero_entries, fraction_nonzero) |>
  write_delim(file.path(plot_dir, "rna_datasets_latex.txt"), delim="&", eol="\\\\\n", quote="none")


