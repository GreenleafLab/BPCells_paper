library(tidyverse)
library(legendry)
library(patchwork)

################################################################################
# Data loading 
################################################################################

repo_root <- "/home/bparks/Sync/BPCells_paper/final_upload_version/"

results_dir <- file.path(repo_root, "/results/data_tables")
plot_dir <- file.path(repo_root, "/results/plots/cellxgene")
dir.create(plot_dir, showWarnings = FALSE)

write_timing <- bind_rows(
  laptop = read_tsv(file.path(results_dir, "cellxgene/01_subset_unique_cells-laptop.tsv")),
  server = read_tsv(file.path(results_dir, "cellxgene/01_subset_unique_cells.tsv")),
  .id="computer"
) |>
  mutate(
    format = factor(format, levels=c("cellmajor", "genemajor", "norm", "raw", "normalized")),
    format_pretty = recode_factor(format, cellmajor="cell-major", genemajor="feature-major", norm="normalizing\ncoefficients", raw="raw", normalized="normalized"),
    tool= recode_factor(tool, bpcells="BPCells", tiledb="TileDB")
  )


slice_timing <- bind_rows(
  laptop = read_tsv(file.path(results_dir, "cellxgene/02_matrix_slicing-laptop.tsv")),
  server = read_tsv(file.path(results_dir, "cellxgene/02_matrix_slicing.tsv")),
  .id="computer"
) |>
  mutate(
    format=replace_na(format, "raw"),
    format = factor(format, levels=c("cellmajor", "genemajor", "cellmajor_norm", "genemajor_norm", "raw", "norm")),
    format_pretty = recode_factor(format, cellmajor="cell-major", genemajor="feature-major", raw="raw", normalized="normalized"),
  )

mean_variance_timing <- bind_rows(
  laptop = read_tsv(file.path(results_dir, "cellxgene/03_mean_variance-laptop.tsv")),
  server = read_tsv(file.path(results_dir, "cellxgene/03_mean_variance.tsv")),
  .id="computer"
) |>
  mutate(
    normalization = replace_na(normalization, "raw"),
    format = replace_na(format, "tiledb"),
    format_pretty = recode_factor(format, cellmajor="cell-major", genemajor="feature-major", tiledb="TileDB"),
    tool= recode_factor(tool, bpcells="BPCells", tiledb="CELLxGENE")
  )

pca_timing <- bind_rows(
  laptop = read_tsv(file.path(results_dir, "cellxgene/04_pca-laptop/timing.tsv")),
  server = read_tsv(file.path(results_dir, "cellxgene/04_pca/timing.tsv")),
  .id="computer"
)

pca_accuracy <- bind_rows(
  laptop = read_tsv(file.path(results_dir, "cellxgene/04_pca-laptop/pca-cor.tsv.gz")),
  server = read_tsv(file.path(results_dir, "cellxgene/04_pca/pca-cor.tsv.gz")),
  .id="computer"
)

n_genes <- min(pca_accuracy$axis_len)
n_cells <- max(pca_accuracy$axis_len)
full_dataset_slice <- slice_timing |>
  filter(slice_size == n_genes) |>
  select(!c(axis, slice_size, slice_type)) |>
  cross_join(tibble(axis=c("gene", "cell"), slice_size=c(n_genes, n_cells))) |>
  cross_join(tibble(slice_type=c("sequential", "random"))) |>
  select(all_of(names(slice_timing))) |>
  mutate(tool= recode_factor(tool, bpcells="BPCells", tiledb="TileDB"))

slice_timing_extended <- bind_rows(slice_timing |> filter(slice_size != n_genes), full_dataset_slice)

# check for same number of nonzeros loaded in each slice
mismatched_entries <- slice_timing_extended |> 
  arrange(axis, slice_size, slice_type, slice_id, computer) |>
  group_by(axis, slice_size, slice_type, slice_id, computer) |>
  filter(n_distinct(entries_loaded) > 1) |>
  nrow()
stopifnot(mismatched_entries == 0)

################################################################################
# Plot styling helpers
################################################################################

color_tools <- RColorBrewer::brewer.pal(n=3, "Set1")[c(2,1)] |>
  set_names(c("CELLxGENE", "BPCells"))

color_formats <- character(0)
color_formats[c("cellmajor", "norm","genemajor")] <- RColorBrewer::brewer.pal(n=6, "Paired")[5:6] |>
  (function(x) colorRampPalette(x)(3))()
color_formats[c("normalized","raw")] <- RColorBrewer::brewer.pal(n=4, "Paired")[1:2]


color_method <- RColorBrewer::brewer.pal(6, "Paired")[c(6,5)] |>
  set_names("exact", "randomized")

################################################################################
# Read / Write bar charts
################################################################################

cellxgene_read_write_style <- list(
  geom_bar(width=1, position=position_dodge2(width=1, preserve = "single"), stat="summary", fun="mean"),
  geom_text(position=position_dodge2(width=1, preserve="single"), vjust = -0.5, stat="summary", fun="max"),
  geom_point(position=position_dodge2(width=1, preserve="single")),
  scale_y_continuous(transform="log10", guide="axis_logticks", expand = expansion(c(0,0.15)), labels=scales::label_number(scale_cut=scales::cut_short_scale())),
  scale_x_discrete(expand=expansion(add=c(0.50,0.3))),
  scale_fill_manual(values=color_formats),
  guides(fill="none"),
  theme_classic(base_size=14, base_family = "Arial")
)

# Fig 3a: Storage size
cellxgene_size <- write_timing |>
  filter(compression_level == 9 | is.na(compression_level), computer=="server") |>
  group_by(tool, format) |>
  summarize(bytes=median(bytes/1e9), format_pretty=unique(format_pretty), .groups = "drop") |>
  ggplot(aes(tool, bytes, fill=format, label=format_pretty)) +
  cellxgene_read_write_style[-3] +
  scale_y_continuous(breaks=scales::breaks_pretty(), expand = expansion(c(0,0.06))) +
  labs(y="Size (GB)", x = " ", subtitle="Storage size")
ggsave(file.path(plot_dir, "cellxgene_size.svg"), cellxgene_size, width=6, height=4, units="in")

# Fig 3b: laptop write speed
cellxgene_laptop_write <- write_timing |>
  filter(compression_level == 9 | is.na(compression_level), computer=="laptop") |>
  ggplot(aes(tool, time_elapsed/60, fill=format, label=format_pretty)) +
  cellxgene_read_write_style +
  scale_fill_manual(values=color_formats) +
  labs(y="Write time (minutes)", x = " ", subtitle="Laptop write speed")
ggsave(file.path(plot_dir, "cellxgene_laptop_write.svg"), cellxgene_laptop_write, width=6, height=4, units="in")

# Fig S4a: server write speed
cellxgene_server_write <- write_timing |>
  filter(compression_level == 9 | is.na(compression_level), computer=="server") |>
  ggplot(aes(tool, time_elapsed/60, fill=format, label=format_pretty)) +
  cellxgene_read_write_style +
  labs(y="Write time (minutes)", x = " ", subtitle="Server write speed")
ggsave(file.path(plot_dir, "cellxgene_server_write.svg"), cellxgene_server_write, width=6, height=4, units="in")

# Fig 3c: laptop read speed
cellxgene_laptop_read <- full_dataset_slice |>
  filter(computer=="laptop", format!= "cellmajor_norm", slice_type=="sequential", axis=="gene") |>
  mutate(
    format = recode(format, norm="normalized", "genemajor_norm"="norm"),
    format_pretty = recode_factor(format, cellmajor="cell-major", genemajor="feature-major", norm="normalized\nfeature-major")) |>
  ggplot(aes(tool, time_elapsed, fill=format, label=format_pretty)) +
  cellxgene_read_write_style + 
  labs(y="Read time (seconds)", x = " ", subtitle="Laptop read speed")
ggsave(file.path(plot_dir, "cellxgene_laptop_read.svg"), cellxgene_laptop_read, width=6, height=4, units="in")

# Fig S4b: server read speed
cellxgene_server_read <- full_dataset_slice |>
  filter(computer=="server", format!= "cellmajor_norm", slice_type=="sequential", axis=="gene") |>
  mutate(
    format = recode(format, norm="normalized", "genemajor_norm"="norm"),
    format_pretty = recode_factor(format, cellmajor="cell-major", genemajor="feature-major", norm="normalized\nfeature-major")) |>
  ggplot(aes(tool, time_elapsed, fill=format, label=format_pretty)) +
  cellxgene_read_write_style + 
  labs(y="Read time (seconds)", x = " ", subtitle="Server read speed")
ggsave(file.path(plot_dir, "cellxgene_server_read.svg"), cellxgene_server_read, width=6, height=4, units="in")





################################################################################
# Subset reads (random and sequential)
################################################################################

slice_breaks <- function(limits) {
  if (max(limits) > n_cells) {
    # c(1e1, 1e3, 1e5, n_cells)
    c(1e0, 1e2, 1e4, 1e6, n_cells)
  } else {
    c(1, 10, 100, 1e3, 1e4, n_genes)
  }
}

slice_timing_style <- list(
  stat_summary(fun="mean", geom="line"),
  stat_summary(fun="mean", fun.min=function(x) quantile(x, .25), fun.max=function(x) quantile(x, .75), geom="pointrange", size=0.1),
  scale_x_continuous(trans="log10", breaks=slice_breaks,  labels=scales::label_number(scale_cut = scales::cut_short_scale(), accuracy=1)),
                     # guide = guide_axis_logticks(mid=0.75, short.theme = element_line(linewidth=unit(0.5, "pt")))),
  scale_y_continuous(trans="log10", guide="axis_logticks"),
  scale_color_manual(values=color_formats, labels=c("genemajor"="BPCells feature-major", "cellmajor"="BPCells cell-major", raw="TileDB")),
  theme_classic(base_size=14, base_family = "Arial"),
  guides(color="none"),
  labs(y="Elapsed time (s)")
)

# Fig 3d
slice_cell_random <- slice_timing_extended |>
  filter(computer == "laptop", slice_type=="random", axis=="cell", is.na(format) | !str_detect(format, "norm")) |>
  ggplot(aes(slice_size, time_elapsed, color=format)) +
  slice_timing_style +
  labs(x = "Cells selected", subtitle="Random Subsets Laptop") 
  
slice_gene_random <- slice_timing_extended |>
  filter(computer == "laptop", slice_type=="random", axis=="gene", is.na(format) | !str_detect(format, "norm")) |>
  ggplot(aes(slice_size, time_elapsed, color=format)) +
  slice_timing_style +
  labs(x = "Genes selected") 

cellxgene_slice_random_laptop <- slice_cell_random + slice_gene_random +
  patchwork::plot_layout(axis_titles = "collect", guides = "collect") &
  theme(legend.position="top")
ggsave(file.path(plot_dir, "cellxgene_slice_random_laptop.svg"), cellxgene_slice_random_laptop, width=7, height=4.125, units="in")


# Fig S4e
slice_cell_random_server <- slice_timing_extended |>
  filter(computer == "server", slice_type=="random", axis=="cell", is.na(format) | !str_detect(format, "norm")) |>
  ggplot(aes(slice_size, time_elapsed, color=format)) +
  slice_timing_style +
  labs(x = "Cells selected") 

slice_gene_random_server <- slice_timing_extended |>
  filter(computer == "server", slice_type=="random", axis=="gene", is.na(format) | !str_detect(format, "norm")) |>
  ggplot(aes(slice_size, time_elapsed, color=format)) +
  slice_timing_style +
  labs(x = "Genes selected") 


# Fig S4g
slice_cell_sequential <- slice_timing_extended |>
  filter(computer == "laptop", slice_type=="sequential", axis=="cell", is.na(format) | !str_detect(format, "norm")) |>
  ggplot(aes(slice_size, time_elapsed, color=format)) +
  slice_timing_style +
  labs(x = "Cells selected") 

slice_gene_sequential <- slice_timing_extended |>
  filter(computer == "laptop", slice_type=="sequential", axis=="gene", is.na(format) | !str_detect(format, "norm")) |>
  ggplot(aes(slice_size, time_elapsed, color=format)) +
  slice_timing_style +
  labs(x = "Genes selected") 


# Fig S4f
slice_cell_sequential_server <- slice_timing_extended |>
  filter(computer == "server", slice_type=="sequential", axis=="cell", is.na(format) | !str_detect(format, "norm")) |>
  ggplot(aes(slice_size, time_elapsed, color=format)) +
  slice_timing_style +
  labs(x = "Cells selected") 

slice_gene_sequential_server <- slice_timing_extended |>
  filter(computer == "server", slice_type=="sequential", axis=="gene", is.na(format) | !str_detect(format, "norm")) |>
  ggplot(aes(slice_size, time_elapsed, color=format)) +
  slice_timing_style +
  labs(x = "Genes selected") 

# Combined layout of S4e,f,g
layout <- "
  33
  12
  66
  45
  99
  78
"
manual_label <- function(l) patchwork::wrap_elements(grid::textGrob(l, x=unit(0, "npc"), y=unit(.95, "npc"), just=c("left", "top"), gp=grid::gpar(fontsize=14, fontfamily="Arial")))
cellxgene_slice_extras <- (slice_cell_random_server + slice_gene_random_server + manual_label("Random Subsets Server") +
                             slice_cell_sequential_server + slice_gene_sequential_server + manual_label("Sequential Subsets Server") +
                             slice_cell_sequential + slice_gene_sequential + manual_label("Sequential Subsets Laptop") &
                             theme(legend.position="bottom") & guides(color="none")) +
  patchwork::plot_layout(design=layout, guides = "collect", heights = rep(grid::unit.c(unit(16, "pt"), unit(1, "null")), 3))
ggsave(file.path(plot_dir, "cellxgene_slice_extras.svg"), cellxgene_slice_extras, width=7, height=12, units="in")

################################################################################
# Gene stats
################################################################################

# Fig S4c
mean_variance_server <- mean_variance_timing |>
  filter(computer == "server", normalization=="raw", axis=="gene", format %in% c("tiledb", "genemajor")) |>
  ggplot(aes(tool, time_elapsed/60, fill=tool)) +
  geom_bar(stat="summary", fun="mean") +
  scale_y_continuous(expand=expansion(c(0,0.05)), breaks=scales::breaks_width(2)) +
  scale_fill_manual(values = color_tools) +
  labs(x="", y="Elapsed time (minutes)", subtitle="Server Gene Stats") +
  theme_classic(base_size=14, base_family = "Arial") +
  geom_point(position=position_dodge2(width=0.9)) +
  guides(fill="none")
ggsave(file.path(plot_dir, "cellxgene_mean_variance_server.svg"), mean_variance_server, width=3, height=4, units="in")

# Fig 3e
mean_variance_laptop <- mean_variance_timing |>
  filter(computer == "laptop", normalization=="raw", axis=="gene", format %in% c("tiledb", "genemajor")) |>
  ggplot(aes(tool, time_elapsed/60, fill=tool)) +
  geom_bar(stat="summary", fun="mean") +
  scale_y_continuous(expand=expansion(c(0,0.05)), breaks=scales::breaks_width(2)) +
  scale_fill_manual(values = color_tools) +
  labs(x="", y="Elapsed time (minutes)", subtitle="Laptop Gene Stats") +
  theme_classic(base_size=14, base_family = "Arial") +
  geom_point(position=position_dodge2(width=0.9)) +
  guides(fill="none")
ggsave(file.path(plot_dir, "cellxgene_mean_variance_laptop.svg"), mean_variance_laptop, width=3, height=4, units="in")

# Fig S4d
cellxgene_mean_variance_cpu_time <- mean_variance_timing |>
  filter(normalization=="raw", axis=="gene", format %in% c("tiledb", "genemajor")) |>
  mutate(label = factor(paste0(computer, ".", tool), unique(paste0(computer, ".", tool))[c(1,3,2,4)])) |>
  ggplot(aes(label, time_cpu/60, fill=tool)) +
  geom_bar(stat="summary", fun="mean") +
  scale_y_continuous(expand=expansion(c(0,0.05))) +
  scale_fill_manual(values = color_tools) +
  labs(x="", y="CPU time (minutes)", subtitle="Gene Stats (CPU time)") +
  theme_classic(base_size=14, base_family = "Arial") +
  geom_point(position=position_dodge2(width=0.9)) +
  guides(x = legendry::guide_axis_nested(key=key_range_auto(sep="\\.")), fill="none")
ggsave(file.path(plot_dir, "cellxgene_mean_variance_cpu_time.svg"), cellxgene_mean_variance_cpu_time, width=3, height=4.25, units="in")


################################################################################
# PCA timing
################################################################################

# Confirm that our y limits won't cut anything off
stopifnot(10 > max(pca_timing$time_elapsed/3600))

# Fig 3f
pca_timing_laptop <- pca_timing |>
  filter(computer=="laptop") |>
  arrange(method, input) |>
  mutate(label=paste0(input, "ed.", method)) |>
  ggplot(aes(label, time_elapsed/3600, fill=method)) +
  geom_bar(stat="summary", fun="mean") +
  scale_y_continuous(expand=expansion(c(0,0)), limits=c(0, 10), breaks=scales::breaks_width(1)) +
  scale_fill_manual(values=color_method) +
  guides(x = legendry::guide_axis_nested(key=key_range_auto(sep="\\.", reverse=TRUE)), fill="none") +
  labs(x="", y="Elapsed time (hours)", subtitle="Laptop PCA speed") +
  geom_point(position=position_dodge2(width=0.9)) +
  theme_classic(base_size=14, base_family = "Arial") 
ggsave(file.path(plot_dir, "cellxgene_pca_timing_laptop.svg"), pca_timing_laptop, width=4, height=4.25, units="in")

# Fig 3h
pca_timing_server <- pca_timing |>
  filter(computer=="server") |>
  arrange(method, input) |>
  mutate(label=paste0(input, "ed.", method)) |>
  ggplot(aes(label, time_elapsed/3600, fill=method)) +
  geom_bar(stat="summary", fun="mean") +
  scale_y_continuous(expand=expansion(c(0,0)), limits=c(0, 10), breaks=scales::breaks_width(1)) +
  scale_fill_manual(values=color_method) +
  guides(x = legendry::guide_axis_nested(key=key_range_auto(sep="\\.", reverse=TRUE)), fill="none") +
  labs(x="", y="Elapsed time (hours)", subtitle="Server PCA speed") +
  geom_point(position=position_dodge2(width=0.9)) +
  theme_classic(base_size=14, base_family = "Arial")
ggsave(file.path(plot_dir, "cellxgene_pca_timing_server.svg"), pca_timing_server, width=4, height=4.25, units="in")

# Fig 3g
pca_accuracy_plot <- pca_accuracy |>
  filter(ref=="compress_exact_rep1", axis=="cell", PC_ref==PC_alt) |>
  mutate(method=if_else(str_detect(alt, "randomized"), "randomized", "exact")) |>
  arrange(desc(method)) |>
  ggplot(aes(PC_ref, abs(pearson), color=method, group=paste0(alt, computer))) +
  geom_point() + 
  geom_line() +
  scale_y_continuous(limits=c(0,1), expand = expansion(c(0, 0.05))) +
  scale_x_continuous(limits=c(1,NA), breaks=c(1,10,20,30)) +
  scale_color_manual(values=color_method) + 
  guides(color="none") +
  labs(x="PC", y="Pearson r (abs. value)", subtitle="rPCA accuracy (iter = 2)") +
  theme_classic(base_size=14, base_family = "Arial") +
  theme(legend.position="inside", legend.position.inside = c(0.05,0.05), legend.justification = c(0,0)) 
  
ggsave(file.path(plot_dir, "cellxgene_pca_accuracy.svg"), pca_accuracy_plot, width=3, height=4, units="in")


################################################################################
# Comparison labels for figures
################################################################################

### Size calculations for Fig 1a diagram
# 44M cells, 60k genes
range(pca_accuracy$axis_len)
# Raw counts on disk is about 120-170GB on disk
write_timing |> filter(format %in% c("cellmajor", "genemajor", "raw"), is.na(compression_level) | compression_level == 9) |>
  mutate(bytes=bytes/1e9) |> pull(bytes) |> range()
# Size of in-memory matrices is about 750GB
max(slice_timing_extended$entries_loaded) * 8/1e9
# BPCells outputs buffer size is about 350MB
max(pca_accuracy$axis_len) * 8/1e6
# Output 50 PCs size is about 8.8 GB (assuming 4-byte floats)
max(pca_accuracy$axis_len) * 50 * 4 / 1e9


# Fig 3b laptop write: cellmajor -> genemajor = 3x, genemajor->raw = 41x  
write_timing |>
  filter(computer=="laptop") |>
  group_by(computer, tool, format, compression_level) |>
  summarize(write_time_cpu=mean(time_cpu), write_time_elapsed=mean(time_elapsed), bytes=mean(bytes), .groups="drop") |>
  mutate(relative_write_cellmajor = write_time_elapsed / write_time_elapsed[format=="cellmajor"],
         relative_write_genemajor = write_time_elapsed / write_time_elapsed[format=="genemajor"])

# Fig 3c tiledb raw is 21x slower than BPCells genemajor
full_dataset_slice |>
  filter(format!= "cellmajor_norm", slice_type=="sequential", axis=="gene", computer=="laptop") |>
  group_by(computer, tool, format) |>
  summarize(read_time_cpu=mean(time_cpu), read_time_elapsed=mean(time_elapsed), .groups="drop") |>
  mutate(relative_read = read_time_elapsed / read_time_elapsed[format=="genemajor"])

# Fig S4a Server write: cellmajor -> genemajor = 3x, genemajor->raw = 63x  
write_timing |>
  filter(computer=="server", is.na(compression_level) | compression_level == 9) |>
  group_by(computer, tool, format, compression_level) |>
  summarize(write_time_cpu=mean(time_cpu), write_time_elapsed=mean(time_elapsed), bytes=mean(bytes), .groups="drop") |>
  mutate(relative_write_cellmajor = write_time_elapsed / write_time_elapsed[format=="cellmajor"],
         relative_write_genemajor = write_time_elapsed / write_time_elapsed[format=="genemajor"])

# Fig S4b tiledb raw is 41x slower than BPCells genemajor
full_dataset_slice |>
  filter(format!= "cellmajor_norm", slice_type=="sequential", axis=="gene", computer=="server") |>
  group_by(computer, tool, format) |>
  summarize(read_time_cpu=mean(time_cpu), read_time_elapsed=mean(time_elapsed), .groups="drop") |>
  mutate(relative_read = read_time_elapsed / read_time_elapsed[format=="genemajor"])

# Fig 3d
#  - Laptop BPCells cell-major vs tiledb 10k cells: 65x faster
#  - Laptop BPCells feature-major vs tiledb 10 genes: 260x faster
# Fig S3e
#  - Server BPCells cell-major vs tiledb 10k cells: 16x faster
#  - Server BPCells feature-major vs tiledb 10 genes:  78x faster
slice_timing_extended |>
  mutate(tool=str_to_lower(tool)) |>
  group_by(computer, tool, format, axis, slice_size, slice_type) |>
  summarize(time_elapsed=mean(time_elapsed), .groups="drop") |>
  filter(format %in% c("cellmajor", "genemajor", "raw")) |>
  pivot_wider(names_from=c("tool", "format"), values_from=time_elapsed) |>
  mutate(relative_genemajor = tiledb_raw / bpcells_genemajor, relative_cellmajor = tiledb_raw / bpcells_cellmajor) |>
  filter((axis == "cell" & slice_size == 10000) | (axis=="gene" & slice_size == 10), slice_type=="random") |>
  arrange(desc(slice_size))

# Fig 3f+h: Uncompressed vs compressed inputs for exact PCA: 1.5x on laptop, 9x on server
pca_timing |> 
  filter(method=="exact") |>
  group_by(computer, input) |>
  summarize(time_elapsed = mean(time_elapsed)) |>
  group_by(computer) |>
  mutate(relative_time = time_elapsed/time_elapsed[input=="compress"])

################################################################################
# Results for text
################################################################################
# Atlas scale gene expression variance is about 10x faster on a laptop (really 11-12x depending on which 
# BPCells format we use as baseline)
left_join(
  mean_variance_timing |>
    filter(normalization=="raw", axis=="gene"),
  write_timing |>
    filter(compression_level == 9 | is.na(compression_level), format %in% c("cellmajor", "genemajor", "raw")) |>
    group_by(tool, format) |>
    summarize(bytes=median(bytes), .groups="drop") |>
    mutate(format=if_else(tool=="TileDB", "tiledb", format)) |>
    select(format, bytes),
  by="format",
  relationship="many-to-one"
) |>
  group_by(computer, tool, format) |>
  summarize(time_elapsed=mean(time_elapsed), time_cpu=mean(time_cpu), bytes=unique(bytes), .groups="drop") |>
  group_by(computer) |>
  mutate(read_GBps=bytes/time_elapsed/1e9,
         relative_time_cpu=time_cpu/min(time_cpu),
         relative_time_elapsed=time_elapsed/min(time_elapsed))


# CELLxGENE file sizes:
#   - TileDB raw: 147GB, normalized: 166GB (using default compression level 9)
#   - BPCells cellmajor 168GB, genemajor 117GB
write_timing |> 
  group_by(tool, format, compression_level) |>
  summarize(size_gb=min(bytes)/1e9)

# Write speed relative to TileDB: 41x-125x faster on laptop, 63-186x faster on server (comparing raw data writes, not normalized)
write_timing |>
  filter(is.na(compression_level) | compression_level == 9) |>
  group_by(computer, tool, format, compression_level) |>
  summarize(write_time_cpu=mean(time_cpu), write_time_elapsed=mean(time_elapsed), bytes=mean(bytes), .groups="drop") |>
  group_by(computer) |>
  mutate(relative_write_raw =  write_time_elapsed[format=="raw"]/ write_time_elapsed,
         relative_write_norm =  write_time_elapsed[format=="normalized"]/ write_time_elapsed)

# Read speed relative to TileDB: 11-21x faster on laptop, 15-41x faster on server
full_dataset_slice |>
  filter(format!= "cellmajor_norm", slice_type=="sequential", axis=="gene") |>
  group_by(computer, tool, format) |>
  summarize(read_time_cpu=mean(time_cpu), read_time_elapsed=mean(time_elapsed), .groups="drop") |>
  group_by(computer) |>
  mutate(relative_read_raw = read_time_elapsed[tool=="TileDB" & format=="raw"] / read_time_elapsed,
         relative_read_norm = read_time_elapsed[tool=="TileDB" & format=="norm"] / read_time_elapsed)

# Relative feature-major and cell-major reads
full_dataset_slice |>
  filter(format %in% c("genemajor", "cellmajor"), slice_type=="sequential", axis=="gene") |>
  group_by(computer, tool, format) |>
  summarize(read_time_cpu=mean(time_cpu), read_time_elapsed=mean(time_elapsed), .groups="drop") |>
  group_by(computer) |>
  mutate(relative_read = max(read_time_elapsed) / read_time_elapsed)

# Storage space of normalization coefficients: 0.35GB
write_timing |> filter(format=="norm") |> mutate(size_gb = bytes/1e9) |> distinct(size_gb)

# On-the-fly normalization: 1.8-2.7x more time on-the-fly, still 11-15x faster than TileDB
full_dataset_slice |>
  filter(format %in% c("genemajor", "genemajor_norm", "norm"), slice_type=="sequential", axis=="gene") |>
  group_by(computer, tool, format) |>
  summarize(read_time_cpu=mean(time_cpu), read_time_elapsed=mean(time_elapsed), .groups="drop") |>
  group_by(computer) |>
  mutate(relative_read_bpcells_norm = read_time_elapsed[format=="genemajor_norm"] / read_time_elapsed,
         relative_read_tiledb_norm = read_time_elapsed[format=="norm"]/read_time_elapsed) 

# 10k cell subset: BPCells cellmajor .9 and 1.1 seconds laptop/server, TileDB 60 and 18 seconds laptop/server
slice_timing_extended |>
  filter(slice_size==10000, slice_type=="random", axis=="cell") |>
  mutate(tool=str_to_lower(tool)) |>
  group_by(computer, tool, format, axis, slice_size, slice_type) |>
  summarize(time_elapsed=mean(time_elapsed), .groups="drop")

# 10 genes: BPCells .5-.7 seconds laptop/server, TileDB 143/54
slice_timing_extended |>
  filter(slice_size==10, slice_type=="random", axis=="gene") |>
  mutate(tool=str_to_lower(tool)) |>
  group_by(computer, tool, format, axis, slice_size, slice_type) |>
  summarize(time_elapsed=mean(time_elapsed), .groups="drop")


# Mean variance: 
#  - BPCells 96s laptop, 32s server
#  - 4.6/5.8x less CPU time on laptop/server compared to CELLxGENE
#  - 11/38x less wall-clock time
mean_variance_timing |>
  filter(normalization=="raw", axis=="gene", format %in% c("tiledb", "genemajor")) |>
  group_by(computer, tool, format) |>
  summarize(time_elapsed=mean(time_elapsed), time_cpu=mean(time_cpu), .groups="drop") |>
  group_by(computer) |>
  mutate(relative_time_cpu=time_cpu/min(time_cpu),
         relative_time_elapsed=time_elapsed/min(time_elapsed))

# PCA timing: 
#  - bitpacked genemajor exact: 6.2 hours on lapotop, 51 minutes on server
pca_timing |>
  group_by(computer, input, method) |>
  summarize(hours_cpu=mean(time_cpu)/(60*60), hours_elapsed=mean(time_elapsed)/(60*60), minutes_elapsed=mean(time_elapsed)/60)

# PCA passes:
# - exact: 167 passes: n_ops = 115 -> passes = (115-32)*2 + 1 = 167 (two passes per multiply, except for a final dense multiply which counts 32 ops in 1 pass)
# - This was done with a manual measurement of the `nops` field on the outputs in cellxgene-census/2024-07-01/pca/compress_exact_*.rds

# Randomized PCA accuracy: above .99 correlation up to component 16, then starts to drop
pca_accuracy |>
  filter(ref=="compress_exact_rep1", axis=="cell", PC_ref==PC_alt, str_detect(alt, "randomized")) |>
  group_by(PC_ref) |>
  summarize(pearson_min=min(abs(pearson)), pearson_max=max(abs(pearson))) |>
  arrange(PC_ref) |>
  filter(PC_ref < 21, PC_ref > 14)

# Data read rates for PCA on:
#  - server 6.35GB/s compressed read, 4.56GB/s uncompressed read
#  - laptop 0.9GB/s compressed, 3.9GB/s uncompressed 
#  Also uncompressed matrix size 765GB
n_ops <- 115 
data_passes <- 2*(n_ops-32) + 1
pca_size_stats <- tibble::tibble(
  input = c("compress", "uncompress"),
  bytes = c(write_timing |> filter(format=="genemajor") |> pull(bytes) |> unique(), unique(full_dataset_slice$entries_loaded)*8)
)
pca_timing |>
  filter(method=="exact") |>
  group_by(computer, input, method) |>
  summarize(time_elapsed=mean(time_elapsed), .groups="drop") |>
  inner_join(pca_size_stats, by="input") |>
  mutate(avg_read_GBps = bytes*data_passes/time_elapsed/1e9) |>
  arrange(computer, method)

# Supplement: 
#  - subsets of 1-100 cells ~0.5 seconds
#  - subsets of 10k cells ~1 seconds
#  - subsets of 100k cells ~3 seconds
slice_timing_extended |>
  filter(slice_size<=100000, slice_type=="random", axis=="cell", format=="cellmajor") |>
  mutate(tool=str_to_lower(tool)) |>
  group_by(computer, tool, format, axis, slice_size, slice_type) |>
  summarize(time_elapsed=mean(time_elapsed), .groups="drop")

# Methods: Using level 1 compressionsaves 4-12% of time at a cost of 10-18% larger files
write_timing |> 
  filter(computer=="server", tool=="TileDB") |>
  group_by(tool, format, compression_level) |>
  summarize(time_cpu = mean(time_cpu), time_elapsed=mean(time_elapsed), bytes=mean(bytes)) |>
  group_by(tool, format) |>
  mutate(percent_faster=100 - 100*time_elapsed/max(time_elapsed), relative_bytes=bytes/min(bytes))



  
