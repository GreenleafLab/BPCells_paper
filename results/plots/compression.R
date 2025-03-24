library(tidyverse)
library(patchwork)

################################################################################
# Data loading 
################################################################################

repo_root <- "/home/bparks/Sync/BPCells_paper/final_upload_version/"

results_dir <- file.path(repo_root, "/results/data_tables")
plot_dir <- file.path(repo_root, "/results/plots/compression")
dir.create(plot_dir, showWarnings = FALSE)

in_memory <- read_tsv(file.path(results_dir, "compression/in-memory-compression.tsv.gz"))

rna_1m_cells <- read_tsv(file.path(results_dir, "compression/rna-1M-cell.tsv"))

fragment_read <- read_tsv(file.path(results_dir, "compression/fragments-read-write/read.tsv.gz"))
fragment_import <- read_tsv(file.path(results_dir, "compression/fragments-read-write/import.tsv.gz"))



# Lzbench doesn't get rid of 8-byte header, but otherwise things should match
in_memory %>%
  distinct(sample, field, input_bytes) %>%
  arrange(field, sample)

in_memory_summary <- in_memory %>%
  mutate(
    filter=replace_na(filter, "na"),
    level=replace_na(level, 0)
  ) %>%
  group_by(sample, codec, runner, filter, level, replicate) %>%
  summarize(
    write_min=sum(write_min),
    read_min=sum(read_min),
    bytes=sum(bytes),
    input_gb=sum(input_bytes)/1e9,
    fields=n(),
    .groups="drop"
  ) %>%
  group_by(sample, codec, filter, level) %>%
  summarize(
    write_speed = median(input_gb/write_min),
    write_lo = min(input_gb/write_min),
    write_hi = max(input_gb/write_min),
    read_speed = median(input_gb/read_min),
    read_lo = min(input_gb/read_min),
    read_hi = max(input_gb/read_min),
    compression_ratio=unique(input_gb*1e9/bytes),
    .groups="drop"
  ) %>%
  mutate(
    codec = recode(codec, lz4fast="lz4", lz4hc="lz4", zstd_fast="zstd", zlib="gzip", bpcells="BPCells"),
    codec=as.factor(paste0(codec, if_else(filter=="byte" & codec != "blosclz", "-blosc", ""))),
    codec = fct_relevel(codec, "BPCells", "blosclz"),
    label = if_else(codec=="BPCells", "BPCells", "")
  )


################################################################################
# Plot styling helpers
################################################################################

set1 <- RColorBrewer::brewer.pal(9, name="Set1")
compression_colors <- RColorBrewer::brewer.pal(12, name="Paired")
names(compression_colors) <- c("zstd-blosc", "zstd", "gzip-blosc", "gzip", "NA", "BPCells",
                               "blosclz", "NA", "lz4-blosc", "lz4", "NA", "memcpy")

style_in_memory_compression <- list(
  geom_point(key_glyph=draw_key_rect),
  geom_linerange(key_glyph=draw_key_rect),
  geom_line(key_glyph=draw_key_rect),
  geom_text(hjust=1.1, size=10/.pt, vjust=0.4),
  scale_color_manual(values=compression_colors),
  scale_y_continuous(labels=scales::label_number(big.mark = ","), limits=c(0, NA), breaks=scales::breaks_width(1), expand = expansion(c(0,0.05))),
  scale_x_continuous(breaks=scales::breaks_width(1), limits=c(1, NA), expand=expansion(c(0,0.05))),
  guides(color="none"),
  labs(x="Compression Ratio"),
  theme_classic(base_size = 14, base_family="Arial")
)

################################################################################
# In-memory compression benchmarks
################################################################################

atac_read <-in_memory_summary %>%
  filter(filter %in% c("na", "byte"), codec != "memcpy", sample=="10k_pbmc") %>%
  ggplot(aes(compression_ratio, read_speed, ymin=read_lo, ymax=read_hi, 
             color=codec, label=label)) +
  style_in_memory_compression +
  labs(y="Read Speed (GB/s)")

atac_write <- in_memory_summary %>%
  filter(filter %in% c("na", "byte"), codec != "memcpy", sample=="10k_pbmc") %>%
  ggplot(aes(compression_ratio, write_speed, ymin=write_lo, ymax=write_hi, 
             color=codec, label=label)) +
  style_in_memory_compression +
  labs(y="Write Speed (GB/s)")

# Fig 2j
atac_in_memory_read_write <- atac_read + atac_write + patchwork::plot_annotation(title="ATAC in-memory compression benchmark", theme=theme(plot.title=element_text(size=14, family="Arial", hjust=0.5)))
ggsave(file.path(plot_dir, "atac_in_memory_read_write.svg"), atac_in_memory_read_write, width=5.6, height=2.8)


rna_transpose_read <-in_memory_summary %>%
  filter(filter %in% c("na", "byte"), codec != "memcpy", sample=="130k_thymus_transpose") %>%
  ggplot(aes(compression_ratio, read_speed, ymin=read_lo, ymax=read_hi, 
             color=codec, label=label)) +
  style_in_memory_compression +
  scale_x_continuous(breaks=scales::breaks_width(2), limits=c(1, NA), expand=expansion(c(0,0.05))) +
  labs(y="Read Speed (GB/s)")

rna_transpose_write <- in_memory_summary %>%
  filter(filter %in% c("na", "byte"), codec != "memcpy", sample=="130k_thymus_transpose") %>%
  ggplot(aes(compression_ratio, write_speed, ymin=write_lo, ymax=write_hi, 
             color=codec, label=label)) +
  style_in_memory_compression +
  scale_x_continuous(breaks=scales::breaks_width(2), limits=c(1, NA), expand=expansion(c(0,0.05))) +
  labs(y="Write Speed (GB/s)")

# Fig 2e
rna_genemajor_in_memory_read_write <- rna_transpose_read + rna_transpose_write + patchwork::plot_annotation(title="RNA in-memory compression benchmark (feature-major)", theme=theme(plot.title=element_text(size=14, family="Arial", hjust=0.5)))
ggsave(file.path(plot_dir, "rna_genemajor_in_memory_read_write.svg"), rna_genemajor_in_memory_read_write, width=5.6, height=2.8)


rna_read <-in_memory_summary %>%
  filter(filter %in% c("na", "byte"), codec != "memcpy", sample=="130k_thymus") %>%
  ggplot(aes(compression_ratio, read_speed, ymin=read_lo, ymax=read_hi, 
             color=codec, label=label)) +
  style_in_memory_compression +
  scale_x_continuous(breaks=scales::breaks_width(2), limits=c(1, NA), expand=expansion(c(0,0.05))) +
  labs(y="Read Speed (GB/s)")

rna_write <- in_memory_summary %>%
  filter(filter %in% c("na", "byte"), codec != "memcpy", sample=="130k_thymus") %>%
  ggplot(aes(compression_ratio, write_speed, ymin=write_lo, ymax=write_hi, 
             color=codec, label=label)) +
  style_in_memory_compression +
  scale_x_continuous(breaks=scales::breaks_width(2), limits=c(1, NA), expand=expansion(c(0,0.05))) +
  labs(y="Write Speed (GB/s)")

# Fig S3a
rna_cellmajor_in_memory_read_write <- rna_read + rna_write + patchwork::plot_annotation(title="RNA in-memory compression benchmark (cell-major)", theme=theme(plot.title=element_text(size=14, family="Arial", hjust=0.5)))
ggsave(file.path(plot_dir, "rna_cellmajor_in_memory_read_write.svg"), rna_cellmajor_in_memory_read_write, width=7, height=3.5)

in_memory_legend <- in_memory_summary %>%
  filter(filter %in% c("na", "byte"), codec != "memcpy", sample=="130k_thymus") %>%
  mutate(codec=fct_relabel(codec, function(x) str_replace(x, "-blosc", " byte-shuffle"))) %>%
  ggplot(aes(compression_ratio, read_speed, ymin=read_lo, ymax=read_hi, 
             color=codec, label=label)) +
  style_in_memory_compression + 
  scale_color_manual(values=compression_colors |> set_names(str_replace(names(compression_colors), "-blosc", " byte-shuffle"))) +
  guides(color=guide_legend(ncol=2)) + 
  theme(legend.key.spacing.y=unit(2, "pt"), legend.byrow = TRUE) +
  labs(color="Compression Algorithms\n(in-memory benchmarks)")
ggsave(file.path(plot_dir, "in_memory_legend.svg"), in_memory_legend, width=5, height=2.8)


in_memory_legend_narrow <- in_memory_summary %>%
  filter(filter %in% c("na", "byte"), codec != "memcpy", sample=="130k_thymus") %>%
  mutate(codec=fct_relabel(codec, function(x) str_replace(x, "-blosc", " byte\nshuffle"))) %>%
  ggplot(aes(compression_ratio, read_speed, ymin=read_lo, ymax=read_hi, 
             color=codec, label=label)) +
  style_in_memory_compression + 
  scale_color_manual(values=compression_colors |> set_names(str_replace(names(compression_colors), "-blosc", " byte\nshuffle"))) +
  guides(color=guide_legend(ncol=1)) + 
  theme(legend.key.spacing.y=unit(1, "pt"), legend.byrow = TRUE) +
  labs(color=NULL)
ggsave(file.path(plot_dir, "in_memory_legend_narrow.svg"), in_memory_legend_narrow, width=5, height=2.8)


################################################################################
# RNA 1M cell example
################################################################################
format_colors <- c("10x" = set1[3], "BPCells"=set1[1], "h5ad"=set1[2], "zarr"=set1[4])
normalize_format_name <- factor(c(
    bpcells="BPCells", "10x"="10x", h5ad="h5ad", h5ad_zarr="zarr"
  ), levels=c("10x", "h5ad", "zarr", "BPCells"))

style_storage_barplot <- list(
  geom_bar = geom_bar(position=position_dodge2(preserve = "single"), stat="summary", fun="mean"),
  geom_point = geom_point(position=position_dodge2(width=0.9), key_glyph="blank"),
  scale_y = scale_y_continuous(expand=expansion(c(0, 0.05))),
  scale_x = scale_x_discrete(guide=guide_axis(angle=30)),
  scale_fill = scale_fill_manual(values=format_colors),
  theme = theme_classic(base_size=14, base_family = "Arial"),
  guides = guides(fill="none"),
  labs = labs(x=NULL)
)

# Fig 2b
rna_size_example <- rna_1m_cells |>
  filter(tool != "cp", format %in% c("10x", "bpcells", "h5ad", "h5ad_zarr")) |>
  mutate(format = normalize_format_name[format]) |>
  group_by(replicate, format) |>
  summarise(read=min(read_seconds_cpu), write=min(write_seconds_cpu), size=unique(size), .groups="drop") |>
  ggplot(aes(format, size/1e9, fill=format)) +
  style_storage_barplot[names(style_storage_barplot) != "geom_point"] +
  labs(y="Size (GB)", subtitle="RNA size (1M cells)")
ggsave(file.path(plot_dir, "rna_size_example.svg"), rna_size_example, width=2.8, height=2.8, units="in")

# Fig 2c
rna_read_example <- rna_1m_cells |>
  filter(tool != "cp", format %in% c("10x", "bpcells", "h5ad", "h5ad_zarr")) |>
  mutate(format = normalize_format_name[format]) |>
  group_by(replicate, format) |>
  summarise(read=min(read_seconds_cpu), write=min(write_seconds_cpu), size=unique(size), .groups="drop") |>
  ggplot(aes(format, read, fill=format)) +
  style_storage_barplot +
  labs(y="Read time (seconds)", subtitle="RNA read (1M cells)")
ggsave(file.path(plot_dir, "rna_read_example.svg"), rna_read_example, width=2.8, height=2.8, units="in")

# Fig 2d
rna_write_example <- rna_1m_cells |>
  filter(tool != "cp", format %in% c("bpcells", "h5ad", "h5ad_zarr")) |>
  mutate(format = normalize_format_name[format]) |>
  group_by(replicate, format) |>
  summarise(read=min(read), write=min(write_seconds_cpu, na.rm = TRUE), size=unique(size), .groups="drop") |>
  ggplot(aes(format, write, fill=format)) +
  style_storage_barplot[-3] + scale_y_continuous(expand=expansion(c(0, 0.15)), breaks=scales::breaks_width(5)) +
  labs(y="Write time (seconds)", subtitle="RNA write (1M cells)")
ggsave(file.path(plot_dir, "rna_write_example.svg"), rna_write_example, width=2.8, height=2.8, units="in")

################################################################################
# Fragment file read + write
################################################################################
tool_colors <- c("BPCells"=set1[1], "ArchR"=set1[2], "10x"=set1[3], "GNU sort"=set1[2], "SnapATAC2"=set1[4], "GNU sort"=set1[2])
normalize_name <- factor(c(
  bpcells="BPCells", archr="ArchR", "snapatac2"="SnapATAC2", "10x"="10x", "unix_sort"="GNU sort"
), levels=c("10x", "ArchR", "GNU sort", "SnapATAC2", "BPCells"))

style_fragment_barplot <- list(
  geom_bar = geom_bar(position=position_dodge2(preserve = "single"), stat="summary", fun="mean"),
  geom_point = geom_point(position=position_dodge2(width=0.9), key_glyph="blank"),
  scale_y = scale_y_continuous(expand=expansion(c(0, 0.05))),
  scale_x = scale_x_discrete(guide=guide_axis(angle=30)),
  scale_fill = scale_fill_manual(values=tool_colors),
  theme = theme_classic(base_size=14, base_family = "Arial"),
  guides = guides(fill="none"),
  labs = labs(x=NULL)
)

# Fig 2h
fragment_read_example <- fragment_read |> 
  filter(dataset == "1m_brain") |>
  group_by(dataset, tool, replicate) |>
  summarize(time_cpu=sum(time_cpu)) |>
  mutate(tool=normalize_name[tool]) |>
  ggplot(aes(tool, time_cpu/60, fill=tool)) +
  style_fragment_barplot +
  labs(y="Read (minutes)", subtitle="ATAC read (1M cells)")
ggsave(file.path(plot_dir, "fragment_read_example.svg"), fragment_read_example, width=2.8, height=2.8, units="in")

# Fig 2i
fragment_write_example <- fragment_import |>
  filter(dataset == "1m_brain", tool != "10x") |>
  group_by(dataset, tool, replicate) |>
  summarize(time_cpu=sum(time_cpu)) |>
  mutate(tool=normalize_name[tool]) |>
  ggplot(aes(tool, time_cpu/(60*60), fill=tool)) +
  style_fragment_barplot +
  scale_y_continuous(expand=expansion(c(0, 0.05)), breaks=scales::breaks_pretty()) +
  labs(y="Import (hours)", subtitle="ATAC import (1M cells)")
ggsave(file.path(plot_dir, "fragment_write_example.svg"), fragment_write_example, width=2.8, height=2.8, units="in")

# Fig 2g
fragment_size_example <- fragment_import |>
  filter(dataset == "1m_brain") |>
  group_by(dataset, tool, replicate) |>
  summarize(bytes=sum(bytes), .groups="drop_last") |>
  summarize(bytes=unique(bytes), .groups="drop") |>
  mutate(tool=normalize_name[tool]) |>
  ggplot(aes(tool, bytes/1e9, fill=tool)) +
  style_fragment_barplot[!(names(style_fragment_barplot) %in% "geom_point")] +
  labs(y="Size (GB)", subtitle="ATAC storage (1M cells)")
ggsave(file.path(plot_dir, "fragment_size_example.svg"), fragment_size_example, width=2.8, height=2.8, units="in")

################################################################################
# Comparison labels for figures
################################################################################
# Fig 2j ATAC in-memory benchmark: BPCells 3.9x faster read, 8.7x faster write
in_memory_summary %>%
  filter(filter %in% c("na", "byte"), codec %in% c("zstd-blosc", "BPCells"), level %in% c(2,0), sample=="10k_pbmc") %>%
  summarize(relative_read = read_speed[codec=="BPCells"]/read_speed[codec != "BPCells"],
            relative_write = write_speed[codec=="BPCells"]/write_speed[codec != "BPCells"])

# Fig 2e RNA feature-major in-memory benchmark: BPCells 4.4x faster read, 7.7x faster write
in_memory_summary %>%
  filter(filter %in% c("na", "byte"), codec %in% c("zstd-blosc", "BPCells"), level %in% c(2,0), sample=="130k_thymus_transpose") %>%
  summarize(relative_read = read_speed[codec=="BPCells"]/read_speed[codec != "BPCells"],
            relative_write = write_speed[codec=="BPCells"]/write_speed[codec != "BPCells"])

# Fig S3a RNA cell-major in-memory benchmark: BPCells 2.9x faster read, 4.1x faster write than lz4-blosc
# (Fastest algorithm with better compression ratio than BPCells)
in_memory_summary %>%
  filter(filter %in% c("na", "byte"), sample=="130k_thymus") |>
  group_by(sample) |>
  mutate(write_slowdown = write_speed[codec=="BPCells"]/write_speed, read_slowdown = read_speed[codec=="BPCells"]/read_speed, size_savings=1 - compression_ratio[codec=="BPCells"]/compression_ratio) |>
  filter(compression_ratio > compression_ratio[codec=="BPCells"]) |>
  filter(write_speed == max(write_speed) | read_speed == max(read_speed))

# Fig 2c + 2d BPCells vs zarr: 3x faster read, 2x faster write
rna_1m_cells |>
  group_by(tool, format) |>
  summarize(size=unique(size), read=mean(read), write=mean(write)) |>
  filter(format %in% c("bpcells", "h5ad_zarr") | (format=="h5ad" & tool=="python")) |>
  ungroup() |>
  mutate(relative_read=read/min(read), relative_write=write/min(write))

# Fig 2h: BPCells 5.0x faster than ArchR for reading
fragment_read |> 
  filter(dataset == "1m_brain") |>
  group_by(dataset, tool, replicate) |>
  summarize(time_cpu=sum(time_cpu), .groups="drop") |>
  mutate(tool=normalize_name[tool]) |>
  group_by(dataset, tool) |>
  summarize(read_time = mean(time_cpu)/60, .groups="drop") |>
  mutate(relative_read_time = read_time/read_time[tool=="BPCells"])

# Fig 2i: BPCells 5.4x faster than SnapATAC2 for import
fragment_import |>
  filter(dataset == "1m_brain") |>
  group_by(dataset, tool, replicate) |>
  summarize(time_cpu=sum(time_cpu), bytes=sum(bytes), .groups="drop") |>
  mutate(tool=normalize_name[tool]) |>
  group_by(dataset, tool) |>
  summarize(write_time = mean(time_cpu)/60, size_gb=unique(bytes)/1e9, .groups="drop") |>
  mutate(write_relative_time=write_time/write_time[tool=="BPCells"])


################################################################################
# Specific results for text
################################################################################

# 10x hdf5 is 7x smaller than h5ad
rna_1m_cells |> filter(format %in% c("10x", "h5ad")) |> distinct(format, size) |>
  mutate(relative_size = size/min(size))

# 10x hdf5 takes over 6x longer to read from disk
rna_1m_cells |>
  filter(tool != "cp", format %in% c("10x", "h5ad")) |>
  group_by(replicate, format) |>
  summarise(read=min(read_seconds_cpu), write=min(write_seconds_cpu), size=unique(size), .groups="drop") |>
  group_by(format) |>
  summarize(read = mean(read)) |>
  mutate(relative_read = read/min(read))

# On 1M RNA dataset, BPCells 4.9GB storage, second only to 10x at 4.2
rna_1m_cells |> distinct(format, size) |>
  mutate(size = size/1e9) |>
  arrange(desc(size))

# 10x 73 seconds to read, zarr 20, uncompressed h5ad 11, bpcells 6.8
# zarr 23 seconds write, h5ad 13 seconds, BPCells 12
rna_1m_cells |>
  filter(tool != "cp") |>
  group_by(replicate, format) |>
  summarise(read=min(read_seconds_cpu), write=min(write_seconds_cpu, na.rm=TRUE), size=unique(size), .groups="drop") |>
  group_by(format) |>
  summarize(read = mean(read), write=mean(write))

# BPCells 1m fragments use 21GB of space, 28% less than next-closest method
fragment_import |>
  filter(dataset == "1m_brain") |>
  group_by(dataset, tool, replicate) |>
  summarize(bytes=sum(bytes), .groups="drop_last") |>
  summarize(size_gb=unique(bytes)/1e9, .groups="drop") |>
  mutate(bpcells_percent_savings = 100 * (size_gb - size_gb[tool=="bpcells"])/size_gb) |>
  arrange(size_gb)

# BPCells 1m fragments read 5x faster than next-closest alternative at 1.2 minutes
fragment_read |> 
  filter(dataset == "1m_brain") |>
  group_by(dataset, tool, replicate) |>
  summarize(time_cpu=sum(time_cpu), .groups="drop") |>
  group_by(dataset, tool) |>
  summarize(read_time = mean(time_cpu)/60, .groups="drop") |>
  mutate(relative_read_time = read_time/read_time[tool=="bpcells"])

# BPCells 1m fragments 5.4x faster import at 40 minutes
fragment_import |>
  filter(dataset == "1m_brain") |>
  group_by(dataset, tool, replicate) |>
  summarize(time_cpu=sum(time_cpu), bytes=sum(bytes), .groups="drop") |>
  group_by(dataset, tool) |>
  summarize(write_time = mean(time_cpu)/60, size_gb=unique(bytes)/1e9, .groups="drop") |>
  mutate(write_relative_time=write_time/write_time[tool=="bpcells"])


# 226 unique compression settings in the supplemental data 
in_memory_summary |> filter(!(codec %in% c("np.copyto", "memcpy", "BPCells"))) |> 
  distinct(codec, filter, level) |>
  nrow()

# 82 unique compression settings on the plot (not counting BPCells)
in_memory_summary %>%
  filter(filter %in% c("na", "byte"), !(codec %in% c("np.copyto", "memcpy", "BPCells")), sample=="10k_pbmc") |> 
  distinct(codec, filter, level) |>
  nrow()

# BPCells reaad speed 7.4-8.6 GB/s, write speed 5.1-5.6
in_memory_summary %>%
  filter(codec == "BPCells")

# Maximal space savings:
#  - feature-major: 27% space savings,  102,500% slowdown
#  - cell-major: 42% space savings, 304,500% slowdown
#  - fragments:  9.6% space savings, 68,900% slowdown
in_memory %>%
  mutate(
    filter=replace_na(filter, "na"),
    level=replace_na(level, 0)
  ) %>%
  group_by(sample, codec, runner, filter, level, replicate) %>%
  summarize(
    write_min=sum(write_min),
    bytes=sum(bytes),
    input_gb=sum(input_bytes)/1e9,
    fields=n(),
    .groups="drop"
  ) %>%
  group_by(sample, codec, runner, filter, level) %>%
  summarize(bytes = unique(bytes), input_gbp = unique(input_gb), write=mean(write_min), .groups="drop") %>%
  group_by(sample) %>%
  mutate(bpcells_bytes = bytes[codec=="bpcells"], space_savings = (bpcells_bytes - bytes)/bpcells_bytes, relative_write=write/write[codec=="bpcells"]) %>%
  slice_min(bytes, n=1)

