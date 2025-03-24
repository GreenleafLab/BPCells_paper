library(tidyverse)
library(patchwork)
library(legendry)

################################################################################
# Data loading 
################################################################################

repo_root <- "/home/bparks/Sync/BPCells_paper/final_upload_version/"

results_dir <- file.path(repo_root, "/results/data_tables")
plot_dir <- file.path(repo_root, "/results/plots/bitwidth-stats")
dir.create(plot_dir, showWarnings = FALSE)

matrix_bitwidth <- read_tsv(file.path(results_dir, "compression/bitwidth-stats/matrix_bitwidth.tsv.gz"))
matrix_storage <- read_tsv(file.path(results_dir, "compression/bitwidth-stats/matrix_storage.tsv"))
matrix_histogram <- read_tsv(file.path(results_dir, "compression/bitwidth-stats/matrix_value_histogram.tsv.gz"))
matrix_shape <- read_tsv(file.path(results_dir, "compression/bitwidth-stats/matrix_shape.tsv"))

fragment_bitwidth <- read_tsv(file.path(results_dir, "compression/bitwidth-stats/fragment_bitwidth.tsv.gz"))
fragment_storage <- read_tsv(file.path(results_dir, "compression/bitwidth-stats/fragment_storage.tsv"))

atac_samples <- read_tsv(file.path(results_dir, "datasets/atac.tsv")) |>
  filter(sample != "merged") |>
  group_by(dataset, sample) |>
  summarize(
    cells = barcodes[filtered],
    barcodes=barcodes[!filtered],
    filtered_fragments = total_fragments[filtered],
    total_fragments = total_fragments[!filtered]
  )
atac_datasets <- read_tsv(file.path(results_dir, "datasets/atac.tsv")) |>
  group_by(dataset) |>
  mutate(samples=n_distinct(sample[sample!="merged"])) |>
  filter(sample == "merged") |>
  summarize(
    cells = barcodes[filtered],
    barcodes = barcodes[!filtered],
    filtered_fragments = total_fragments[filtered],
    total_fragments = total_fragments[!filtered],
    median_fragments = median_fragments[filtered]
  ) |>
  arrange(cells)

rna_datasets <- read_tsv(file.path(results_dir, "datasets/rna.tsv"))


################################################################################
# Plot styling helpers
################################################################################
colors_tab10 <- palette.colors(palette="Tableau 10")
field_colors <- colors_tab10[c(1:5, 7)]
names(field_colors) <- c("index", "cell", "value", "total", "end", "start")

################################################################################
# ATAC fragment storage
################################################################################

# Fig S2h: bits per fragment within 1m_brain samples
total_bits_fragments <- fragment_bitwidth %>%
  group_by(dataset, sample, field, filtered) %>%
  summarize(bitwidth = weighted.mean(bitwidth, chunk_count)) %>%
  group_by(dataset, sample, filtered) %>%
  summarize(bitwidth = sum(bitwidth)) %>%
  mutate("field" = "total") %>%
  inner_join(atac_samples, by=c("dataset", "sample"))

atac_bits_plot_brain <- fragment_bitwidth %>%
  group_by(dataset, sample, field, filtered) %>%
  summarize(bitwidth = weighted.mean(bitwidth, chunk_count), .groups="drop") %>%
  inner_join(atac_samples, by=c("dataset", "sample")) %>%
  filter(dataset == "1m_brain", filtered==TRUE) %>%
  ggplot(aes(x=cells, y=bitwidth, color=field, fill=field)) +
  scale_color_manual(values = field_colors, breaks = c("cell", "start", "end", "total")) + 
  scale_fill_manual(values = field_colors, breaks = c("cell", "start", "end", "total")) +
  guides(color = guide_legend(byrow = TRUE)) +
  geom_point(key_glyph=draw_key_rect) +
  geom_point(data=filter(total_bits_fragments, dataset=="1m_brain", filtered==TRUE), key_glyph=draw_key_rect) + 
  scale_x_continuous(transform="log10", labels=scales::label_number(scale_cut=scales::cut_short_scale())) +
  labs(y="Mean bits per fragment", x="cells", subtitle="ATAC fragment storage") +
  theme_classic(base_size = 14, base_family="Arial") +
  theme(legend.title = element_blank(), legend.key.spacing.y = unit(2, "pt"))

ggsave(file.path(plot_dir, "atac_bits_plot_brain.svg"), atac_bits_plot_brain, width=4.5, height=4)


# Fig S2i: Bar plot of bits per fragment for different field usage
atac_filter_bar_plot <- fragment_storage %>%
  mutate(merged = sample == "merged") %>%
  group_by(dataset, merged, filtered, field) %>%
  summarize(bytes=sum(bytes), elements=sum(elements)) %>%
  mutate(
    filtered = as.factor(if_else(filtered, "filtered", "all")),
    merged = as.factor(if_else(merged, "merged", "separate")),
    label = factor(paste0(filtered, ".", merged), 
                   levels=c("all.separate", "filtered.separate", "all.merged", "filtered.merged")),
    field = recode_factor(field, cell="cell", start="start", end="end")
  ) %>%
  filter(dataset == "35k_hematopoiesis", field %in% c("cell", "end", "start")) %>%
  ggplot(aes(label, bytes*8/elements, fill=field)) +
  geom_col() +
  scale_fill_manual(values=field_colors) +
  scale_y_continuous(expand=expansion(c(0, NA)), limits = c(0,NA)) +
  labs(y = "bits per fragment", x="", subtitle="Fragment storage\n(35k hematopoiesis)") +
  guides(x = legendry::guide_axis_nested(), fill="none") +
  theme_classic(base_size = 14, base_family="Arial") +
  theme(legend.title = element_blank(), axis.text.x = element_text(margin=margin(b=0.5)))
ggsave(file.path(plot_dir, "atac_filter_bar_plot.svg"), atac_filter_bar_plot, width=4, height=4.25)

# Fig S2g: Distribution of start coordinate bit widths in 10k PBMC dataset
atac_start_distribution <- fragment_bitwidth %>%
  filter(dataset == "10k_pbmc", field == "start", chunk_count > 0) %>%
  group_by(field) %>%
  mutate(chunk_fraction = chunk_count / sum(chunk_count)) %>%
  ggplot(aes(bitwidth, chunk_fraction)) +
  geom_col(fill=field_colors["start"]) +
  scale_x_continuous(expand=expansion(0), limits = c(0,12)) +
  scale_y_continuous(expand=expansion(c(0, 0.05)), labels = scales::label_percent()) +
  labs(x = "Bit width", y="Fraction of chunks", subtitle="Start coordinates (10k PBMC)") +
  theme_classic(base_size = 14, base_family="Arial")
ggsave(file.path(plot_dir, "atac_start_distribution.svg"), atac_start_distribution, width=4, height=4)


################################################################################
# RNA & ATAC matrix storage order
################################################################################

# Fig S2b: Bar plot showing effect of RNA matrix ordering on storage size (1M neurons)
rna_ordering_1m_neuron <- matrix_storage %>%
  filter(matrix == "rna", dataset == "1m_neurons", field %in% c("val", "index"), ordering %in% c("original", "mean")) %>%
  mutate(
    cell_major = if_else(cell_major, "cell-major", "feature-major"),
    field = dplyr::recode(field, val="value"),
    dataset = factor(dataset, levels=rev(rna_datasets$dataset)),
    ordering = factor(ordering, levels=c("original", "mean")),
    label = factor(
      paste0(if_else(ordering == "mean", "sorted", "raw"), ".", cell_major),
      levels=c("raw.cell-major", "sorted.cell-major", "raw.feature-major", "sorted.feature-major")
    ),
    compression_ratio = bytes*8 / elements
  ) %>%
  ggplot(aes(label, compression_ratio, shape=cell_major, line_type=ordering, fill=field)) +
  geom_col() +
  scale_y_continuous(breaks=scales::breaks_pretty(), expand = expansion(mult=c(0,0.05))) +
  scale_fill_manual(values=field_colors) +
  labs(y="bits per non-zero entry", x="", subtitle="RNA storage order (1M Neurons)") +
  guides(x = legendry::guide_axis_nested(key=key_range_auto(sep="\\."))) +
  theme_classic(base_size = 14, base_family="Arial") +
  theme(legend.position = "inside", legend.position.inside = c(1,1), legend.justification = c(1,1), legend.title = element_blank()) +
  theme(axis.text.x = element_text(margin=margin(b=0.5)), legend.key.spacing.y = unit(1, "pt"))

ggsave(file.path(plot_dir, "rna_ordering_1m_neuron.svg"), rna_ordering_1m_neuron, width=4, height=4.25)


# Fig S2d: Bar plot of space usage with ATAC peak matrix
atac_peak_ordering <- matrix_storage %>%
  filter(matrix=="peaks", dataset=="35k_hematopoiesis", 
         field %in% c("val", "index"), ordering %in% c("original", "mean")) %>%
  mutate(
    cell_major = if_else(cell_major, "cell-major", "feature-major"),
    field = dplyr::recode(field, val="value"),
    ordering = factor(ordering, levels=c("original", "mean")),
    label = factor(
      paste0(if_else(ordering == "mean", "sorted", "raw"), ".", cell_major),
      levels=c("raw.cell-major", "sorted.cell-major", "raw.feature-major", "sorted.feature-major")
    ),    compression_ratio = bytes*8 / elements
  ) %>%
  ggplot(aes(label, compression_ratio, shape=cell_major, line_type=ordering, fill=field)) +
  geom_col() +
  scale_y_continuous(breaks=scales::breaks_pretty(), expand = expansion(mult=c(0,0)), limits=c(0, 15)) +
  scale_fill_manual(values=field_colors) +
  labs(y="bits per non-zero entry", x="", subtitle="ATAC peak matrix\n(35k hematopoiesis)") +
  guides(x = legendry::guide_axis_nested(key=key_range_auto(sep="\\."))) +
  theme_classic(base_size = 14, base_family="Arial") +
  theme(legend.position = "inside", legend.position.inside = c(1,1.15), legend.justification = c(1,1), legend.title = element_blank()) +
  theme(axis.text.x = element_text(margin=margin(b=0.5)), legend.key.spacing.y = unit(1, "pt"))

ggsave(file.path(plot_dir, "atac_peak_ordering_35k_hematopoiesis.svg"), atac_peak_ordering, width=4, height=4.25)




# Fig S2a: Boxplot of compression ratio by storage order (RNA)
rna_ordering_boxplot <- matrix_storage %>%
  filter(matrix=="rna", field %in% c("val", "index"), ordering %in% c("original", "mean")) %>%
  group_by(dataset, cell_major, ordering, elements) %>%
  summarize(bytes=sum(bytes), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(
    cell_major = if_else(cell_major, "cell-major", "feature-major"),
    dataset = factor(dataset, levels=rev(rna_datasets$dataset)),
    ordering = factor(ordering, levels=c("original", "mean")),
    label = factor(
      paste0(if_else(ordering == "mean", "sorted", "raw"), ".", cell_major),
      levels=c("raw.cell-major", "sorted.cell-major", "raw.feature-major", "sorted.feature-major")
    ),
    compression_ratio = elements*8 / bytes
  ) %>%
  ggplot(aes(label, compression_ratio)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_dodge2(width=0.2)) + 
  scale_y_continuous(breaks=scales::breaks_pretty()) +
  scale_color_brewer(palette="Set1") +
  guides(color="none", shape="none", x=legendry::guide_axis_nested(key=key_range_auto(sep="\\."))) +
  labs(x=" ", y="Compression Ratio", subtitle="RNA Counts Matrix Compression") +
  theme_classic(base_size = 14, base_family="Arial") +
  theme(legend.title = element_blank(), axis.text.x = element_text(margin=margin(b=0.5)))
ggsave(file.path(plot_dir, "rna_ordering_boxplot.svg"), rna_ordering_boxplot, width=4, height=4.25)


# Fig S2c: Boxplot of compression ratio by storage order (ATAC matrices)
atac_ordering_boxplot <- matrix_storage %>%
  inner_join(atac_datasets, by="dataset") %>%
  filter(field %in% c("val", "index"), ordering %in% c("original", "mean"), 
         elements>0,
        matrix %in% c("peaks", "tiles")) %>%
  group_by(dataset, cell_major, ordering, matrix, elements) %>%
  summarize(bytes=sum(bytes), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(
    cell_major = if_else(cell_major, "cell-major", "feature-major"),
    dataset = factor(dataset, levels=rev(atac_datasets$dataset)),
    ordering = factor(ordering, levels=c("original", "mean")),
    label = factor(
      paste0(if_else(ordering == "mean", "sorted", "raw"), ".", cell_major),
      levels=c("raw.cell-major", "sorted.cell-major", "raw.feature-major", "sorted.feature-major")
    ),
    compression_ratio = elements*8 / bytes
  ) %>%
  ggplot(aes(label, compression_ratio, color=matrix, fill=matrix)) +
  geom_boxplot(outlier.shape = NA, mapping = aes(fill=NULL), key_glyph=draw_key_rect) +
  geom_point(position=ggplot2::position_jitterdodge(jitter.width=0.1, jitter.height=0, seed=125124), key_glyph=draw_key_rect) +
  scale_y_continuous(breaks=scales::breaks_pretty(3), limits=c(NA, 6.4)) +
  scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
  labs(x="", y="Compression Ratio", subtitle="ATAC Counts Matrix Compression", color="", fill="") +
  guides(color = guide_legend(byrow = TRUE), x=legendry::guide_axis_nested(key=key_range_auto(sep="\\."))) +
  theme_classic(base_size = 14, base_family="Arial") +
  theme(legend.position="inside", legend.position.inside = c(0.01,1), legend.justification = c(0,.75)) +
  theme(axis.text.x = element_text(margin=margin(b=0.5)), legend.key.spacing.y = unit(2, "pt"))

ggsave(file.path(plot_dir, "atac_ordering_boxplot.svg"), atac_ordering_boxplot, width=4, height=4.25)

################################################################################
# Matrix storage scaling relationships
################################################################################

# Fig S2f: Dot plot of matrix value field storage scaling
matrix_p99 <- matrix_histogram %>%
  group_by(dataset, matrix) %>%
  arrange(value, count) %>%
  summarize(p99 = which(cumsum(count/sum(count)) <= 0.99) |> tail(1))

matrix_value_scaling <- matrix_storage %>%
  inner_join(matrix_p99, by=c("dataset", "matrix")) %>%
  filter(field == "val", cell_major==TRUE, ordering=="original") %>%
  mutate(
    cell_major = if_else(cell_major, "cell-major", "feature-major"),
    matrix=recode_factor(matrix, rna="RNA", peaks="peaks", tiles="tiles")
  ) %>%
  ggplot(aes(log2(p99), bytes/elements * 8, shape=matrix)) +
  geom_point(color=field_colors["value"]) +
  scale_shape_discrete(solid=FALSE) +
  scale_x_continuous(labels=scales::label_math(2^.x)) +
  scale_y_continuous(breaks=1:8) +
  theme_classic(base_size = 7, base_family="Arial") +
  theme(legend.position = "inside", legend.position.inside = c(1,0.01), legend.justification = c(1,0)) +
  theme(legend.key.size = unit(8, "points")) +
  labs(x="99th percentile non-zero value", y="Mean bits per entry", subtitle="Matrix value storage", shape="")
ggsave(file.path(plot_dir, "matrix_scaling_value.svg"), matrix_value_scaling, width=2, height=2)



# Fig S2e: Dot plot of matrix index field storage scaling
matrix_index_scaling <- matrix_storage |>
  inner_join(matrix_shape) |>
  filter(field == "index", cell_major==TRUE, ordering=="original") %>%
  mutate(
    cell_major = if_else(cell_major, "cell-major", "feature-major"),
    matrix=recode_factor(matrix, rna="RNA", peaks="peaks", tiles="tiles"),
    fraction_nonzero = elements / (rows*cols)
  ) %>%
  ggplot(aes(fraction_nonzero, bytes/elements * 8, shape=matrix)) +
  geom_point(color=field_colors["index"]) +
  scale_shape_discrete(solid=FALSE) +
  scale_x_continuous(transform = "log10", guide="axis_logticks", labels=scales::label_log()) +
  scale_y_continuous(breaks=scales::breaks_pretty(4)) +
  theme_classic(base_size = 7, base_family="Arial") +
  theme(legend.position="inside", legend.position.inside = c(1,1), legend.justification = c(1,1)) +
  theme(legend.key.size = unit(8, "points")) +
  labs(x="Fraction non-zero", y="Mean bits per entry", subtitle="Matrix index storage", shape="")
ggsave(file.path(plot_dir, "matrix_scaling_index.svg"), matrix_index_scaling, width=2, height=2)



################################################################################
# Specific results for text
################################################################################

# The majority of non-zero values are 1 in an RNA-seq matrix
matrix_histogram |>
  filter(matrix == "rna") |>
  group_by(dataset) |>
  summarize(count[value==1]/sum(count))

# Typically more than 99% of values in RNA matrix are less than 50
matrix_histogram %>%
  filter(matrix == "rna") %>% 
  group_by(dataset, matrix) %>%
  arrange(value, count) %>%
  summarize(p99 = value[which(cumsum(count/sum(count)) > 0.99)[1]])

# RNA space savings of ~4x for benchmark datasets in cell-major, 6-7x when using different ordering
matrix_storage %>%
  filter(matrix=="rna", field %in% c("val", "index"), ordering %in% c("original", "mean")) %>%
  group_by(dataset, cell_major, ordering, elements) %>%
  summarize(bytes=sum(bytes), .groups = "drop") %>%
  mutate(compression_ratio = elements*8 / bytes) %>%
  group_by(cell_major, ordering) |>
  summarize(median_compression_ratio = median(compression_ratio))

# For ATAC peak and tile matrices 4-6x space savings
matrix_storage %>%
  inner_join(atac_datasets, by="dataset") %>%
  filter(field %in% c("val", "index"), ordering %in% c("original", "mean"), 
         elements>0,
         matrix %in% c("peaks", "tiles")) %>%
  group_by(dataset, cell_major, ordering, matrix, elements) %>%
  summarize(bytes=sum(bytes), .groups = "drop") %>%
  mutate(compression_ratio = elements*8 / bytes) %>%
  group_by(cell_major, ordering, matrix) %>%
  summarize(median_compression_ratio = median(compression_ratio))


# Filtering out non-cell barcodes reduces per-fragment storage size by ~18%
fragment_storage %>%
  mutate(merged = sample == "merged") %>%
  filter(dataset=="35k_hematopoiesis") %>%
  group_by(dataset, merged, filtered, field) %>%
  summarize(bytes=sum(bytes), elements=sum(elements), .groups = "drop_last") %>%
  summarize(bits_per_frag=sum(bytes)*8/unique(elements)) %>%
  group_by(dataset, merged) %>%
  summarize(
    bits_all = bits_per_frag[!filtered],
    bits_filt = bits_per_frag[filtered],
    size_reduction = (bits_all - bits_filt)/bits_all
  )
  
# Size of the 22m_pansci dataset is 30GB  
matrix_storage |> 
  filter(dataset=="22m_pansci", ordering=="original", cell_major) |> 
  pull(bytes) |>
  sum()

# Methods: Matrix storage metadata typically <2% of total bytes
matrix_storage |> group_by(dataset, matrix, cell_major, ordering) |>
  summarize(frac_metadata=bytes[field=="metadata"]/sum(bytes)) |>
  group_by(dataset, matrix) |>
  summarize(min_metadata_frac=min(frac_metadata), max_metadata_frac=max(frac_metadata)) |>
  arrange(max_metadata_frac)

# Methods: Fragment metadata is at worst 5%, typically under 1%
fragment_storage |>
  group_by(dataset, filtered) |>
  summarize(frac_metadata = sum(bytes[field %in% c("end_max", "metadata")])/sum(bytes))

