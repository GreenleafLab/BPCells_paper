library(tidyverse)
library(patchwork)
library(svglite)

################################################################################
# Data loading 
################################################################################
repo_root <- "/home/bparks/Sync/BPCells_paper/final_upload_version/"

results_dir <- file.path(repo_root, "/results/data_tables")
plot_dir <- file.path(repo_root, "/results/plots/atac-timing")
dir.create(plot_dir, showWarnings = FALSE)

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
    
peak_tile_timing <- read_tsv(file.path(results_dir, "atac-timing/peak-tile-timing.tsv.gz"))


fragment_merge <- read_tsv(file.path(results_dir, "atac-timing/merge-fragments.tsv"))

################################################################################
# Plot styling helpers
################################################################################

colors_tab10 <- palette.colors(palette="Tableau 10")
field_colors <- colors_tab10[c(1:5, 7)]
names(field_colors) <- c("index", "cell", "value", "total", "end", "start")

set1 <- RColorBrewer::brewer.pal(4, name="Set1")
tool_colors <- c("BPCells"=set1[1], "ArchR"=set1[2], "10x"=set1[3], "GNU sort"=set1[2], "SnapATAC2"=set1[4], "GNU sort"=set1[2])
normalize_name <- factor(c(
  bpcells="BPCells", archr="ArchR", "snapatac2"="SnapATAC2", "10x"="10x", "unix_sort"="GNU sort"
), levels=c("10x", "ArchR", "GNU sort", "SnapATAC2", "BPCells"))

cell_axis <- atac_datasets |>
  mutate(label=str_split_fixed(dataset, "_", 2)[,1]) %>%
  mutate(label=str_replace_all(label, "m", "M")) %>%
  arrange(filtered_fragments) %>%
  {sec_axis(transform="identity", breaks=.$filtered_fragments, labels=.$label, guide=guide_axis(angle=90))}


# Slight concern of B being confused as 10^12 for European readers, but I find B much clearer than G for US audiences 
label_si_custom <- function(x) {
  ret <- scales::label_number(scale_cut=scales::cut_short_scale())(x)
  return(ret)
}

style_peak_tile <- list(
  stat_summary(fun="mean", geom="line", key_glyph=draw_key_rect),
  stat_summary(fun.min=min, fun.max=max, fun="mean", geom="pointrange", size=0.25),
  scale_x_continuous(transform="log10", guide="axis_logticks", sec.axis=cell_axis, labels = label_si_custom, breaks=c(1e7, 1e8, 1e9, 1e10)), 
  scale_y_continuous(transform="log10", guide="axis_logticks", labels = label_si_custom),
  scale_color_manual(values=tool_colors),
  labs(x = "Total fragments", color="Tool"),
  theme_classic(base_size=14, base_family="Arial"),
  theme(panel.border = element_rect(color="black", linewidth = rel(1), fill=NA)),
  guides(color="none")
)

################################################################################
# Peak & Tile matrix timing
################################################################################

# Fig 1f: peak matrix creation time
peaks_time <- peak_tile_timing |>
  filter(region_type=="peaks-100000") |>
  group_by(dataset, tool, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)) |>
  inner_join(atac_datasets, by="dataset") |>
  mutate(tool=normalize_name[tool]) |>
  ggplot(aes(filtered_fragments, time_cpu, color=tool)) +
  style_peak_tile +
  labs(y="CPU Time (seconds)", subtitle="Peak matrix creation")
ggsave(file.path(plot_dir, "peaks_time.svg"), peaks_time, width=3.375, height=3.75, units="in")

# Fig S1f: peak matrix creation memory
peaks_memory <- peak_tile_timing |>
  filter(region_type=="peaks-100000") |>
  group_by(dataset, tool, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)) |>
  inner_join(atac_datasets, by="dataset") |>
  mutate(tool=normalize_name[tool]) |>
  ggplot(aes(filtered_fragments, max_rss/1e9, color=tool)) +
  style_peak_tile +
  scale_y_continuous(limits = c(0,NA), expand=expansion(c(0,0.05)), labels=scales::label_number(prefix="    ")) +
  labs(y="Memory (GB)",  subtitle="Peak matrix creation")
ggsave(file.path(plot_dir, "peaks_memory.svg"), peaks_memory, width=3.5, height=3.75, units="in")

# Fig 1e: tile matrix creation time
tiles_time <- peak_tile_timing |>
  filter(region_type=="tiles") |>
  group_by(dataset, tool, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)) |>
  inner_join(atac_datasets, by="dataset") |>
  mutate(tool=normalize_name[tool]) |>
  ggplot(aes(filtered_fragments, time_cpu, color=tool)) +
  style_peak_tile +
  labs(y="CPU Time (seconds)", subtitle="Tile matrix creation")
ggsave(file.path(plot_dir, "tiles_time.svg"), tiles_time, width=3.375, height=3.75, units="in")

# Fig S1e: tile matrix creation memory
tiles_memory <- peak_tile_timing |>
  filter(region_type=="tiles") |>
  group_by(dataset, tool, replicate) |>
  summarize(time_cpu=sum(time_cpu), max_rss=max(max_rss)) |>
  inner_join(atac_datasets, by="dataset") |>
  mutate(tool=normalize_name[tool]) |>
  ggplot(aes(filtered_fragments, max_rss/1e9, color=tool)) +
  style_peak_tile +
  scale_y_continuous(limits = c(0,NA), expand=expansion(c(0,0.05)), breaks=scales::breaks_pretty(), labels=scales::label_number(prefix="  ")) +
  labs(y="Memory (GB)", subtitle="Tile matrix creation")
ggsave(file.path(plot_dir, "tiles_memory.svg"), tiles_memory, width=3.5, height=3.75, units="in")

# Fig 1g: Peak subset creatio time
peaks_subset <- peak_tile_timing |>
  filter(str_detect(region_type, "peaks"), dataset=="1m_brain") |>
  group_by(dataset, tool, replicate, peak_count) |>
  summarize(time_cpu=sum(time_cpu)) |>
  inner_join(atac_datasets, by="dataset") |>
  mutate(tool=normalize_name[tool]) |>
  ggplot(aes(peak_count, time_cpu, color=tool)) +
  style_peak_tile[-3] +
  scale_x_continuous(transform="log10", guide="axis_logticks", labels=label_si_custom, breaks= c(10, 1e3, 1e5)) +
  labs(x="Peaks", y="CPU Time (seconds)", subtitle="Peak subsets\n(1M cells)\n")
ggsave(file.path(plot_dir, "peaks_subset.svg"), peaks_subset, width=3.375, height=3.75, units="in")


################################################################################
# Fragment merge
################################################################################

# Fig S3c: Fragment merge time + memory
merge_time <- fragment_merge |>
  inner_join(atac_datasets, by="dataset") |>
  mutate(tool=normalize_name[tool]) |>
  ggplot(aes(filtered_fragments, time_elapsed/60, color=tool)) +
  style_peak_tile + 
  scale_x_continuous(transform="log10", guide="axis_logticks", sec.axis=cell_axis, labels = scales::label_log(), breaks=c(1e7, 1e8, 1e9, 1e10)) +
  labs(y="Elapsed time (minutes)")

merge_memory <- fragment_merge |>
  inner_join(atac_datasets, by="dataset") |>
  mutate(tool=normalize_name[tool]) |>
  ggplot(aes(filtered_fragments, max_rss/1e9, color=tool)) +
  style_peak_tile + 
  scale_y_continuous(limits = c(0,NA), expand=expansion(c(0,0.05)), labels=scales::label_number(prefix="   ")) +
  scale_x_continuous(transform="log10", guide="axis_logticks", sec.axis=cell_axis, labels = scales::label_log(), breaks=c(1e7, 1e8, 1e9, 1e10)) +
  labs(y="Memory (GB)")

merge_combined <- merge_memory + merge_time + patchwork::plot_annotation(title="Fragment merge", theme=theme(plot.title=element_text(size=14, family="Arial", hjust=0.5)))
ggsave(plot=merge_combined, file.path(plot_dir, "merge_combined.svg"), width=7, height=3.75)


################################################################################
# Comparison labels for figures
################################################################################

# Fig 1f SnapATAC2 is 52x slower on 500k cell, 100k peaks
# Fig 1g ArchR is 245x slower for 1M cells, 10 peaks 
peak_tile_timing |>
  filter((region_type=="peaks-100000"&dataset=="500k_heart") | (region_type=="peaks-10"&dataset=="1m_brain")) |>
  group_by(dataset, tool, replicate, peak_count) |>
  summarize(time_cpu=sum(time_cpu), .groups="drop") |>
  group_by(tool, peak_count, dataset) |>
  summarize(time_cpu=mean(time_cpu), .groups="drop") |>
  group_by(peak_count, dataset) |>
  mutate(relative_time = time_cpu/time_cpu[tool=="bpcells"]) |>
  arrange(peak_count, dataset)

################################################################################
# Results for text
################################################################################

# BPCells tile matrix 10x faster than ArchR and 3x faster than SnapATAC2
peak_tile_timing |>
  filter(region_type=="tiles",dataset=="500k_heart") |>
  group_by(dataset, tool, replicate, peak_count) |>
  summarize(time_cpu=sum(time_cpu), .groups="drop") |>
  group_by(tool, peak_count, dataset) |>
  summarize(time_cpu=mean(time_cpu), .groups="drop") |>
  group_by(peak_count, dataset) |>
  mutate(relative_time = time_cpu/time_cpu[tool=="bpcells"]) |>
  arrange(peak_count, dataset)

# BPCells peak matrix 52x and 70x faster than SnapATAC2 and ArchR respectively
peak_tile_timing |>
  filter(region_type=="peaks-100000",dataset=="500k_heart") |>
  group_by(dataset, tool, replicate, peak_count) |>
  summarize(time_cpu=sum(time_cpu), .groups="drop") |>
  group_by(tool, peak_count, dataset) |>
  summarize(time_cpu=mean(time_cpu), .groups="drop") |>
  group_by(peak_count, dataset) |>
  mutate(relative_time = time_cpu/time_cpu[tool=="bpcells"]) |>
  arrange(peak_count, dataset)

# 1M cells, 10 peaks: BPCells 6 seconds, ArchR 25 minutes, SnapATAC2 31 minutes
peak_tile_timing |>
  filter(region_type=="peaks-10",dataset=="1m_brain") |>
  group_by(dataset, tool, replicate, peak_count) |>
  summarize(time_cpu=sum(time_cpu), .groups="drop") |>
  group_by(tool, peak_count, dataset) |>
  summarize(time_cpu=mean(time_cpu), .groups="drop") |>
  group_by(peak_count, dataset) |>
  mutate(relative_time = time_cpu/time_cpu[tool=="bpcells"]) |>
  mutate(time_cpu_minutes = time_cpu/60) |> 
  arrange(peak_count, dataset)

# Fragment merge is geomean 38x faster
fragment_merge |>
  group_by(dataset, tool) |> 
  summarize(time_elapsed=mean(time_elapsed), max_rss=mean(max_rss)) |>
  inner_join(atac_datasets, by="dataset") |>
  mutate(fragments_per_sec = filtered_fragments/time_elapsed) |>
  group_by(dataset) |>
  mutate(relative_time=time_elapsed/time_elapsed[tool=="bpcells"]) |>
  ungroup() |>
  filter(tool != "bpcells") |>
  pull(relative_time) |>
  log() |> mean() |> exp()

# Overlaps per second is 33M for peaks, 9M for tiles
peak_tile_timing |>
  filter(tool=="bpcells", region_type %in% c("tiles", "peaks-100000")) |>
  group_by(region_type) |>
  summarize(overlaps_per_sec=sum(matrix_sum)/sum(time_cpu))


# Latex-formatted dataset summary
atac_citation <- tibble::tribble(
  ~dataset, ~citation, ~link,
  "3k_pbmc", "(\\href{https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0}{10x Genomics})", "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz",
  "10k_pbmc", "(\\href{https://www.10xgenomics.com/datasets/10k-human-pbmcs-atac-v2-chromium-controller-2-standard}{10x Genomics})", "https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz",
  "35k_hematopoiesis", "\\cite{granja2019a}", "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139369&format=file&id=32%2C33%2C34%2C35%2C36%2C37%2C38%2C39%2C40%2C41",
  "45k_pancreas", "(ENCODE \\cite{moore2020a})", "https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&assay_title=snATAC-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=pancreas",
  "500k_heart", "(ENCODE \\cite{moore2020a})", "https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&assay_title=snATAC-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=heart+right+ventricle&biosample_ontology.term_name=heart+left+ventricle&limit=100&biosample_ontology.term_name=left+cardiac+atrium&biosample_ontology.term_name=right+cardiac+atrium&biosample_ontology.term_name=heart&biosample_ontology.term_name=left+ventricle+myocardium+superior&biosample_ontology.term_name=left+ventricle+myocardium+inferior&biosample_ontology.term_name=Right+ventricle+myocardium+superior&biosample_ontology.term_name=Right+ventricle+myocardium+inferior&replicates.library.biosample.life_stage=adult&lab.title=Michael+Snyder%2C+Stanford",
  "1m_brain", "\\cite{li2022a}", "http://catlas.org/catlas_downloads/humanbrain/bedpe/"
)

atac_datasets |> 
  arrange(cells) |>
  inner_join(atac_citation, by="dataset") |>
  mutate(
    filtered_fragments = scales::label_comma()(filtered_fragments),
    total_fragments = scales::label_comma()(total_fragments),
    cells = scales::label_comma()(cells),
    barcodes = scales::label_comma()(barcodes),
    median_fragments = scales::label_comma()(median_fragments),
    dataset=str_c(str_replace_all(dataset, "_", " "), " ", citation),
    link=str_c("\\href{", link, "}{\\LinkSymbol}")
  ) |>
  select(dataset, link, cells, filtered_fragments, median_fragments, barcodes, total_fragments) |>
  write_delim(file.path(plot_dir, "atac_datasets_latex.txt"), delim="&", eol="\\\\\n", quote="none")

