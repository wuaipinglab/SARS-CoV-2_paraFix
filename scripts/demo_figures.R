library(sitePath)
library(ggplot2)
suppressPackageStartupMessages(library(ggtree))

output_plot_dir <- file.path("plots", "demo_figures")
tree_height_ratio <- 1.5

if (dir.exists(output_plot_dir)) {
    unlink(output_plot_dir, recursive = TRUE)
}
dir.create(output_plot_dir, showWarnings = FALSE)

#===============================================================================

data(h3n2_align)
data(h3n2_tree)

paths <- addMSA(h3n2_tree, alignment = h3n2_align)
paths <- lineagePath(paths, 0.04)

min_entropy <- sitesMinEntropy(paths)
fixed_sites <- fixationSites(min_entropy)
para_sites <- parallelSites(min_entropy, minSNP = 10, mutMode = "exact")

#===============================================================================

p_paths <- plotMutSites(paths)

ggsave(
    filename = file.path(output_plot_dir, "paths.svg"),
    plot = p_paths,
    device = "svg",
    width = 5,
    height = 5
)

p_single_fixation <- plot(extractSite(fixed_sites, 211)) +
    ggtitle(label = NULL, subtitle = NULL) +
    theme(legend.position = "none")

ggsave(
    filename = file.path(output_plot_dir, "singleFixation.svg"),
    plot = p_single_fixation,
    device = "svg",
    width = 1.5,
    height = 2
)

p_parallel_mutation <- plotSingleSite(paths, 280, showPath = FALSE) +
    ggtitle(label = NULL, subtitle = NULL) +
    theme(legend.position = "none")

ggsave(
    filename = file.path(output_plot_dir, "parallelMutation.svg"),
    plot = p_parallel_mutation,
    device = "svg",
    width = 1.5,
    height = 2
)

p_mutations <- plotSingleSite(para_sites, 192, showPath = TRUE) +
    ggtitle(label = NULL, subtitle = NULL) +
    theme(legend.position = "none")

ggsave(
    filename = file.path(output_plot_dir, "mutations.svg"),
    plot = p_mutations,
    device = "svg",
    width = 3.5,
    height = 3.5 * tree_height_ratio
)
