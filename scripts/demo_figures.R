library(sitePath)
library(ggplot2)
suppressPackageStartupMessages(library(ggtree))

output_plot_dir <- file.path("plots", "demo_figures")
tree_height_ratio <- 1.5

if (dir.exists(output_plot_dir)) {
    unlink(output_plot_dir, recursive = TRUE)
}
dir.create(output_plot_dir, showWarnings = FALSE)

plot_path_comparison <- function(paths,
                                 select1,
                                 select2,
                                 path_size = NULL,
                                 branch_size = 1) {
    if (is.null(path_size)) {
        path_size <- 6
    }
    path_nodes <- unique(unlist(paths[c(select1, select2)]))

    group <- rep(1, times = n_nodes)
    group[path_nodes] <- 0
    group <- factor(group)

    group_colors <- rep(2, times = n_nodes)
    group_colors[path_nodes] <- 1
    group_colors[unlist(paths[select1])] <- 0
    group_colors <- factor(group_colors)

    size <- rep(1, times = n_nodes)
    size[path_nodes] <- 2

    p <- ggtree(tree, aes(
        color = group_colors,
        linetype = group,
        size = size
    )) +
        scale_size(range = c(branch_size, path_size)) +
        scale_color_manual(values = c("red", "blue", "gainsboro")) +
        theme(legend.position = "none")
    p
}

plot_sub_paths <- function(paths,
                           select = NULL,
                           path_color = "blue",
                           path_size = 2,
                           branch_size = 1) {
    tree <- ape::as.phylo(paths)
    if (is.null(select)) {
        path_nodes <- unique(unlist(paths))
    } else {
        path_nodes <- unique(unlist(paths[select]))
    }
    group <- rep(1, times = n_nodes)
    group[path_nodes] <- 0
    group <- factor(group)
    size <- rep(1, times = n_nodes)
    size[path_nodes] <- 2
    p <- ggtree(tree, aes(color = group, linetype = group, size = size)) +
        scale_size(range = c(branch_size, path_size)) +
        scale_color_manual(values = c(path_color, "gainsboro")) +
        theme(legend.position = "none")
    p
}

#===============================================================================

data(h3n2_align)
data(h3n2_tree)

paths <- addMSA(h3n2_tree, alignment = h3n2_align)

tree <- as.phylo(paths)
n_nodes <- ape::Nnode(tree, internal.only = FALSE)
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

# =======================================================================

p_lineages <- plot_sub_paths(paths, path_color = "blue", branch_size = 0.2)

ggsave(
    filename = file.path(output_plot_dir, "lineages.svg"),
    plot = p_lineages,
    device = "svg",
    width = 3,
    height = 4.5
)

path_select <- list("0" = NULL, "1" = 1, "2" = c(1, 2), "3" = 2, "4" = 5)

for (n in names(path_select)) {
    p <- plot_sub_paths(
        paths,
        select = path_select[[n]],
        path_color = "blue",
        path_size = 1,
        branch_size = 0.05
    )
    ggsave(
        filename = file.path(output_plot_dir, paste0("pathway", n, ".svg")),
        plot = p,
        device = "svg",
        width = 1,
        height = 1.3
    )
}

# =======================================================================

select1 <- c(2, 4)

p <- plot_path_comparison(
    paths,
    select1,
    seq_along(paths)[-select1],
    path_size = 1.5,
    branch_size = 0.5
)
p

ggsave(
    filename = file.path(output_plot_dir, "pathComparison0.svg"),
    plot = p,
    device = "svg",
    width = 3.5,
    height = 3.5 * tree_height_ratio
)

n <- 0

for (select2 in seq_along(paths)[-select1]) {
    n <- n + 1
    p <- plot_path_comparison(paths, select1, select2, 1, 0.05)
    ggsave(
        filename = file.path(
            output_plot_dir,
            paste0("pathComparison", n, ".svg")
        ),
        plot = p,
        device = "svg",
        width = 1,
        height = 1 * tree_height_ratio
    )
}

# =======================================================================

all_path_nodes <- table(unlist(paths))
most_recent_nodes <- as.integer(names(which(all_path_nodes == 1)))

path_nodes <- unique(unlist(paths))

group <- rep(1, times = n_nodes)
group[path_nodes] <- 0
group <- factor(group)

group_colors <- rep(2, times = n_nodes)
group_colors[path_nodes] <- 1
group_colors[most_recent_nodes] <- 0
group_colors <- factor(group_colors)

size <- rep(1, times = n_nodes)
size[path_nodes] <- 2

p_parallelity <- ggtree(tree, aes(
    color = group_colors,
    linetype = group,
    size = size
)) +
    scale_size(range = c(0.5, 1.5)) +
    scale_color_manual(values = c("green", "black", "gainsboro")) +
    theme(legend.position = "none")

ggsave(
    filename = file.path(output_plot_dir, "parallelity.svg"),
    plot = p_parallelity,
    device = "svg",
    width = 3.5,
    height = 3.5 * tree_height_ratio
)
