library(jsonlite)
library(ggplot2)

PROTEIN_NAME <- "Spike"
SITES_PREVALENCE_FILE <- "Data/sitesPrevalence.json"
BACKGROUND_NUM_FILE <- "Data/background_num.json"
MUTATION_NUM_FILE <- "Data/mutation_num.json"
PREVALENCE_INTO_FILE <- "Data/prevalenceInfo.csv"

BASELINE_FILE <- "Data/nextstrain_trees_unused.csv"
THRESHOLD_FILE <- "Data/nextstrain_trees_unused_0.1.csv"
PARAFIXSITES_FILE <- "Data/nextstrain_sitePath_results.csv"
HOMOPLASYFINDER_RES_FILE <- "Data/nextstrain_homoplasyFinder.csv"
HYPHY_RES_FILE <- "Data/nextstrain_hyphy_results.csv"
CONSERVED_SITES_FILE <- "Data/nextstrain_conserved_sites.csv"
DATES_FILE <- "Data/nextstrain_dates.json"
METRIC_PLOT_1 <- "Output/nextstrain_specificity.pdf"
METRIC_PLOT_2 <- "Output/nextstrain_sensitivity.pdf"

# BASELINE_FILE <- "Data/sampled_trees_unused.csv"
# THRESHOLD_FILE <- "Data/sampled_trees_unused_0.1.csv"
# PARAFIXSITES_FILE <- "Data/sampled_sitePath_results.csv"
# HOMOPLASYFINDER_RES_FILE <- "Data/sampled_homoplasyFinder.csv"
# HYPHY_RES_FILE <- "Data/sampled_hyphy_results.csv"
# CONSERVED_SITES_FILE <- "Data/sampled_conserved_sites.csv"
# DATES_FILE <- "Data/sampled_dates.json"
# METRIC_PLOT_1 <- "Output/sampled_specificity.pdf"
# METRIC_PLOT_2 <- "Output/sampled_sensitivity.pdf"



background <- read_json(BACKGROUND_NUM_FILE, simplifyVector = TRUE)
mutationsNum <- read_json(MUTATION_NUM_FILE, simplifyVector = TRUE)
allDates <- read_json(DATES_FILE, simplifyVector = TRUE)

sitesPrevalence <-
    read_json(SITES_PREVALENCE_FILE, simplifyVector = TRUE)

conserved <- read.csv(CONSERVED_SITES_FILE)
conserved <- conserved[which(conserved$product == PROTEIN_NAME), ]
conserved <- split(conserved$aaPos, conserved$date)
# allMutSites <- read.csv(PREVALENCE_INTO_FILE)

allPrevalentSites <- integer()
nonPrevalentSites <- integer()
for (prevalency in sitesPrevalence) {
    allPrevalentSites <- c(allPrevalentSites, prevalency[["prevalent"]])
    nonPrevalentSites <-
        c(nonPrevalentSites, prevalency[["non_prevalent"]])
}
allPrevalentSites <- unique(allPrevalentSites)
nonPrevalentSites <- unique(nonPrevalentSites)


# names(allPrevalentSites) <- allPrevalentSites
# allPrevalentSites <- sapply(allPrevalentSites, function(site) {
#     site
# })

baselineSites <- read.csv(BASELINE_FILE)
baselineSites <-
    baselineSites[which(baselineSites$product == PROTEIN_NAME), ]
baselineSites <-
    split(baselineSites[, c("aaPos", "mutRatio")], baselineSites$date)

thresholdSites <- read.csv(THRESHOLD_FILE)
thresholdSites <-
    thresholdSites[which(thresholdSites$product == PROTEIN_NAME), ]
thresholdSites <-
    split(thresholdSites[, c("aaPos", "mutRatio")], thresholdSites$date)


paraFixSites <- read.csv(PARAFIXSITES_FILE)
paraFixSites <-
    paraFixSites[which(paraFixSites$product == PROTEIN_NAME), ]
paraFixSites <-
    split(paraFixSites[, c("aaPos", "type")], paraFixSites$date)


hyphySites <- read.csv(HYPHY_RES_FILE)
hyphySites <-
    split(hyphySites[, c("site", "posterior")], hyphySites$date)


homoplasySites <- read.csv(HOMOPLASYFINDER_RES_FILE)
homoplasySites <-
    homoplasySites[which(homoplasySites$product == PROTEIN_NAME), ]
homoplasySites <-
    split(homoplasySites[, c("aaPos", "consistencyIndex")], homoplasySites$date)


getMetrics <- function(pred,
                       true_positive,
                       true_negative,
                       test_date,
                       method_name) {
    tp <- length(intersect(pred, true_positive))
    fp <- length(setdiff(pred, true_positive))
    
    tn <- length(setdiff(true_negative, pred))
    fn <- length(setdiff(true_positive, pred))
    return(
        data.frame(
            "specificity" = tn / (tn + fp),
            "sensitivity" = tp / (tp + fn),
            "method" = method_name,
            "date" = as.Date(test_date)
        )
    )
}


allMetrics <- do.call(rbind, lapply(allDates, function(d) {
    true_negative <- setdiff(nonPrevalentSites, conserved[[d]])
    bs_res <- getMetrics(baselineSites[[d]][["aaPos"]],
                         allPrevalentSites,
                         true_negative,
                         d,
                         "MSA")
    th_res <- getMetrics(thresholdSites[[d]][["aaPos"]],
                         allPrevalentSites,
                         true_negative,
                         d,
                         "MSA_0.1")
    sp_res <- getMetrics(paraFixSites[[d]][["aaPos"]],
                         allPrevalentSites,
                         true_negative,
                         d,
                         "sitePath")
    hp_res <- getMetrics(hyphySites[[d]][["site"]],
                         allPrevalentSites,
                         true_negative,
                         d,
                         "hyphy_fubar")
    pa_res <- getMetrics(unique(homoplasySites[[d]][["aaPos"]]),
                         allPrevalentSites,
                         true_negative,
                         d,
                         "homoplasy")
    do.call(rbind, list(bs_res, th_res, sp_res, hp_res, pa_res))
}))


p <- ggplot(allMetrics,
            aes(
                x = date,
                y = specificity,
                group = method,
                color = method
            )) +
    geom_line() + geom_point() + ggtitle("specificity")
p

ggsave(
    filename = METRIC_PLOT_1,
    plot = p,
    device = "pdf",
    width = 6,
    height = 6
)

p <- ggplot(allMetrics,
            aes(
                x = date,
                y = sensitivity,
                group = method,
                color = method
            )) +
    geom_line() + geom_point() + ggtitle("sensitivity")
p
ggsave(
    filename = METRIC_PLOT_2,
    plot = p,
    device = "pdf",
    width = 6,
    height = 6
)
