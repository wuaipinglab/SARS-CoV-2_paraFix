library(jsonlite)
library(ggplot2)

PROTEIN_NAME <- "Spike"
# PROTEIN_NAME <- "N"

SAMPLING_METHOD <- "nextstrain"
# SAMPLING_METHOD <- "sampled"

BACKGROUND_NUM_FILE <- "Data/background_num.json"

SITES_PREVALENCE_FILE <-
    paste0("Data/sitesPrevalence_", PROTEIN_NAME, ".json")
MUTATION_NUM_FILE <-
    paste0("Data/mutation_num_", PROTEIN_NAME, ".json")
PREVALENCE_INTO_FILE <-
    paste0("Data/prevalenceInfo_", PROTEIN_NAME, ".csv")


BASELINE_FILE <-
    paste0("Data/", SAMPLING_METHOD, "_trees_unused.csv")
THRESHOLD_FILE <-
    paste0("Data/", SAMPLING_METHOD, "_trees_unused_0.1.csv")
PARAFIXSITES_FILE <-
    paste0("Data/", SAMPLING_METHOD, "_sitePath_results.csv")
HOMOPLASYFINDER_RES_FILE <-
    paste0("Data/", SAMPLING_METHOD, "_homoplasyFinder.csv")
HYPHY_RES_FILE <-
    paste0("Data/", SAMPLING_METHOD, "_hyphy_results.csv")
CONSERVED_SITES_FILE <-
    paste0("Data/", SAMPLING_METHOD, "_conserved_sites.csv")
DATES_FILE <- paste0("Data/", SAMPLING_METHOD, "_dates.json")

PLOT_TITLE <- paste0(SAMPLING_METHOD, "_", PROTEIN_NAME)
METRIC_PLOT <- paste0("Output/", PLOT_TITLE, ".pdf")



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
    return(data.frame(
        "metric" = c("specificity", "sensitivity"),
        "value" = c(tn / (tn + fp), tp / (tp + fn)),
        "method" = method_name,
        "date" = as.Date(test_date)
    ))
}


allMetrics <- do.call(rbind, lapply(allDates, function(d) {
    true_positive <- setdiff(allPrevalentSites, conserved[[d]])
    true_negative <- setdiff(nonPrevalentSites, conserved[[d]])
    # bs_res <- getMetrics(baselineSites[[d]][["aaPos"]],
    #                      true_positive,
    #                      true_negative,
    #                      d,
    #                      "MSA")
    th_res <- getMetrics(thresholdSites[[d]][["aaPos"]],
                         true_positive,
                         true_negative,
                         d,
                         "MSA_only")
    sp_res <- getMetrics(paraFixSites[[d]][["aaPos"]],
                         true_positive,
                         true_negative,
                         d,
                         "sitePath")
    hp_res <- getMetrics(hyphySites[[d]][["site"]],
                         true_positive,
                         true_negative,
                         d,
                         "hyphy_fubar")
    pa_res <- getMetrics(unique(homoplasySites[[d]][["aaPos"]]),
                         true_positive,
                         true_negative,
                         d,
                         "homoplasy")
    do.call(rbind, list(th_res, sp_res, hp_res, pa_res))
}))


p <- ggplot(allMetrics,
            aes(
                x = date,
                y = value,
                group = method,
                color = method
            )) +
    geom_line() +
    geom_point() +
    facet_grid(cols = vars(metric)) +
    ylim(0, 1) +
    ggtitle(PLOT_TITLE)
p

ggsave(
    filename = METRIC_PLOT,
    plot = p,
    device = "pdf",
    width = 10,
    height = 6
)
