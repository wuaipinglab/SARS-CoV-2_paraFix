library(jsonlite)

PROTEIN_NAME <- "Spike"
SITES_PREVALENCE_FILE <- "Data/sitesPrevalence.json"
BACKGROUND_NUM_FILE <- "Data/background_num.json"
MUTATION_NUM_FILE <- "Data/mutation_num.json"


# BASELINE_FILE <- "Data/nextstrain_trees_unused_results.csv"
# PARAFIXSITES_FILE <- "Data/nextstrain_sitePath_results.csv"
# HOMOPLASYFINDER_RES_FILE <- "Data/nextstrain_homoplasyFinder.csv"
# HYPHY_RES_FILE <- "Data/nextstrain_hyphy_results.csv"
# DATES_FILE <- "Data/nextstrain_dates.json"

BASELINE_FILE <- "Data/sampled_trees_unused_results.csv"
PARAFIXSITES_FILE <- "Data/sampled_sitePath_results.csv"
HOMOPLASYFINDER_RES_FILE <- "Data/sampled_homoplasyFinder.csv"
HYPHY_RES_FILE <- "Data/sampled_hyphy_results.csv"
DATES_FILE <- "Data/sampled_dates.json"


sitesPrevalence <- read_json(SITES_PREVALENCE_FILE, simplifyVector = TRUE)
background <- read_json(BACKGROUND_NUM_FILE, simplifyVector = TRUE)
mutationsNum <- read_json(MUTATION_NUM_FILE, simplifyVector = TRUE)
allDates <- read_json(DATES_FILE, simplifyVector = TRUE)


baselineSites <- read.csv(BASELINE_FILE)
baselineSites <- baselineSites[which(baselineSites$product == PROTEIN_NAME), ]
baselineSites <- split(baselineSites[, c("aaPos", "mutRatio")], baselineSites$date)


paraFixSites <- read.csv(PARAFIXSITES_FILE)
paraFixSites <- paraFixSites[which(paraFixSites$product == PROTEIN_NAME), ]
paraFixSites <- split(paraFixSites[, c("aaPos", "type")], paraFixSites$date)


hyphySites <- read.csv(HYPHY_RES_FILE)
hyphySites <- split(hyphySites[, c("site", "posterior")], hyphySites$date)


homoplasySites <- read.csv(HOMOPLASYFINDER_RES_FILE)
homoplasySites <- homoplasySites[which(homoplasySites$product == PROTEIN_NAME), ]
homoplasySites <- split(homoplasySites[, c("aaPos", "consistencyIndex")], homoplasySites$date)


for (d in allDates) {
    bs_res <- baselineSites[[d]][["aaPos"]]
    sp_res <- paraFixSites[[d]][["aaPos"]]
    hp_res <- hyphySites[[d]][["site"]]
    pa_res <- unique(homoplasySites[[d]][["aaPos"]])
    break
}
