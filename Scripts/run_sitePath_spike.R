library(jsonlite)
library(ggplot2)
library(sitePath)


PROTEIN_NAME <- "Spike"
SITESMAPPING_FILE <- "Data/sitesMapping.csv"

TREES_DIR <- "Data/nextstrain_trees_with_MSA/"
SITEPATH_RES_DIR <- "Data/nextstrain_sitePath_results/"
PARAFIXSITES_FILE <-
    paste0("Data/nextstrain_sitePath_results_", PROTEIN_NAME, ".csv")
DATES_FILE <- "Data/nextstrain_dates.json"
OUTPUT_DIR <- "Data/SARS-CoV-2_spike/"

dir.create(OUTPUT_DIR, showWarnings = FALSE)


allDates <- read_json(DATES_FILE, simplifyVector = TRUE)
sitesMapping <- read.csv(SITESMAPPING_FILE, row.names = 1)
targetInfo <-
    sitesMapping[which(sitesMapping$product == PROTEIN_NAME),]
targetInfo[["dnaPos"]] <- seq_len(nrow(targetInfo))
targetSites <- unique(targetInfo$peptidePos)


options(list("cl.cores" = 20))

fn <- as.character(max(as.Date(allDates)))

outFile <- file.path(OUTPUT_DIR, "minEntropy.rds")
if (!file.exists(outFile)) {
    inputDir <- file.path(TREES_DIR, fn)
    tree <- read.tree(file.path(inputDir, paste0(fn, ".nwk")))
    paths <- addMSA(
        tree = tree,
        msaPath = file.path(inputDir, paste0(fn, "_aa.fasta")),
        msaFormat = "fasta"
    )
    paths <- setSiteNumbering(paths, usedSites = targetSites)
    minEntropy <- sitesMinEntropy(paths)
    saveRDS(minEntropy, file = outFile)
} else {
    minEntropy <- readRDS(outFile)
}

outFile <- file.path(OUTPUT_DIR, paste0("fixation.rds"))
if (!file.exists(outFile)) {
    fixed <- fixationSites(minEntropy)
    saveRDS(fixed, file = outFile)
}

outFile <- file.path(OUTPUT_DIR, paste0("parallel.rds"))
if (!file.exists(outFile)) {
    para <- parallelSites(minEntropy)
    saveRDS(para, file = outFile)
}

res <- list("fixation" = as.integer(allSitesName(fixed)),
            "parallel" = as.integer(allSitesName(para)))

write_json(res, "Data/SARS-CoV-2_spike/targetSites.json")


p <- plot(fixed) + theme(legend.position = "none")
print(p)
ggsave(
    filename = "Output/sitePath_fixation.svg",
    plot = p,
    device = "svg",
    width = 10,
    height = 6
)

p <- plot(para, y = F)
print(p)
ggsave(
    filename = "Output/sitePath_parallel.svg",
    plot = p,
    device = "svg",
    width = 10,
    height = 6
)
