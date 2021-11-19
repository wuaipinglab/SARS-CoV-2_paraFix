library(jsonlite)
library(ggplot2)
library(sitePath)


PROTEIN_NAME <- "Spike"
SITESMAPPING_FILE <- "Data/sitesMapping.csv"

TREES_DIR <- "Data/latest_trees_with_MSA/"
OUTPUT_DIR <- "Data/SARS-CoV-2_spike/"
TARGETSITES_FILE <- "Data/SARS-CoV-2_spike/targetSites.json"

dir.create(OUTPUT_DIR, showWarnings = FALSE)


sitesMapping <- read.csv(SITESMAPPING_FILE, row.names = 1)
targetInfo <-
    sitesMapping[which(sitesMapping$product == PROTEIN_NAME),]
targetInfo[["dnaPos"]] <- seq_len(nrow(targetInfo))
targetSites <- unique(targetInfo$peptidePos)

longBranch <- c("EPI_ISL_6125398", "EPI_ISL_2431853")

options(list("cl.cores" = 20))

fn <- as.character(max(as.Date(dir(TREES_DIR))))

outFile <- file.path(OUTPUT_DIR, "minEntropy.rds")
if (!file.exists(outFile)) {
    inputDir <- file.path(TREES_DIR, fn)
    tree <- read.tree(file.path(inputDir, paste0(fn, ".nwk")))
    tree <- ape::drop.tip(tree, longBranch)
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
} else {
    fixed <- readRDS(outFile)
}

outFile <- file.path(OUTPUT_DIR, paste0("parallel.rds"))
if (!file.exists(outFile)) {
    para <- parallelSites(minEntropy)
    saveRDS(para, file = outFile)
} else {
    para <- readRDS(outFile)
}

fixation <- as.integer(allSitesName(fixed))
parallel <- as.integer(allSitesName(para))

res <- list("fixation" = fixation,
            "paraFix" = intersect(fixation, parallel))

write_json(res, TARGETSITES_FILE)
