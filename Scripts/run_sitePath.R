library(jsonlite)
library(sitePath)


SITESMAPPING_FILE <- "Data/sitesMapping.csv"

# TREES_DIR <- "Data/sampled_trees_with_MSA/"
# SITEPATH_RES_DIR <- "Data/sampled_sitePath_results/"
# PARAFIXSITES_FILE <- "Data/sampled_sitePath_results.csv"
# DATES_FILE <- "Data/sampled_dates.json"

TREES_DIR <- "Data/nextstrain_trees_with_MSA/"
SITEPATH_RES_DIR <- "Data/nextstrain_sitePath_results/"
PARAFIXSITES_FILE <- "Data/nextstrain_sitePath_results.csv"
DATES_FILE <- "Data/nextstrain_dates.json"


cat("Run sitePath for each date\n")

allDates <- read_json(DATES_FILE, simplifyVector = TRUE)

options(list("cl.cores" = 20))

for (fn in allDates) {
    print(fn)
    outputDir <- file.path(SITEPATH_RES_DIR, fn)
    dir.create(outputDir, showWarnings = FALSE)
    
    outFile <- file.path(outputDir, paste0(fn, ".rds"))
    if (!file.exists(outFile)) {
        inputDir <- file.path(TREES_DIR, fn)
        # f_base_name <- paste0("sample_com_", gsub("-", "", fn))
        # tree <-
        #     read.tree(file.path(inputDir, paste0(f_base_name, ".fasta.treefile")))
        # paths <- addMSA(
        #     tree = tree,
        #     msaPath = file.path(inputDir, paste0(f_base_name, "_aa.fasta")),
        #     msaFormat = "fasta"
        # )
        tree <-
            read.tree(file.path(inputDir, paste0(fn, ".nwk")))
        paths <- addMSA(
            tree = tree,
            msaPath = file.path(inputDir, paste0(fn, "_aa.fasta")),
            msaFormat = "fasta"
        )
        minEntropy <- sitesMinEntropy(paths)
        saveRDS(minEntropy, file = outFile)
    } else {
        minEntropy <- readRDS(outFile)
    }
    
    outFile <- file.path(outputDir, "fixation.rds")
    if (!file.exists(outFile)) {
        fixed <- fixationSites(minEntropy)
        saveRDS(fixed, file = outFile)
    }
    
    outFile <- file.path(outputDir, "parallel.rds")
    if (!file.exists(outFile)) {
        para <- parallelSites(minEntropy)
        saveRDS(para, file = outFile)
    }
    
    outFile <- file.path(outputDir, "paraFix.rds")
    if (!file.exists(outFile)) {
        paraFix <- paraFixSites(minEntropy, mutMode = "all")
        saveRDS(paraFix, file = outFile)
    }
}


cat("Summary the results...")

sitesMapping <- read.csv(SITESMAPPING_FILE, row.names = 1)

allSites <- lapply(allDates, function(collectionDate) {
    print(collectionDate)
    
    resDir <- file.path(SITEPATH_RES_DIR, collectionDate)
    fixedSites <- readRDS(file.path(resDir, "fixation.rds"))
    paraSites <- readRDS(file.path(resDir, "parallel.rds"))
    collectionDate <- as.Date(collectionDate)
    
    fixed <- allSitesName(fixedSites)
    para <- allSitesName(paraSites)
    paraFix <- intersect(fixed, para)
    fixedOnly <- setdiff(fixed, paraFix)
    paraOnly <- setdiff(para, paraFix)
    data.frame(
        "site" = as.integer(c(fixedOnly, paraOnly, paraFix)),
        "type" = c(
            rep("fixation", length(fixedOnly)),
            rep("parallel", length(paraOnly)),
            rep("paraFix", length(paraFix))
        ),
        "date" = rep(collectionDate, length(c(
            fixedOnly, paraOnly, paraFix
        )))
    )
    
})

paraFixSites <- merge(
    x = do.call(rbind, allSites),
    y = unique(sitesMapping[, c("gene", "product", "aaPos", "peptidePos", "aa")]),
    by.x = "site",
    by.y = "peptidePos",
    all.x = TRUE
)

write.csv(paraFixSites, PARAFIXSITES_FILE, row.names = FALSE)
