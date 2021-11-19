library(sitePath)


SITESMAPPING_FILE <- "Data/sitesMapping.csv"
TREES_DIR <- "Data/latest_trees_with_MSA/"
PARAFIXSITES_FILE <- "Data/latest_sitePath_results.csv"

options(list("cl.cores" = 20))

sitesMapping <- read.csv(SITESMAPPING_FILE, row.names = 1)

allSites <- data.frame(
    "site" = integer(),
    "type" = character(),
    "date" = vector()
)
for (fn in dir(TREES_DIR)) {
    print(fn)
    outputDir <- file.path(TREES_DIR, fn)
    outFile <- file.path(outputDir, paste0(fn, ".rds"))
    if (!file.exists(outFile)) {
        inputDir <- file.path(TREES_DIR, fn)
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
        fixedSites <- fixationSites(minEntropy)
        saveRDS(fixedSites, file = outFile)
    } else {
        fixedSites <- readRDS(outFile)
    }
    
    outFile <- file.path(outputDir, "parallel.rds")
    if (!file.exists(outFile)) {
        paraSites <- parallelSites(minEntropy)
        saveRDS(paraSites, file = outFile)
    } else {
        paraSites <- readRDS(outFile)
    }
    
    fixed <- allSitesName(fixedSites)
    para <- allSitesName(paraSites)
    paraFix <- intersect(fixed, para)
    fixedOnly <- setdiff(fixed, paraFix)
    paraOnly <- setdiff(para, paraFix)
    res <- data.frame(
        "site" = as.integer(c(fixedOnly, paraOnly, paraFix)),
        "type" = c(
            rep("fixation", length(fixedOnly)),
            rep("parallel", length(paraOnly)),
            rep("paraFix", length(paraFix))
        ),
        "date" = rep(as.Date(fn), length(c(
            fixedOnly, paraOnly, paraFix
        )))
    )
    allSites <- rbind(allSites, res)
}


paraFixSites <- merge(
    x = allSites,
    y = unique(sitesMapping[, c("gene", "product", "aaPos", "peptidePos", "aa")]),
    by.x = "site",
    by.y = "peptidePos",
    all.x = TRUE
)

write.csv(paraFixSites, PARAFIXSITES_FILE, row.names = FALSE)
