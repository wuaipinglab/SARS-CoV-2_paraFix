#!/usr/bin/env Rscript

library(phangorn)

PROTEIN <- "Spike"
SITESMAPPING_FILE <- "data/sitesMapping.csv"

HOMOPLASYFINDER_RES_DIR <- "output/nextstrain_homoplasyFinder_results/"
HOMOPLASYFINDER_RES_FILE <- "output/nextstrain_homoplasyFinder.csv"

sitesMapping <- read.csv(SITESMAPPING_FILE, row.names = 1)

res <- data.frame(
    "site" = integer(),
    "consistencyIndex" = double(),
    "date" = character()
)
comparisonRes <- logical()

for (fn in list.files(HOMOPLASYFINDER_RES_DIR)) {
    print(fn)
    workingDirectory <- file.path(HOMOPLASYFINDER_RES_DIR, fn)
    f_base_name <- file.path(workingDirectory,
                             paste0(fn, "_genome"))
    fastaFile <- paste0(f_base_name, ".fasta")
    treeFile <- paste0(f_base_name, ".newick")
    
    seqs <- read.phyDat(fastaFile,
                        format = "fasta",
                        type = "DNA", )
    tree <- read.tree(treeFile)
    consistencyIndex <- CI(tree, seqs, sitewise = TRUE)
    print(length(consistencyIndex))
    sites <- which(!is.na(consistencyIndex) & consistencyIndex < 1)
    results <- data.frame(
        "site" = sites,
        "consistencyIndex" = consistencyIndex[sites],
        "date" = fn
    )
    
    results <- merge(
        x = results,
        y = sitesMapping,
        by.x = "site",
        by.y = "genomePos",
        all.y = TRUE
    )
    
    results <- results[which(complete.cases(results)), ]
    
    res <- rbind(res, results)
    write.csv(results, paste0(f_base_name, "_consistency_index.csv"))
}

write.csv(res, HOMOPLASYFINDER_RES_FILE, row.names = FALSE)
