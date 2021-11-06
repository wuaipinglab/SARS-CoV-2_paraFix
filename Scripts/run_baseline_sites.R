library(jsonlite)
library(seqinr)

SITESMAPPING_FILE <- "Data/sitesMapping.csv"
MUTATION_THRESHOLD <- 0.1

# TREES_DIR <- "Data/nextstrain_trees_with_MSA/"
# BASELINE_FILE <- "Data/nextstrain_trees_unused_results.csv"
# DATES_FILE <- "Data/nextstrain_dates.json"

TREES_DIR <- "Data/sampled_trees_with_MSA/"
BASELINE_FILE <- "Data/sampled_trees_unused_results.csv"
DATES_FILE <- "Data/sampled_dates.json"

sitesMapping <- read.csv(SITESMAPPING_FILE, row.names = 1)
allDates <- read_json(DATES_FILE, simplifyVector = TRUE)

sitesMappingAA <-
    unique(sitesMapping[, c("product", "peptidePos", "aa", "aaPos")])

res <- lapply(allDates, function(fn) {
    print(fn)
    
    seqFile <- file.path(TREES_DIR, fn, paste0(fn, "_aa.fasta"))
    seqs <- read.alignment(seqFile, "fasta")
    seqs <- toupper(seqs$seq)
    nsites <- nchar(seqs[[1]])
    mutRatios <- do.call(rbind, lapply(seq_len(nsites), function(i) {
        aaSummary <- table(substr(seqs, i, i))
        refInfo <-
            sitesMappingAA[which(sitesMappingAA$peptidePos == i), ]
        refAA <- refInfo[["aa"]]
        
        refNum <- aaSummary[[refAA]]
        mutNum <- sum(aaSummary[which(names(aaSummary) != refAA &
                                          names(aaSummary) %in% sitePath:::AA_UNAMBIGUOUS)])
        data.frame(
            "product" = refInfo[["product"]],
            "aaPos" = refInfo[["aaPos"]],
            "mutRatio" = mutNum / refNum,
            "date" = fn
        )
    }))
    mutRatios[which(mutRatios$mutRatio > MUTATION_THRESHOLD), ]
})

res <- do.call(rbind, res)

write.csv(res, BASELINE_FILE, row.names = FALSE)
