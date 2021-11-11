library(jsonlite)
library(seqinr)

# SAMPLING_METHOD <- "nextstrain"
SAMPLING_METHOD <- "sampled"

MUTATION_THRESHOLD <- 0.1

TREES_DIR <- paste0("Data/", SAMPLING_METHOD, "_trees_with_MSA/")
DATES_FILE <- paste0("Data/", SAMPLING_METHOD, "_dates.json")
CONSERVED_SITES_FILE <-
    paste0("Data/", SAMPLING_METHOD, "_conserved_sites.csv")
BASELINE_FILE <-
    paste0("Data/", SAMPLING_METHOD, "_trees_unused.csv")
THRESHOLD_FILE <- paste0("Data/",
                         SAMPLING_METHOD,
                         "_trees_unused_",
                         MUTATION_THRESHOLD,
                         ".csv")

# TREES_DIR <- "Data/sampled_trees_with_MSA/"
# DATES_FILE <- "Data/sampled_dates.json"
# CONSERVED_SITES_FILE <- "Data/sampled_conserved_sites.csv"
# BASELINE_FILE <- "Data/sampled_trees_unused.csv"
# THRESHOLD_FILE <- "Data/sampled_trees_unused_0.1.csv"

sitesMapping <- read.csv("Data/sitesMapping.csv", row.names = 1)
allDates <- read_json(DATES_FILE, simplifyVector = TRUE)

sitesMappingAA <-
    unique(sitesMapping[, c("product", "peptidePos", "aa", "aaPos")])

mutRatioSummary <- lapply(allDates, function(fn) {
    print(fn)
    
    seqFile <- file.path(TREES_DIR, fn, paste0(fn, "_aa.fasta"))
    seqs <- read.alignment(seqFile, "fasta")
    seqs <- toupper(seqs$seq)
    nsites <- nchar(seqs[[1]])
    mutRatios <-
        do.call(rbind, lapply(seq_len(nsites), function(i) {
            aaSummary <- table(substr(seqs, i, i))
            refInfo <-
                sitesMappingAA[which(sitesMappingAA$peptidePos == i), ]
            refAA <- refInfo[["aa"]]
            
            refNum <- aaSummary[[refAA]]
            mutNum <-
                sum(aaSummary[which(names(aaSummary) != refAA &
                                        names(aaSummary) %in% sitePath:::AA_UNAMBIGUOUS)])
            data.frame(
                "product" = refInfo[["product"]],
                "aaPos" = refInfo[["aaPos"]],
                "mutRatio" = mutNum / refNum,
                "date" = fn
            )
        }))
})

res <- do.call(rbind, lapply(mutRatioSummary, function(mutRatios) {
    mutRatios[which(mutRatios$mutRatio > MUTATION_THRESHOLD), ]
}))

write.csv(res, THRESHOLD_FILE, row.names = FALSE)

res2 <- do.call(rbind, lapply(mutRatioSummary, function(mutRatios) {
    mutRatios[which(mutRatios$mutRatio > 0), ]
}))
write.csv(res2, BASELINE_FILE, row.names = FALSE)


conserved <-
    do.call(rbind, lapply(mutRatioSummary, function(mutRatios) {
        mutRatios[which(mutRatios$mutRatio == 0),]
    }))
write.csv(conserved, CONSERVED_SITES_FILE, row.names = FALSE)
