library(phangorn)
library(homoplasyFinder)

PROTEIN <- "Spike"
SITESMAPPING_FILE <- "Data/sitesMapping.csv"

# HOMOPLASYFINDER_RES_DIR <-
#     "Data/nextstrain_homoplasyFinder_results/"
# HOMOPLASYFINDER_RES_FILE <- "Data/nextstrain_homoplasyFinder.csv"
# # HOMOPLASYFINDER_RES_AA_FILE <-
# #     paste0("Data/nextstrain_homoplasyFinder_", PROTEIN, ".csv")

HOMOPLASYFINDER_RES_DIR <- "Data/sampled_homoplasyFinder_results/"
HOMOPLASYFINDER_RES_FILE <- "Data/sampled_homoplasyFinder.csv"
# HOMOPLASYFINDER_RES_AA_FILE <-
#     paste0("Data/sampled_homoplasyFinder_", PROTEIN, ".csv")

CURRENT_DIR <- getwd()

sitesMapping <- read.csv(SITESMAPPING_FILE, row.names = 1)

res <- data.frame(
    "site" = integer(),
    "consistencyIndex" = double(),
    "date" = character()
)
comparisonRes <- logical()

for (fn in list.files(HOMOPLASYFINDER_RES_DIR)) {
    print(fn)
    workingDirectory <-
        file.path(CURRENT_DIR, HOMOPLASYFINDER_RES_DIR, fn)
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
    
    # tryCatch({
    #     setwd(workingDirectory)
    #     unlink(
    #         list.files(
    #             ".",
    #             "sequences_noInconsistentSites_|annotatedNewickTree_|consistencyIndexReport_*"
    #         )
    #     )
    #     invisible(runHomoplasyFinderJarTool(treeFile, fastaFile))
    #
    #     results2 <- read.table(
    #         list.files(".", "consistencyIndexReport_*"),
    #         header = TRUE,
    #         sep = "\t",
    #         stringsAsFactors = FALSE
    #     )
    # },
    # finally = setwd(CURRENT_DIR))
    # test <-
    #     merge(results, results2, by.x = "site", by.y = "Position")
    # comparisonRes[fn] <-
    #     identical(test$ConsistencyIndex, test$consistencyIndex)
}

write.csv(res, HOMOPLASYFINDER_RES_FILE)

# print(all(comparisonRes))


# res_aa <- data.frame(
#     "site" = integer(),
#     "consistencyIndex" = double(),
#     "protein" = character(),
#     "date" = character()
# )
#
# for (fn in list.files(HOMOPLASYFINDER_RES_DIR)) {
#     print(fn)
#     workingDirectory <-
#         file.path(CURRENT_DIR, HOMOPLASYFINDER_RES_DIR, fn)
#     f_base_name <- file.path(workingDirectory,
#                              paste0(fn, "_", PROTEIN))
#     fastaFile <- paste0(f_base_name, ".fasta")
#     treeFile <- paste0(f_base_name, ".newick")
#
#     lvls <- tolower(sitePath:::AA_SHORT_NAMES)
#
#     seqs <- read.phyDat(
#         fastaFile,
#         format = "fasta",
#         type = "USER",
#         levels = lvls
#     )
#     # seqs <- read.phyDat(fastaFile,
#     #                     format = "fasta",
#     #                     type = "AA")
#     tree <- read.tree(treeFile)
#     consistencyIndex <- CI(tree, seqs, sitewise = TRUE)
#     print(length(consistencyIndex))
#     sites <- which(!is.na(consistencyIndex) & consistencyIndex < 1)
#     results <- data.frame(
#         "site" = sites,
#         "consistencyIndex" = consistencyIndex[sites],
#         "protein" = PROTEIN,
#         "date" = fn
#     )
#     res_aa <- rbind(res_aa, results)
#     write.csv(results, paste0(f_base_name, "_consistency_index.csv"))
#
# }
#
# write.csv(res, HOMOPLASYFINDER_RES_AA_FILE)
