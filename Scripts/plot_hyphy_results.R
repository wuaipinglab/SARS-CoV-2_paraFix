library(jsonlite)

PROTEIN_NAME <- "Spike"
HYPHY_DIR <- "Data/sampled_hyphy_results/"
HYPHY_RES_FILE <- "Data/sampled_hyphy_results.csv"

res <- data.frame(
    "site" = integer(),
    "posterior" = double(),
    "protein" = character(),
    "date" = character()
)

for (fn in list.files(HYPHY_DIR)) {
    fp <-
        file.path(HYPHY_DIR, fn, paste0(fn, "_", PROTEIN_NAME, ".nexus.FUBAR.json"))
    hyphyFubar <- read_json(fp)
    resolved <- hyphyFubar[["MLE"]][["content"]][[1]]
    names(resolved) <- seq_along(resolved)
    
    contentNames <- sapply(hyphyFubar[["MLE"]][["headers"]], "[[", 1)
    resolved <- lapply(resolved, function(site) {
        names(site) <- contentNames
        site
    })
    selectedSites <- resolved[which(sapply(resolved, function(site) {
        if (site[["Prob[alpha<beta]"]] >= 0.9) {
            return(TRUE)
        }
        return(FALSE)
    }))]
    
    selectedSites <-
        do.call(rbind, lapply(names(selectedSites), function(siteName) {
            site <- selectedSites[[siteName]]
            posterior <- site[["Prob[alpha<beta]"]]
            data.frame(
                "site" = as.integer(siteName),
                "posterior" = posterior,
                "protein" = PROTEIN_NAME,
                "date" = fn
            )
        }))
    res <- rbind(res, selectedSites)
}

write.csv(res, HYPHY_RES_FILE)


# hyphySlac <- read_json("Data/sampled_hyphy_results/2021-07-01/2021-07-01_Spike.nexus.SLAC.json")
#
# resolved <- hyphySlac[["MLE"]][["content"]][[1]][["by-site"]][["RESOLVED"]]
# names(resolved) <- seq_along(resolved)
#
# contentNames <- sapply(hyphySlac[["MLE"]][["headers"]], "[[", 1)
# selectedSites <- which(sapply(resolved, function(site) {
#     if (site[[3]] + site[[4]] != 0) {
#         return(site[[9]] < 0.05 | site[[10]] < 0.05)
#     }
#     return(FALSE)
# }))
# selectedSites <- resolved[selectedSites]
#
# selectedSites <- do.call(rbind, lapply(names(selectedSites), function(siteName) {
#     site <- selectedSites[[siteName]]
#     pValues <- c(site[[9]], site[[10]])
#     i <- which(pValues < 0.05)
#     data.frame(
#         "site" = as.integer(siteName),
#         "test" = contentNames[9:10][[i]],
#         "pValue" = pValues[[i]]
#     )
# }))
#
# selectedSites <- selectedSites[which(selectedSites$test == "P [dN/dS > 1]"), ]
