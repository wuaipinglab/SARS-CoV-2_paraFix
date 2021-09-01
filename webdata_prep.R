library(sitePath)

DATA_DIR <- file.path("webdata", "sitePath_results")

options(list(
  "cl.cores" = 5
))

for (fn in list.files(DATA_DIR)) {
  print(fn)
  minEntropy <- readRDS(file.path(DATA_DIR, fn, paste0(fn, ".rds")))
  
  fixed <- fixationSites(minEntropy)
  saveRDS(fixed, file = file.path(DATA_DIR, fn, "fixation.rds"))
  
  para <- parallelSites(minEntropy)
  saveRDS(para, file = file.path(DATA_DIR, fn, "parallel.rds"))

  paraFix <- paraFixSites(minEntropy, mutMode = "all")
  saveRDS(paraFix, file = file.path(DATA_DIR, fn, "paraFix.rds"))
  
}
