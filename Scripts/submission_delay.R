library(ggplot2)

SURVEILLANCE_FILE <- "Data/variant_surveillance.tsv"

df <- read.csv(SURVEILLANCE_FILE, sep = "\t")
df[["Continent"]] <- sapply(strsplit(df$Location, " / "), "[[", 1)


datesInfo <- data.frame(
    "collection" = as.Date(df$Collection.date),
    "submission" = as.Date(df$Submission.date)
)

datesInfo <- datesInfo[which(complete.cases(datesInfo)), ]
datesInfo <-
    datesInfo[which(datesInfo$collection > as.Date("2019-11-30")), ]

submission_delay <-
    as.integer(datesInfo$submission - datesInfo$collection)


hist(x = submission_delay[which(submission_delay < 60 &
                                    submission_delay >= 0)],
     main = "Submission delay",
     xlab = "days")


plot(x = table(df[["Continent"]]),
     main = "Sampling total",
     ylab = "#Sequences")
