library(ape)
library(jsonlite)
library(sitePath)
suppressPackageStartupMessages(library(trackViewer))

DATA_DIR <- "Data"
PLOTS_DIR <- "Output"

PROTEIN_DOMAIN = list(
    "SARS-CoV-2" = list(
        "spike" = list(
            "RBD" = 319:541,
            "NTD" = 13:303,
            "S2" = 686:1273,
            "FP" = 816:855,
            "HR1" = 920:970,
            "HR2" = 1163:1202
        )
    )
)

SITES_COLOR <- c(
    "fixation" = "gold",
    "parallel" = "lightblue", 
    "paraFix" = "red"
)
DOMAIN_COLOR <- c(
    "Seq" = "#BCBDBD",
    "NTD" = "#9BCAC8",
    "RBD" = "#C6716B",
    "S2" = "#9EC9A1",
    "FP" = "#BAD3E1",
    "HR1" = "#F0D0CE",
    "HR2" = "#F0D0CE"
    
    #     "Seq" = "black",
    #     "NTD" = "blue",
    #     "RBD" = "magenta",
    #     "S2" = "green",
    #     "Fusion peptide" = "orange",
    #     "HR1" = "purple",
    #     "HR2" = "purple"
)


assignInNamespace(
    x = "plotFeatures",
    value = function(feature.splited, LINEH, bottomHeight, label_on_feature=FALSE){
        feature.height <- 0
        for(n in seq_along(feature.splited)){
            this.feature.height <- 
                max(c(feature.splited[[n]]$height/2, 
                      .0001)) + 0.2 * LINEH
            feature.height <- feature.height + this.feature.height
            ##baseline
            grid.lines(x=c(0, 1), y=c(bottomHeight+feature.height, 
                                      bottomHeight+feature.height))
            for(m in seq_along(feature.splited[[n]])){
                this.dat <- feature.splited[[n]][m]
                color <- if(is.list(this.dat$color)) this.dat$color[[1]] else 
                    this.dat$color
                if(length(color)==0) color <- "black"
                fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else 
                    this.dat$fill
                if(length(fill)==0) fill <- "white"
                this.cex <- if(length(this.dat$cex)>0) this.dat$cex[[1]][1] else 1
                lwd <- if(length(this.dat$lwd)>0) this.dat$lwd[[1]][1] else 1
                this.feature.height.m <- 
                    if(length(this.dat$height)>0) 
                        this.dat$height[[1]][1] else 
                            2*this.feature.height
                grid.roundrect(x=start(this.dat)-.1, y=bottomHeight+feature.height, 
                          width=width(this.dat)-.8, 
                          height=this.feature.height.m,
                          just="left", gp=gpar(col=color, fill=fill, lwd=lwd), 
                          default.units = "native", r=unit(0.2, "snpc"))
                if(label_on_feature & !is.na(names(this.dat)[1])){
                  grid.text(x=(start(this.dat)+end(this.dat))/2, 
                            y=bottomHeight+feature.height,
                            just = "centre",
                            label = names(this.dat)[1],
                            gp= gpar(list(cex=this.cex * 
                                            this.feature.height.m/
                                            this.feature.height,
                                          color=color)), 
                            default.units = "native")
                }
            }
            feature.height <- feature.height + this.feature.height
        }
        feature.height
    },
    ns = asNamespace("trackViewer")
)

domain4Lolli <- function(domains, seqLen) {
    domainHeight <- 0.1
    seqHelight <- 0.02
    res <- list()
    
    blockStarts <- c(0)
    blockEnds <- c(seqLen)
    blockNames <- c("Seq")
    blockHeights <- c(seqHelight)

    for (domainName in names(domains)) {
        domainSites <- sort(domains[[domainName]])
        site <- domainSites[1]
        cSites <- c(site, site + 1)

        blockStarts <- c(blockStarts, site)
        blockEnd <- site + 1
        blockNames <- c(blockNames, domainName)
        blockHeights <- c(blockHeights, domainHeight)

        pSites <- cSites

        for (site in domainSites[-1]) {
            cSites <- c(site - 1, site, site + 1)
            if (length(intersect(pSites, cSites)) == 0) {
                blockEnds <- c(blockEnds, blockEnd)
                blockStarts <- c(blockStarts, site)
                blockNames <- c(blockNames, domainName)
                blockHeights <- c(blockHeights, domainHeight)
            }
            pSites <- cSites
            blockEnd <- site + 1
        }
        blockEnds <- c(blockEnds, blockEnd)
    }
    
    res[["blockStarts"]] <- blockStarts
    res[["blockWidths"]] <- blockEnds - blockStarts
    res[["blockNames"]] <- blockNames
    res[["blockHeights"]] <- blockHeights
    
    return(res)
}



virus <- "SARS-CoV-2"
protein <- "spike"

inputDir <- file.path(DATA_DIR, paste0(virus, "_", protein))

# The length of the target protein
maxLen <- 1274
domainSites <- PROTEIN_DOMAIN[[virus]][[protein]]

minEntropy <- readRDS(file.path(inputDir, "minEntropy.rds"))


fixed <- fixationSites(minEntropy)

# Collect the fixation sites without gap
df <- data.frame(
    "prevAA" = character(),
    "refSite" = integer(),
    "fixedAA" = character(),
    "color" = character(),
    "fixationDepth" = numeric()
)
paths <- attr(fixed, "paths")
tree <- attr(paths, "tree")
edgeLengths <- node.depth.edgelength(tree)
for (site in allSitesName(fixed)) {
    aaPos <- as.integer(site)
    for (mp in extractSite(fixed, site)) {
        qualified <- TRUE
        for (i in seq_along(mp)[-1]) {
            tips <- mp[[i]]
            AA <- attr(tips, "AA")
            if (AA == '-') { qualifed <- FALSE }
            mrca <- getMRCA(phy = tree, tip = tips)
            fixationDepth <- edgeLengths[mrca]
            prevAA <- attr(mp[[i - 1]], "AA")
        }
        if (qualified) {
            fullNameAA <- sitePath:::AA_FULL_NAMES[[tolower(AA)]]
            newRow <- data.frame(
                # Map the site to the reference
                "prevAA" = prevAA,
                "refSite" = aaPos,
                "fixedAA" = AA,
                "mutName" = paste0(prevAA, aaPos, AA),
                "color" = sitePath:::AA_COLORS[[fullNameAA]],
                "fixationDepth" = fixationDepth
            )
            df <- rbind(df, newRow)
        }
    }
}

df <- do.call(rbind, lapply(split(df, df[["mutName"]]), function(mut) {
    mut[which.min(mut[["fixationDepth"]]), ]
}))


# Prepare for the lolliplot
seqname <- paste("SARS-CoV-2", protein, sep = "_")
gr <- GRanges(seqname, IRanges(df[["refSite"]],
                               width = 1,
                               names = paste0(df[["prevAA"]], df[["refSite"]], df[["fixedAA"]])))
gr$color <- as.character(df[["color"]])
gr$score <- (df[["fixationDepth"]] - min(df[["fixationDepth"]])) * 5000
gr$label <- as.character(df[["fixedAA"]])
gr$label.col <- "#DBD9DA"

# Prepare the domains
domainInfo <- domain4Lolli(domainSites, maxLen)
features <- GRanges(seqname, IRanges(
    "start" = domainInfo$blockStarts,
    "width" = domainInfo$blockWidths,
    "names" = domainInfo$blockNames
))
features$fill <- DOMAIN_COLOR[domainInfo$blockNames]
features$color <- 0
features$height <- domainInfo$blockHeights

pdf(
    file = file.path(PLOTS_DIR, paste0(protein, "_lolliplot_fixed.pdf")),
    width = 16,
    height = 4,
    bg = "transparent"
)
lolliplot(
    SNP.gr = gr,
    features = features,
    ranges = GRanges(seqname, IRanges(0, maxLen)),
    ylab = "Distance to root",
    yaxis = FALSE
)
invisible(dev.off())

lolliplot(
    SNP.gr = gr,
    features = features,
    ranges = GRanges(seqname, IRanges(0, maxLen)),
    ylab = "Distance to root",
    yaxis = FALSE
)


paraSites <- parallelSites(minEntropy)

# Collect the parallel sites without gap
df <- data.frame(
    "prevAA" = character(),
    "refSite" = integer(),
    "fixedAA" = character(),
    "color" = character(),
    "fixationDepth" = numeric()
)
paths <- attr(fixed, "paths")
tree <- attr(paths, "tree")
edgeLengths <- node.depth.edgelength(tree)
for (site in allSitesName(paraSites)) {
    if (!site %in% allSitesName(fixed)) {
        next
    }
    aaPos <- as.integer(site)
    mutTips <- as.data.frame(sapply(
        X = extractTips(paraSites[[site]]),
        FUN = function(mutTips) {
            attr(mutTips, "mutName")[4]
        }
    ))
    
    colnames(mutTips) <- "mutName"
    
    mutDist <- sapply(split(rownames(mutTips), mutTips$mutName), function(tips) {
        min(sapply(tips, function(tip) {
            t <- which(tree$tip.label == tip)
            if (length(t)) {
                a <- tree$edge[which(tree$edge[, 2] == t), 1]
            } else {
                a <- as.integer(tip)
            }
            edgeLengths[a]
        }))
    })
    for (mut in names(mutDist)) {
        prevAA <- substr(mut, 1, 1)
        AA <- substr(mut, nchar(mut), nchar(mut))
        fullNameAA <- sitePath:::AA_FULL_NAMES[[tolower(AA)]]
        newRow <- data.frame(
            # Map the site to the reference
            "prevAA" = prevAA,
            "refSite" = aaPos,
            "fixedAA" = AA,
            "color" = sitePath:::AA_COLORS[[fullNameAA]],
            "fixationDepth" = mutDist[[mut]]
        )
        df <- rbind(df, newRow)
    }
}

df <- df[which(!duplicated(df[, c("prevAA", "refSite", "fixedAA")])), ]


# Prepare for the lolliplot
seqname <- paste("SARS-CoV-2", protein, sep = "_")
gr <- GRanges(seqname, IRanges(df[["refSite"]],
                               width = 1,
                               names = paste0(df[["prevAA"]], df[["refSite"]], df[["fixedAA"]])))
gr$color <- as.character(df[["color"]])
gr$score <- (df[["fixationDepth"]] - min(df[["fixationDepth"]])) * 5000
gr$label <- as.character(df[["fixedAA"]])
gr$label.col <- "#DBD9DA"

# Prepare the domains
domainInfo <- domain4Lolli(domainSites, maxLen)
features <- GRanges(seqname, IRanges(
    "start" = domainInfo$blockStarts,
    "width" = domainInfo$blockWidths,
    "names" = domainInfo$blockNames
))
features$fill <- DOMAIN_COLOR[domainInfo$blockNames]
features$color <- 0
features$height <- domainInfo$blockHeights

pdf(
    file = file.path(PLOTS_DIR, paste0(protein, "_lolliplot_para.pdf")),
    width = 16,
    height = 4,
    bg = "transparent"
)
lolliplot(
    SNP.gr = gr,
    features = features,
    ranges = GRanges(seqname, IRanges(0, maxLen)),
    ylab = "Distance to root",
    yaxis = FALSE
)
invisible(dev.off())

lolliplot(
    SNP.gr = gr,
    features = features,
    ranges = GRanges(seqname, IRanges(0, maxLen)),
    ylab = "Distance to root",
    yaxis = FALSE
)


p <- plot(fixed) + theme(legend.position = "none")
print(p)
ggsave(
    filename = "Output/sitePath_fixation.svg",
    plot = p,
    device = "svg",
    width = 10,
    height = 6
)

p <- plot(paraSites, y = F)
print(p)
ggsave(
    filename = "Output/sitePath_parallel.svg",
    plot = p,
    device = "svg",
    width = 10,
    height = 6
)
