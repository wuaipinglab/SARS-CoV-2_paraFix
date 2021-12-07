library(jsonlite)
library(ggplot2)
library(lemon)

PROTEIN_NAMES <- c("Spike", "N")
SAMPLING_METHODS <- c("nextstrain", "sampled")

BACKGROUND_NUM_FILE <- "Data/background_num.json"

background <- read_json(BACKGROUND_NUM_FILE, simplifyVector = TRUE)

pdf(file = "Output/metrics.pdf", width = 8, height = 4)
for (protein in PROTEIN_NAMES) {
    for (s_method in SAMPLING_METHODS) {
        sites_prevalency_file <-
            paste0("Data/sitesPrevalence_", protein, ".json")
        mutation_num_file <-
            paste0("Data/mutation_num_", protein, ".json")
        prevalenceInfoFile <-
            paste0("Data/prevalenceInfo_", protein, ".csv")
        
        
        baseline_file <-
            paste0("Data/", s_method, "_trees_unused.csv")
        threshold_file <-
            paste0("Data/", s_method, "_trees_unused_0.1.csv")
        paraFixSites_file <-
            paste0("Data/", s_method, "_sitePath_results.csv")
        homoplasy_res_file <-
            paste0("Data/", s_method, "_homoplasyFinder.csv")
        hyphy_res_file <-
            paste0("Data/", s_method, "_hyphy_results.csv")
        conserved_sites_file <-
            paste0("Data/", s_method, "_conserved_sites.csv")
        dates_file <- paste0("Data/", s_method, "_dates.json")
        
        plot_title <- paste0(protein, " (", s_method, ")")
        metric_plot <-
            paste0("Output/", s_method, "_", protein, ".svg")
        
        print(plot_title)
        
        allDates <- read_json(dates_file, simplifyVector = TRUE)
        mutationsNum <-
            read_json(mutation_num_file, simplifyVector = TRUE)
        sitesPrevalence <-
            read_json(sites_prevalency_file, simplifyVector = TRUE)
        
        conserved <- read.csv(conserved_sites_file)
        conserved <-
            conserved[which(conserved$product == protein),]
        conserved <- split(conserved$aaPos, conserved$date)
        
        allPrevalentSites <- integer()
        nonPrevalentSites <- integer()
        for (prevalency in sitesPrevalence) {
            allPrevalentSites <- c(allPrevalentSites, prevalency[["prevalent"]])
            nonPrevalentSites <-
                c(nonPrevalentSites, prevalency[["non_prevalent"]])
        }
        allPrevalentSites <- unique(allPrevalentSites)
        nonPrevalentSites <- unique(nonPrevalentSites)
        
        
        baselineSites <- read.csv(baseline_file)
        baselineSites <-
            baselineSites[which(baselineSites$product == protein),]
        baselineSites <-
            split(baselineSites[, c("aaPos", "mutRatio")], baselineSites$date)
        
        thresholdSites <- read.csv(threshold_file)
        thresholdSites <-
            thresholdSites[which(thresholdSites$product == protein),]
        thresholdSites <-
            split(thresholdSites[, c("aaPos", "mutRatio")], thresholdSites$date)
        
        
        paraFixSites <- read.csv(paraFixSites_file)
        paraFixSites <-
            paraFixSites[which(paraFixSites$product == protein),]
        paraFixSites <-
            split(paraFixSites[, c("aaPos", "type")], paraFixSites$date)
        
        
        hyphySites <- read.csv(hyphy_res_file)
        hyphySites <-
            split(hyphySites[, c("site", "posterior")], hyphySites$date)
        
        
        homoplasySites <- read.csv(homoplasy_res_file)
        homoplasySites <-
            homoplasySites[which(homoplasySites$product == protein),]
        homoplasySites <-
            split(homoplasySites[, c("aaPos", "consistencyIndex")], homoplasySites$date)
        
        
        getMetrics <- function(pred,
                               true_positive,
                               true_negative,
                               test_date,
                               method_name) {
            tp <- length(intersect(pred, true_positive))
            fp <- length(setdiff(pred, true_positive))
            
            tn <- length(setdiff(true_negative, pred))
            fn <- length(setdiff(true_positive, pred))
            return(
                data.frame(
                    "metric" = c("Specificity", "Sensitivity"),
                    "value" = c(tn / (tn + fp), tp / (tp + fn)),
                    "method" = method_name,
                    "date" = as.Date(test_date)
                )
            )
        }
        
        
        allMetrics <- do.call(rbind, lapply(allDates, function(d) {
            true_positive <- setdiff(allPrevalentSites, conserved[[d]])
            true_negative <-
                setdiff(nonPrevalentSites, conserved[[d]])
            # bs_res <- getMetrics(baselineSites[[d]][["aaPos"]],
            #                      true_positive,
            #                      true_negative,
            #                      d,
            #                      "MSA")
            th_res <- getMetrics(thresholdSites[[d]][["aaPos"]],
                                 true_positive,
                                 true_negative,
                                 d,
                                 "MSA_only")
            sp_res <- getMetrics(paraFixSites[[d]][["aaPos"]],
                                 true_positive,
                                 true_negative,
                                 d,
                                 "sitePath")
            hp_res <- getMetrics(hyphySites[[d]][["site"]],
                                 true_positive,
                                 true_negative,
                                 d,
                                 "hyphy_fubar")
            pa_res <-
                getMetrics(unique(homoplasySites[[d]][["aaPos"]]),
                           true_positive,
                           true_negative,
                           d,
                           "homoplasy")
            do.call(rbind, list(th_res, sp_res, hp_res, pa_res))
        }))
        
        # lineType <-
        #     factor(as.integer(allMetrics$method != "sitePath"))
        
        p <- ggplot(allMetrics,
                    aes(
                        x = date,
                        y = value,
                        group = method,
                        color = method
                    )) +
            geom_line() +
            geom_point(data = subset(allMetrics, method == "sitePath")) +
            facet_rep_grid( ~ metric) +
            theme_bw() +
            theme(
                strip.background = element_rect(colour = "white",
                                                fill = "white"),
                strip.text = element_text(size = 12, face = "bold"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(),
                # axis.line = element_blank(),
                panel.border = element_blank(),
                axis.text.x = element_text(
                    size = 10,
                    face = "bold",
                    angle = 60,
                    vjust = 0.75,
                    hjust = 0.75
                ),
                axis.text.y = element_text(size = 10,
                                           face = "bold"),
                axis.title = element_blank()
            ) +
            coord_capped_cart(bottom = 'both', left = 'both') +
            ylim(0, 1) +
            scale_color_manual(
                values = c(
                    "sitePath" = "purple",
                    "hyphy_fubar" = "green",
                    "homoplasy" = "orange",
                    "MSA_only" = "grey"
                )
            ) +
            scale_x_date(date_labels = "%b %y", date_breaks  = "1 month") +
            # scale_color_viridis(discrete = TRUE,) +
            ggtitle(protein)
        
        print(p)
        
        ggsave(
            filename = metric_plot,
            plot = p,
            device = "svg",
            width = 8,
            height = 4
        )
    }
}
dev.off()
