#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(fs)
library(shiny)
library(sitePath)

options(list("shiny.host" = "0.0.0.0",
             "shiny.port" = 9000))

DATA_DIR <- "webdata"
SITEPATH_RES_DIR <- file.path(DATA_DIR, "sitePath_results")

availableDate <- list.files(SITEPATH_RES_DIR)

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("Nextstrain SARS-CoV-2 sites"),
    
    # Sidebar with a slider input for date
    sidebarLayout(sidebarPanel(
        selectInput(
            inputId = "date",
            label = "Date:",
            choices = availableDate,
            selected = as.character(max(as.Date(availableDate)))
        ),
        selectInput(
            inputId = "site",
            label = "Site:",
            choices = "all"
        )
    ),
    
    # Show a plot of the generated plot
    mainPanel(h4(textOutput(
        "date"
    )),
    plotOutput("distPlot")))
)

# Define server logic required to draw the plot
server <- function(input, output, session) {
    fixedSites <- reactive({
        # read RDS file from ui.R
        readRDS(file.path(SITEPATH_RES_DIR, input$date, "fixation.rds"))
    })
    output$distPlot <- renderCachedPlot({
        fixed <- fixedSites()
        # draw the plot based on the site
        site <- input$site
        if (site == "all" || !site %in% allSitesName(fixed)) {
            plot(fixed) + ggplot2::theme(legend.position = "none")
        } else {
            plotSingleSite(fixed, site)
        }
    }, cacheKeyExpr = { c(input$date, input$site) })
    observe({
        fixed <- fixedSites()
        sites <- allSitesName(fixed)
        selected <- "all"
        if (input$site %in% sites) {
            selected <- input$site
        }
        sites <- sort(as.integer(sites))
        updateSelectInput(
            session = session,
            inputId = "site",
            choices = c("all", sites),
            selected = selected
        )
    })
    output$date <- renderText(input$date)
}

# Run the application
shinyApp(ui = ui, server = server)
