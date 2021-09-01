library(fs)
library(shiny)
library(sitePath)

options(list("shiny.host" = "0.0.0.0",
             "shiny.port" = 10000))

DATA_DIR <- "webdata"
SITEPATH_RES_DIR <- file.path(DATA_DIR, "sitePath_results")
AVAILABLE_DATES <- list.files(SITEPATH_RES_DIR)

ui <- fluidPage(
    # Application title
    titlePanel("Nextstrain SARS-CoV-2 sites"),
    # Sidebar with a slider input for date
    sidebarLayout(
        sidebarPanel(
            selectInput(
                inputId = "date",
                label = "Date:",
                choices = AVAILABLE_DATES,
                selected = as.character(max(as.Date(AVAILABLE_DATES)))
            ),
            selectInput(
                inputId = "site",
                label = "Site:",
                choices = "all"
            )
        ),
        # Show a plot of the generated plot
        mainPanel(
            h4(textOutput("date")),
            h5("Fixed sites"),
            plotOutput("fixedSitesPlot"),
            h5("Parallel sites"),
            plotOutput("paraSitePlot")
        ),
    )
)

# Define server logic required to draw the plot
server <- function(input, output, session) {
    fixedSites <- reactive({
        readRDS(file.path(SITEPATH_RES_DIR, input$date, "fixation.rds"))
    })
    paraSites <- reactive({
        readRDS(file.path(SITEPATH_RES_DIR, input$date, "parallel.rds"))
    })
    output$fixedSitesPlot <- renderCachedPlot({
        fixed <- fixedSites()
        # draw the plot based on the site
        site <- input$site
        if (site == "all" || !site %in% allSitesName(fixed)) {
            plot(fixed) + ggplot2::theme(legend.position = "none")
        } else {
            plotSingleSite(fixed, site)
        }
    }, cacheKeyExpr = {
        c(input$date, input$site)
    })
    output$paraSitePlot <- renderCachedPlot({
        para <- paraSites()
        # draw the plot based on the site
        plot(para)
    }, cacheKeyExpr = {
        c(input$date, input$site)
    })
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
