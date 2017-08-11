# devtools::install_github("johansent/pmartRmems")

###############################################################################
# GLOBAL #
###############################################################################

library(shiny)
# library(pmartRqc)
library(pmartRmems)
library(ggplot2)
library(dplyr)

data("soil")
source('sample_library_sizes.R')

soil2 <- group_designation(soil, c("Time", "Treatment"))
oldgdf <- attributes(soil2)$group_DF
gdf <- oldgdf %>% filter(Treatment != "aerobic", Time != "t9d")
samps <- as.character(gdf$RNA)
(attributes(soil2)$group_DF <- gdf)
soil2$e_data <- soil$e_data %>% select(OTU, one_of(samps))

groups <- list(2*1:5, 2*1:5 + 1)

###############################################################################
# UI #
###############################################################################

ui <- fluidPage(
   shinyjs::useShinyjs(),
   titlePanel("Simulate Soil Data"),

   sidebarLayout(
      sidebarPanel(
        numericInput("effect_size", "Effect Size", value = 10),
        numericInput("nTP", "Number of True Positives in Simdata", value = 500),
        radioButtons("scaling_function", "Scaling Function", choices = c("i^2", "log(i+1)")),
        # checkboxInput("shift_zeros", "Shift zeros (by effect_size/2)"),
        checkboxInput("rerun_sim", "Simulate Base Data Each Time"),
        radioButtons("zero_method", "Method of dealing with zeros",
                     choices = c("Basic" = "shift_basic",
                                 "Shift zeros by effect_size/2" = "shift_zeros",
                                 "Add signal after transformation" = "shift_transform"),
                     selected = "shift_basic"),
        actionButton("go", "Simulate Data")
      ),


      mainPanel(
        # plotOutput("my_ggplot", brush = "brush"),
        column(6,
          plotOutput("my_ggplot", height = 600, click = "click")),
        column(6,
        tableOutput("brushtable"),
        verbatimTextOutput("OTU_info"))
      )
   )
)

###############################################################################
# SERVER #
###############################################################################

server <- function(input, output, session) {
  # reactive data
  sim <- reactiveValues(list = NULL)

  # update simdata
  observeEvent(input$go, {
    shinyjs::html("go", "Loading...")

    isolate({
      funtext <- input$scaling_function
      rerun_sim <- input$rerun_sim
      nTP <- input$nTP
      effect_size <- input$effect_size
      zero_method <- input$zero_method
    })

    shift <- ifelse(zero_method == "shift_zeros", TRUE, FALSE)
    transfun <- ifelse(funtext == "i^2", function(i) {i^2}, function(i) {log(i + 1)})
    base_data <- if (rerun_sim == TRUE) {NULL} else {sim$list$sim_data_base}

    sim$list <- simulate_count_data(soil2$e_data[,-1], treatment_cols = 2*1:5,
                                   nTP = nTP,
                                   effect_size = effect_size,
                                   transfun = transfun,
                                   shift = shift,
                                   base_data = base_data)
  })

  edata <- reactiveValues(list = NULL)

  pdat <- reactive({
    simdata <- req(sim$list)

    ## base data
    soil_temp0 <- soil2
    soil_temp0$e_data[,-1] <- simdata$sim_data_base
    filt <- imd_filter(soil_temp0)
    soil_temp <- applyFilt(filt, soil_temp0)
    soil_temp <- transform_data(soil_temp)
    pvals_base <- ttest_rrna(soil_temp$e_data, groups = groups)

    zero_method <- isolate(input$zero_method)
    is_trans_method <- zero_method == "shift_transform"
    if (is_trans_method) {

      isolate({
        funtext <- input$scaling_function
        rerun_sim <- input$rerun_sim
        nTP <- input$nTP
        effect_size <- input$effect_size
      })
      transfun <- ifelse(funtext == "i^2", function(i) {i^2}, function(i) {log(i + 1)})
      base_data <- soil_temp$e_data[,-1]

      simdata <- simulate_count_data(base_data, treatment_cols = 2*1:5,
                                      nTP = nTP,
                                      effect_size = effect_size,
                                      transfun = transfun,
                                      shift = FALSE,
                                      base_data = base_data)
    }

    ## signal data
    soil30 <- soil2
    if (is_trans_method) {
      soil30$e_data <- as.data.frame(cbind(soil_temp$e_data[[1]], simdata$sim_data))
      names(soil30$e_data)[1] <- "OTU"
    } else {
      soil30$e_data[,-1] <- as.data.frame(simdata$sim_data)
    }
    # store vector of significant OTUs
    sig_otus <- soil30$e_data[[1]][simdata$significant_rows]
    if (!is_trans_method) {
      filt <- imd_filter(soil30)
      soil3 <- applyFilt(filt, soil30)
    } else {
      soil3 <- soil30
    }
    # throw out OTUs that aren't in the base data
    soil3$e_data <- soil3$e_data[soil3$e_data[[1]] %in% soil_temp$e_data[[1]],]
    # significant rows in filtered data
    sigvec <- soil3$e_data[[1]] %in% sig_otus

    if (!is_trans_method) {
      soil3 <- transform_data(soil3)
    }

    # abs_val <- ifelse(is_trans_method, TRUE, FALSE)
    pvals <- ttest_rrna(soil3$e_data, groups = groups)#, abs_val = abs_val)

    # store edata
    edata$list <- list(soil_temp0 = soil_temp0$e_data, soil_temp = soil_temp$e_data,
                       soil30 = soil30$e_data, soil3 = soil3$e_data)

    # combine data
    pdat <- data.frame(OTU = soil3$e_data[[1]], pvals = pvals, base_pvals = pvals_base,
                       sim_sig = sigvec)
    pdat <- mutate(pdat, obs_sig = pvals < 0.05, base_sig = base_pvals < 0.05)

    return(pdat)
  })

  output$my_ggplot <- renderPlot({
    plotdat <- pdat()
    req(input$go)
    p <- ggplot(plotdat) + geom_point(aes(base_pvals, pvals, col = sim_sig)) +
      scale_y_continuous(breaks = 0.1*0:10) +
      scale_x_continuous(breaks = 0.1*0:10) +
      labs(title = "P-values for each OTU,\nbefore and after adding signal",
           x = "Before Signal",
           y = "After Signal",
           col = "Signal Added") +
      theme_grey(base_size = 14) + theme(legend.position = "bottom")
    shinyjs::html("go", "Simulate Data")
    return(p)
  })

  point_clicked <- reactiveValues(pt = NULL)

  output$brushtable <- renderTable({
    res <- nearPoints(pdat(), input$click)
    if (nrow(res) == 0)
      return()
    point_clicked$pt <- res[1,]
    return(point_clicked$pt)
  })

  output$OTU_info <- renderPrint({
    otu_compare <- req(point_clicked$pt)
    edata_list <- req(edata$list)
    wl <- list(base_data = edata_list$soil_temp0 %>% filter(OTU == otu_compare[[1]]),
                base_data_trans = edata_list$soil_temp %>% filter(OTU == otu_compare[[1]]),
                signal_data = edata_list$soil30 %>% filter(OTU == otu_compare[[1]]),
                signal_data_trans = edata_list$soil3 %>% filter(OTU == otu_compare[[1]]))

    df_list <- lapply(wl, function(x) {
      vec <- as.numeric(x[1,-1])
      df <- data.frame(control = vec[2*1:5 - 1], treatment = vec[2*1:5]) %>%
        mutate(diff = treatment - control)
      return(df)
    })

    return(df_list)
  })

}

# Run the application
shinyApp(ui = ui, server = server)

