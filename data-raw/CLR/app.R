
library(shiny)
library(pmartRmems)
library(ALDEx2)
library(dplyr)
library(tidyr)
library(ggplot2)

data("rRNA_obj")
data("rRNA_obj2")
data("selex")
mice <- rRNA_obj
soil <- rRNA_obj2

# choose OTUs with missing values
soil$e_data <- rbind(slice(soil$e_data, 801:806),
                     slice(soil$e_data, c(1:800, 807:nrow(soil$e_data))))

# get selex in pmartR format and choose OTUs with missing values
temp <- selex %>% mutate(OTU=row.names(selex)) %>% select(OTU, 1:3, 8:10, everything())
temp <- rbind(slice(temp, 91:96),
              slice(temp, c(1:90, 97:nrow(temp))))
selex <- list()
selex$e_data <- temp

theme_set(theme_grey(base_size = 18))

# breaks
base_breaks <- function(n = 2){
  function(x) {
    r <- range(x)
    if (r[1] <= 0) {r[1] <- 1}
    return(axisTicks(r, log = TRUE, nint = n))
  }
}

# plotting function
plotfun <- function(dat, ID, raw = TRUE) {
  # facet labels
  ## Missing
  dat$Missing <- factor(dat$Missing, labels = c("NonMissing", "Missing"))
  ## sample size
  newdat <- dat %>% group_by_(ID, "Missing") %>% summarise(n=n()) %>%
    mutate(nshow = paste("n =", n))

  p <- ggplot(dat) + geom_density(aes(Reads, fill = Missing), alpha = 0.7) +
    facet_grid(formula(paste("Missing ~", ID)), scales = "free_y") +
    scale_fill_discrete(guide = FALSE) +
    labs(x = "Transformed Reads") +
    geom_text(data = newdat, aes(x=Inf, y=Inf, hjust = 1, vjust = 1, label = nshow), size = 5)
  if (raw) {
    p <- p + scale_x_continuous(trans = "log1p", breaks = base_breaks()) +
    labs(x = "log(Reads + 1)")
  }
  return(p)
}
###############################################################################

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectInput("choose_data", "Select Data", choices = c("mice", "soil", "selex")),
      selectInput("choose_var", "Visualize by", choices = c("Sample", "OTU")),
      sliderInput("shift", "Shift", min=0.05, max=2, value=0.5, step=0.05,
                  animate = animationOptions(interval = 1200)),
      fluidRow(
        column(4, numericInput("shift_min", "Min. shift", value = 0.05)),
        column(4, numericInput("shift_max", "Max. shift", value = 2)),
        column(4, numericInput("shift_step", "Step size", value = 0.05))
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Raw Data", plotOutput("rawplot", height = "600px")),
        tabPanel("CLR With MC", plotOutput("mcplot", height = "600px")),
        tabPanel("CLR Without MC", plotOutput("nomcplot", height = "600px"))
      )
    )
  )
)

###############################################################################
server <- function(input, output, session) {

  # update shift slider
  observe({
    req(input$shift_min, input$shift_max, input$shift_step)
    mn <- input$shift_min
    mx <- input$shift_max
    step <- input$shift_step
    updateSliderInput(session, "shift", min = mn, max = mx, step = step)
  })

  # retrieve dataset
  getdata <- reactive({
    datname <- input$choose_data
    ret <- get(datname)
    return(ret)
  })


  # raw data tab
  output$rawplot <- renderPlot({
    mydata <- getdata()
    tempdata <- mydata$e_data

    if (input$choose_var == "Sample") {
      dat <- tempdata[,1:7] %>%
        gather("Sample", "Reads", 2:7) %>%
        mutate(Missing = (Reads == 0))
      p <- plotfun(dat, "Sample")

    } else {
      dat <- tempdata[1:6,] %>%
        gather("Sample", "Reads", 2:ncol(tempdata)) %>%
        mutate(Missing = (Reads == 0)) %>%
        select(OTU = 1, everything())
      p <- plotfun(dat, "OTU")
    }

    return(p)

  })


  # MC tab
  output$mcplot <- renderPlot({
    mydata <- getdata()
    tempdata <- mydata$e_data
    ind_zero <- tempdata == 0

    tempdata[,-1] <- apply(tempdata[,-1], 2, function(x) log2(gtools::rdirichlet(1, x + input$shift)))
    tempdata[,-1] <- apply(tempdata[,-1], 2, function(x){x - mean(x, na.rm = TRUE)})

    if (input$choose_var == "Sample") {
      dat <- tempdata[,1:7] %>%
        gather("Sample", "Reads", 2:7)
      dat <- ind_zero[,1:7] %>%
        as.data.frame() %>%
        gather("key", "Missing", 2:7) %>%
        select(Missing) %>%
        cbind(dat, .)
      p <- plotfun(dat, "Sample", FALSE)

    } else {
      dat <- tempdata[1:6,] %>%
        gather("Sample", "Reads", 2:ncol(tempdata))
      dat <- ind_zero[1:6,] %>%
        as.data.frame() %>%
        gather("key", "Missing", 2:ncol(tempdata)) %>%
        select(Missing) %>%
        cbind(dat, .) %>%
        select(OTU = 1, everything())
      p <- plotfun(dat, "OTU", FALSE)
    }

    return(p)
  })


  # No MC tab
  output$nomcplot <- renderPlot({
    mydata <- getdata()
    tempdata <- mydata$e_data
    ind_zero <- tempdata == 0

    tempdata[,-1] <- log2(tempdata[,-1] + input$shift)
    tempdata[,-1] <- apply(tempdata[,-1], 2, function(x){x - mean(x, na.rm = TRUE)})

    if (input$choose_var == "Sample") {
      dat <- tempdata[,1:7] %>%
        gather("Sample", "Reads", 2:7)
      dat <- ind_zero[,1:7] %>%
        as.data.frame() %>%
        gather("key", "Missing", 2:7) %>%
        select(Missing) %>%
        cbind(dat, .)
      p <- plotfun(dat, "Sample", FALSE)

    } else {
      dat <- tempdata[1:6,] %>%
        gather("Sample", "Reads", 2:ncol(tempdata))
      dat <- ind_zero[1:6,] %>%
        as.data.frame() %>%
        gather("key", "Missing", 2:ncol(tempdata)) %>%
        select(Missing) %>%
        cbind(dat, .) %>%
        select(OTU = 1, everything())
      p <- plotfun(dat, "OTU", FALSE)
    }

    return(p)
  })



}

# Run the application
shinyApp(ui = ui, server = server)

