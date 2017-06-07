
library(shiny)
library(pmartRmems)
data("rRNA_obj")
data("rRNA_obj2")
mice <- rRNA_obj
dirt <- rRNA_obj2

###############################################################################

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectInput("choose_data", "Select Data", choices = c("mice", "dirt")),
      selectInput("choose_var", "Visualize by:", choices = c("Sample", "OTU")),
      sliderInput("shift", "Shift", min=0.01, max=2, value=0.5, step=0.01),
      checkboxInput("missing", "Include Missing Data in Plots", value = TRUE)
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
server <- function(input, output) {

  getdata <- reactive({
    if(input$choose_data == "mice") {
      ret <- mice
    } else {
      ret <- dirt
    }

    return(ret)
  })


  shrinkdata <- reactive({
    mydata <- getdata()

    if(input$choose_var == "Sample") {
      ret <- mydata$e_data[,2:7]
    } else {
      ret <- mydata$e_data[1:6,-1]
    }

    return(ret)
  })


  output$rawplot <- renderPlot({
    mydata <- shrinkdata()

    par(mfrow=c(2,3))
    i <- ifelse(input$choose_var == "Sample", 2, 1)
    if(input$missing) {
      x <- apply(mydata, i, hist)
    } else {
      x <- apply(mydata, i, function(x) hist(x[x>0]))
    }
  })


  output$mcplot <- renderPlot({
    mydata <- getdata()

    tempdata <- mydata$e_data[,-1]
    # ind_notzero <- lapply(tempdata, function(x) which(x!=0))
    ind_notzero <- tempdata != 0
    dirdata <- apply(tempdata, 2, function(x) log2(gtools::rdirichlet(1, x + input$shift)))
    temp <- apply(dirdata, 2, function(x){x - mean(x, na.rm = TRUE)})

    if(input$choose_var == "Sample") {
      ret <- temp[,1:6]
      i <- 2
      ind_notzero <- ind_notzero[,1:6]

      par(mfrow=c(2,3))
      if(input$missing) {
        x <- apply(ret, i, hist)
      } else {
        x <- sapply(1:6, function(i) {hist(ret[ind_notzero[,i],i])})
      }

    } else {
      ret <- temp[1:6,]
      i <- 1
      # ind_notzero <- lapply(ind_notzero, function(x) x[x<=6])
      ind_notzero <- ind_notzero[1:6,]

      par(mfrow=c(2,3))
      if(input$missing) {
        x <- apply(ret, i, hist)
      } else {
        # print(sapply(1:length(ind_notzero), function(i) {ret[ind_notzero[[i]],i]}))
        # print(length(ind_notzero))
        x <- sapply(1:6, function(i) {hist(ret[i,ind_notzero[i,]])})
      }
    }

    # par(mfrow=c(2,3))
    # if(input$missing) {
    #   x <- apply(ret, i, hist)
    # } else {
    #
    #   x <- sapply(1:length(ind_notzero), function(i) {hist(ret[ind_notzero[[i]],i])})
    # }
  })


  output$nomcplot <- renderPlot({
    mydata <- getdata()

    log2data <- log2(mydata$e_data[,-1] + input$shift)
    temp <- apply(log2data, 2, function(x){x - mean(x, na.rm = TRUE)})

    if(input$choose_var == "Sample") {
      ret <- temp[,1:6]
      i <- 2
    } else {
      ret <- temp[1:6,]
      i <- 1
    }

    par(mfrow=c(2,3))
    if(input$missing) {
      x <- apply(ret, i, hist)
    } else {
      x <- apply(ret, i, function(x) hist(x[x>min(x)]))
    }
  })



}

# Run the application
shinyApp(ui = ui, server = server)

