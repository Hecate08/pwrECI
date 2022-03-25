#########################
# reading in libaries   #
#########################
# libraries for R shiny
library(shiny)
#library(DT)
library(shinythemes)
library(shinyBS)
#library(shinyWidgets)
library(shinyhelper)
library(magrittr)
library(ipc)
library(future)
library(promises)
plan(multiprocess)


# libraries for ECI power function
library(truncnorm)
library(matrixStats)
library(limma)
library(ECEA)
source("ECI_power.R", local = TRUE)
source("ECI_functions_power.R", local = TRUE)

# libraries for the plots and the tables
library(ggplot2)
library(ggalt)


##################
# functions
##################
plotPower <- function(df){
  df2 <- data.frame(x = rep(df$x,2), Power = c(df$Power, df$adjPower),
                    type = rep(c("Power", "adjusted Power"), each = dim(df)[1]))

  p <- ggplot(df2, aes(x,Power, color = type)) +
    geom_point(color = "black") +
    geom_line(size = 0.75) +
    geom_hline(yintercept=0.8, linetype="dashed", color = "darkgrey") +
    xlab("Experiment") +
    theme_bw() + scale_x_continuous(breaks=df2$x,labels=as.factor(df2$x)) +
    scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0, 1)) +
    theme(panel.grid.minor = element_blank())

  return(p)
}



##################
# R shiny part  #
#################

ui = navbarPage(
  title = "pwrECI", # title for app
  theme = shinytheme("cerulean"), # sets color theme
  tabPanel(
    "Main",
    sidebarLayout(
      sidebarPanel(
        bsCollapsePanel(
          "Study 1",
          numericInput("conSD1", "sd of control group", value = 1, min = 0),
          numericInput("caseSD1", "sd of case group", value = 1, min = 0),
          numericInput("conMean1", "mean of control group", value = 2.5),
          numericInput("MeanDiff1", "difference of mean between groups", value = 0.5)
        ) %>% helper(icon = "question", content = "Set characteristics for study 1", type = "inline"),
        bsCollapsePanel(
          "Study 2",
          numericInput("conSD2", "sd of control group", value = 1, min = 0),
          numericInput("caseSD2", "sd of case group", value = 1, min = 0),
          numericInput("conMean2", "mean of control group", value = 3),
          numericInput("MeanDiff2", "difference of mean between groups", value = 1)
        ) %>% helper(icon = "question", content = "Set characteristics for study 2", type = "inline"),

        checkboxInput("unbal", "unbalanced sample size", FALSE),
        uiOutput("listORrangechoice1"),
        uiOutput("listORrange"),
        conditionalPanel(
          condition = "input.unbal",
          uiOutput("listORrangechoice2"),
          uiOutput("listORrange2")
          #textOutput("text")
        ),

        bsCollapsePanel(
          "FDR adjusted Power calculation",
          numericInput("alphaU", "alpha", value = 0.05, min = 0, max = 1),
          numericInput("m1", "number of prognostic genes", value = 300, min = 0),
          numericInput("m", "total number of genes m", value = 20000, min = 0),
          numericInput("f", "FDR level", value = 0.05, min = 0, max = 1)
        ) %>% helper(icon = "question", content = "variables regarding FDR adjusted power calculation", type = "inline"),

        bsCollapsePanel(
          "Other",
          numericInput("seed", "starting seed", value = 654654, min = 0),
          selectInput("ngenes", "number of iterations for power calculation",
                      c("500" = "500", "1000" = "1000", "10,000" = "10000"), selected = "1000")
        ) %>% helper(icon = "question", content = "other variables to modify", type = "inline"),
        actionButton("go", "Go")
        #,
        #actionButton('cancel', 'Cancel')

      ),
      mainPanel(
        textOutput("gS"),
        plotOutput("plot"),
        tableOutput("table2")

      )
    )
  ),
  tabPanel(
    "Instructions",
    includeMarkdown("instructions.md")

  ),
  tabPanel(
    "About",
    includeMarkdown("about.md")

  )
)


# modifying the output
server = function(input, output) {
  observe_helpers(withMathJax = TRUE)


  output$listORrangechoice1 <- renderUI({
    if(input$unbal){
      radioButtons("radio", label = "Group sizes for study1",
                   choices = list("Range" = "range", "List" = "list"), selected = "range")
    } else {
      radioButtons("radio", label = "Group sizes",
                   choices = list("Range" = "range", "List" = "list"), selected = "range")
    }
      })
  output$listORrange <- renderUI({
    if(is.null(input$radio)) return(NULL)
    if(input$radio == "list"){
      textInput("myList", "list input", "10 20 30")
    } else {
      splitLayout(
        numericInput("from1", "from", value = 10, min = 1),
        numericInput("to1", "to", value = 50, min = 1),
        numericInput("by1", "by", value = 20, min = 1)
      )
    }
  })
  output$listORrangechoice2 <- renderUI({
    radioButtons("radio2", label = "Group sizes for study 2",
                 choices = list("Range" = "range2", "List" = "list2"), selected = "range2")
  })
  output$listORrange2 <- renderUI({
    if(is.null(input$radio2)) return(NULL)
    if(input$radio2 == "list2"){
      textInput("myList2", NULL, "10 20 30")
    } else if(input$radio2 == "range2"){
      splitLayout(
        numericInput("from2", "from", value = 20, min = 1),
        numericInput("to2", "to", value = 60, min = 1),
        numericInput("by2", "by", value = 20, min = 1)
      )
    } else {}
  })

  gs1 <- reactive({
    if(is.null(input$radio)) return(NULL)
    gS <- c()
    if(input$radio == "list"){
      gs <- as.numeric(unlist(strsplit(as.character(input$myList), " ")))
    } else {
      if(is.null(input$from1) | is.null(input$to1) | is.null(input$by1)) return(NULL)
      gs <- seq(input$from1,input$to1, input$by1)
    }
  })
  gs2 <- reactive({
    if(is.null(input$radio2)) return(NULL)
    gS <- c()
    if(input$radio2 == "list2"){
      gs <- as.numeric(unlist(strsplit(as.character(input$myList2), " ")))
    } else {
      if(is.null(input$from2) | is.null(input$to2) | is.null(input$by2)) return(NULL)
      gs <- seq(input$from2,input$to2, input$by2)
    }
  })
  groupSize <- reactive({
    if(input$unbal){
      if(length(gs1()) != length(gs2())) return(NULL)
      groupS <- cbind(gs1(),gs2())
    } else {
      groupS <- gs1()
    }
    groupS
  })
  output$gS <- renderText({
    if(is.null(groupSize())) return("Please provide equal range for group sizes for both studies")
    else if(length(gs1()) > 50) return("Please test less than 50 experiments at once")
  })


  powerVals <- eventReactive(input$go,{
    if(is.null(groupSize())) return(NULL)
    else if(length(gs1()) > 50) return(NULL)

    progress <- shiny::Progress$new()
    progress$set(message = "Power calculation: ", value = 0)
    on.exit(progress$close())

    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        progress$set(detail = detail)
      } else {
        progress$inc(amount = value, detail = detail)
      }
    }


    data <- ECI_power(
      alphaU = input$alphaU,
      control_data_sd1 = input$conSD1,
      control_data_sd2 = input$conSD2,
      case_data_sd1 = input$caseSD1,
      case_data_sd2 = input$caseSD2,
      meanC1 = input$conMean1,
      meanC2 = input$conMean2,
      meanDiff1 = input$MeanDiff1,
      meanDiff2 = input$MeanDiff2,
      seed = input$seed,
      ngenes = as.numeric(input$ngenes),
      sizeG = groupSize(),
      m = input$m,
      m1 = input$m1,
      f = input$f,
      updateProgress = updateProgress,
      unbalanced = input$unbal)

    if(input$unbal){
      df <- data.frame(x = seq(1:length(groupSize()[,1])), gs = groupSize()[,1], gs2 = groupSize()[,2], Power = data[,1], adjPower = data[,2])
    } else {
      df <- data.frame(x = seq(1:length(groupSize())), gs = groupSize(), gs2 = groupSize(), Power = data[,1], adjPower = data[,2])
    }
    df
  })


  ###################
  # output
  # plot output
  output$plot <- renderPlot({
    if(is.null(powerVals())) return(NULL)
    plotPower(powerVals())
  })

  # table output
  output$table2 <- renderTable({
    if(is.null(powerVals())) return(NULL)


    table <- powerVals()
    colnames(table) <- c("Experiment","Group size study 1", "Group size study 2", "Power", "FDR adjusted Power")
    table[,2] <- as.integer(table[,2])
    table[,3] <- as.integer(table[,3])

    table
  }, digits = 3)

  # # other tabs
  # output$text <- renderText({
  #   "test"
  # })




}


shinyApp(ui = ui, server = server)
