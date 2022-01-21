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
library(coxed)
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
  p <- ggplot(df, aes(x,Power)) +
    geom_point(color = "black") +
    geom_line(size = 0.75) +
    #geom_xspline(spline_shape = -0.5, size = 0.75) +
    geom_hline(yintercept=0.8, linetype="dashed", color = "darkgrey") +
    xlab("group size") + ylim(0,1) +
    theme_bw() +
    scale_x_continuous(breaks=df$x,
                     labels=df$gs)

  return(p)
}



##################
# R shiny part  #
#################

ui = navbarPage(
  title = "ECI power calculation", # title for app
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
          "Other",
          #numericRangeInput("equiRange", label = "Range for equivalence", value = c(1,2.5), min = 0),
          numericInput("seed", "starting seed", value = 654654, min = 0),
          numericInput("ngenes", "number of iterations for power calculation", value = 1000, min = 1),
          numericInput("alphaU", "alpha", value = 0.05, min = 0, max = 1),
          checkboxInput("filter", "Filter ECI values for differential expression?", TRUE)
        ) %>% helper(icon = "question", content = "other variables to modify", type = "inline"),
        actionButton("go", "Go"),
        actionButton('cancel', 'Cancel')

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
    textOutput("Developers")

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
    #validate(need(!is.null(input$radio), "please select an input type"))
    if(input$radio == "list"){
      textInput("myList", "list input", "10 20 30")
    } else {
      splitLayout(
        numericInput("from1", "from", value = 10, min = 1),
        numericInput("to1", "to", value = 50, min = 1),
        numericInput("by1", "by", value = 20, min = 1)
      )
      #textInput("myRange", "range input (from, to, by)", "10, 50, 20")
    }
  })
  output$listORrangechoice2 <- renderUI({
    radioButtons("radio2", label = "Group sizes for study 2",
                 choices = list("Range" = "range2", "List" = "list2"), selected = "range2")
  })
  output$listORrange2 <- renderUI({
    if(is.null(input$radio2)) return(NULL)
    #validate(need(!is.null(input$radio), "please select an input type"))
    if(input$radio2 == "list2"){
      textInput("myList2", NULL, "10 20 30")
    } else if(input$radio2 == "range2"){
      splitLayout(
        numericInput("from2", "from", value = 20, min = 1),
        numericInput("to2", "to", value = 60, min = 1),
        numericInput("by2", "by", value = 20, min = 1)
      )
      #textInput("myRange2", "study 2 range input (from, to, by)", "10, 50, 20")
    } else {}
  })

  gs1 <- reactive({
    if(is.null(input$radio)) return(NULL)
    gS <- c()
    if(input$radio == "list"){
      gs <- as.numeric(unlist(strsplit(as.character(input$myList), " ")))
    } else {
      if(is.null(input$from1) | is.null(input$to1) | is.null(input$by1)) return(NULL)
      #val <- as.numeric(unlist(strsplit(as.character(input$myRange), ", ")))
      #gs <- seq(val[1],val[2],val[3])
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
      #val <- as.numeric(unlist(strsplit(as.character(input$myRange2), ", ")))
      #gs <- seq(val[1],val[2],val[3])
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
  })
  # output$gS <- renderText({
  #   c(input$conSD1,
  #     input$conSD2,
  #   input$caseSD1,
  #   input$caseSD2,
  #   input$conMean1,
  #   input$conMean2,
  #   input$MeanDiff1,
  #   input$MeanDiff2,
  #   min(input$equiRange),
  #   max(input$equiRange),
  #   input$seed)
  # })

  ###################
  #receiving data
  # interruptor <- AsyncInterruptor$new()    # To signal STOP to the future
  # powerVals <- reactiveVal()
  # running <- reactiveVal(FALSE)
  #
  # observeEvent(input$go,{
  #
  #   if(is.null(groupSize())) return(NULL)
  #
  #   if(running()) return(NULL)
  #   running(TRUE)
  #
  #   # Create new progress bar
  #   progress <- AsyncProgress$new(message="Complex analysis")
  #
  #   powerVals(NULL)
  #
  #   progress$set(message = "Power calculation: ", value = 0)
  #   #on.exit(progress$close())
  #
  #   updateProgress <- function(value = NULL, detail = NULL) {
  #     if (is.null(value)) {
  #       progress$set(detail = detail)
  #     } else {
  #       progress$inc(amount = value, detail = detail)
  #     }
  #   }
  #
  #   fut <- future({
  #     data <- ECI_power(
  #       alphaU = input$alphaU,
  #       control_data_sd1 = input$conSD1,
  #       control_data_sd2 = input$conSD2,
  #       case_data_sd1 = input$caseSD1,
  #       case_data_sd2 = input$caseSD2,
  #       meanC1 = input$conMean1,
  #       meanC2 = input$conMean2,
  #       meanDiff1 = input$MeanDiff1,
  #       meanDiff2 = input$MeanDiff2,
  #       seed = input$seed,
  #       ngenes = input$ngenes,
  #       sizeG = groupSize(),
  #       filter = input$filter,
  #       updateProgress = updateProgress,
  #       unbalanced = input$unbal,
  #       progressMonitor = function(i) interruptor$execInterrupts())
  #
  #     if(input$unbal){
  #       df <- data.frame(gs = groupSize()[,1], gs = groupSize()[,2], Power = data, x = seq(1:length(groupSize()[,1])))
  #     } else {
  #       df <- data.frame(gs = groupSize(), gs2 = groupSize(), Power = data, x = seq(1:length(groupSize())))
  #     }
  #     df
  #   }) %...>% powerVals
  #
  #   # Show notification on error or user interrupt
  #   fut <- catch(fut,
  #                function(e){
  #                  powerVals(NULL)
  #                  print(e$message)
  #                  showNotification(e$message)
  #                })
  #   # When done with analysis, remove progress bar
  #   fut <- finally(fut, function(){
  #     progress$close()
  #     running(FALSE) # Declare done with run
  #   })
  #
  #   # Return something other than the future so we don't block the UI
  #   NULL
  #
  # })
  #
  # # Send interrupt signal to future
  # observeEvent(input$cancel,{
  #   if(running())
  #     interruptor$interrupt("User Interrupt")
  # })
  #

  powerVals <- eventReactive(input$go,{
    #withProgress(message = 'Power calculation, please wait', value = 0, {
    if(is.null(groupSize())) return(NULL)

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
      #lowRange = min(input$equiRange),
      #highRange = max(input$equiRange),
      seed = input$seed,
      ngenes = input$ngenes,
      sizeG = groupSize(),
      filter = input$filter,
      updateProgress = updateProgress,
      unbalanced = input$unbal)

    if(input$unbal){
      df <- data.frame(gs = groupSize()[,1], gs = groupSize()[,2], Power = data, x = seq(1:length(groupSize()[,1])))
    } else {
      df <- data.frame(gs = groupSize(), gs2 = groupSize(), Power = data, x = seq(1:length(groupSize())))
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


    table <- powerVals()[,1:3]
    colnames(table) <- c("Group size study 1", "Group size study 2","Power")
    table[,1] <- as.integer(table[,1])
    table[,2] <- as.integer(table[,2])

    table
  }, digits = 3)

  # other tabs
  output$text <- renderText({
    "test"
  })
  output$Developers <- renderText({
    "Lisa Neums"
  })



}


shinyApp(ui = ui, server = server)
