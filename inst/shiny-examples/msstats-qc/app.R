library(shiny)
library(shinyBS)
library(shinyjs)
library(plotly)
library(RecordLinkage)
library(hash)
library(gridExtra)
library(ggExtra)
library(markdown)
library(shinythemes)

packageBaseDir <- system.file(package = "MSstatsQC")
appDir <- paste0(packageBaseDir,"/shiny-examples/msstats-qc")
source(paste0(appDir,"/plot-functions.R"))
source(paste0(appDir,"/data-validation.R"))
source(paste0(appDir,"/helper-functions.R"))
source(paste0(appDir,"/QCMetrics.R"))
runner <- function(custom_data) {
  shinyApp(
    ui = fluidPage(
      shinyjs::useShinyjs(),
      titlePanel(title=p(strong("MSstatsQC"),align = "center",style="color:#0A4476;",style="font-size:170%;",
                         style="font-family:arial;"),windowTitle = "MSstatsQC"),
      navbarPage(h4("Longitudinal system suitability monitoring and quality control for targeted proteomic experiments",style="color:darkblue;"),
                 ##888888

                 #################################################################################################################
                 tabPanel("Home",
                          tags$img(src=paste0(appDir,'/www/logo.png'), height=220, width=220, style = "float: right"),
                          tags$img(src=paste0(appDir,'/www/home.png'), height=200, width=500, style = "float: left"),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          p("MSstatsQC is an open-source web-based software which provides longitudinal
                            system suitability monitoring tools (control charts) for proteomic experiments."),
                          h5(strong("Metrics you can monitor")),
                          p("MSstatsQC uses control charts to monitor the instrument performance by tracking system
                            suitability metrics including total peak area, retention time and full width at half maximum (FWHM) and peak assymetry.
                            Additional metrics can also be analyzed by including them to the input file."),
                          h5(strong("Statistical functionalities")),
                          p("This framework includes simultaneous monitoring tools for mean and dispersion of suitability metrics and presents
                            alternative methods of monitoring such as time weighted control charts to ensure that various types
                            of process disturbances are detected effectively. Simultaneous control charts used in this framework
                            can be classified into two groups: individual-moving range (XmR) control charts and mean and dispersion
                            cumulative sum (CUSUM) control charts. To successfully identify the time of change, change point analysis
                            is also included in this framework. Experiment specific control limits are provided with the control
                            charts to distinguish between random noise and systematic error. MSstatsQC can also help user on decision making.
                            Decision regions (red, yellow and blue) can be designed with 'Create Decision Rules' tab and results are available in 'Metric Summary' tab."),
                          h5(strong("Using MSstatsQC")),
                          p("The steps for generating results are as follows:"),
                          p("INCLUDE CHEATSHEET WORKFLOW!!!"),
                          ("1) Import your SST/QC data "),
                          br(),
                          ("2)	Determine the guide set to estimate metric mean and variance "),
                          br(),
                          ("3)	Select specific precursor(s) or select all"),
                          br(),
                          ("4)	Design decision rules"),
                          br(),
                          ("5) Run and generate control charts"),
                          br(),
                          ("6)	Check with heatmaps, metric summary plots and change point analysis for better reasoning"),
                          br(),
                          ("7)	Navigate results and download them for your QC reports"),
                          br(),
                          br(),

                          br(),
                          br(),
                          h5 ("Project Team: "),
                          h5("Eralp Dogu,",span("eralp.dogu@gmail.com",style = "color:blue")),
                          h5("Sara Taheri,",span("mohammadtaheri.s@husky.neu.edu",style = "color:blue")),
                          h5("Olga Vitek,",span("o.vitek@neu.edu",style = "color:blue")),
                          br(),
                          br(),
                          ("Olga Vitek Lab"),
                          br(),
                          ("College of Science"),
                          br(),
                          ("College of Computer and Information Science"),
                          br(),
                          ("360 Huntington Ave"),
                          br(),
                          ("Boston, Massachusetts 02115"),
                          br(),
                          br(),
                          br()

                          ),
                 ########################################################################################################
                 tabPanel("Data Import and Selection",
                          tabsetPanel(
                            tabPanel("Data Import",
                                     sidebarLayout(

                                       sidebarPanel(
                                         wellPanel(
                                           p("Upload your data (Comma-separated (*.csv) QC file format)"),

                                           p("To see acceptable example data, look at", strong("Help"),"tab"),

                                           fileInput("filein", "Upload file")
                                         ),

                                         wellPanel(
                                           p("If you want to run", strong("MSstatsQC"), "with example data file, click this button"),
                                           actionButton("sample_button", "Run with example data")
                                           #bsTooltip("sample_button","If you want to run MSstatsQC with example data file, click this button", placement = "bottom", trigger = "hover",
                                           #options = NULL)
                                         ),

                                         wellPanel(
                                           p("If you want to clean existing results, click this button"),

                                           actionButton("clear_button", "Clear data and plots")

                                           #bsTooltip("clear_button","click this button to clear your data and all the tables and plots from the system.", placement = "bottom", trigger = "hover",
                                           #options = NULL)
                                         ),


                                         tags$style("body{background-color:linen; color:black}")


                                       ),
                                       mainPanel(

                                         tabPanel("Data",
                                                  dataTableOutput('prodata_table'))
                                       ),
                                       position = "left")
                            ),

                            tabPanel("Options",
                                     p(strong("Select QC metrics for all the subsequence analyses:")),

                                     wellPanel(
                                       fluidRow(
                                         column(10,
                                                uiOutput("metricSelection")
                                         )
                                       )
                                     ),
                                     wellPanel(
                                       radioButtons("selectGuideSetOrMeanSD",

                                                    "If you want to select mean and standard deviation yourself select them here. Otherwise choose the guide set button.",
                                                    #"define mean and standard deviation",
                                                    choices = c("Mean and standard deviation estimated by the user","Mean and standard deviation estimated from guide set")
                                       ),
                                       conditionalPanel(
                                         condition = "input.selectGuideSetOrMeanSD == 'Mean and standard deviation estimated by the user'",
                                         p("Select the mean and standard deviation"),
                                         uiOutput("selectMeanSD")
                                       ),
                                       conditionalPanel(
                                         condition = "input.selectGuideSetOrMeanSD == 'Mean and standard deviation estimated from guide set'",
                                         p("Select a guide set to estimate control limits"),

                                         uiOutput("selectGuideSet")
                                       )
                                     ),
                                     wellPanel(
                                       p("Select a precursor or select all"),
                                       uiOutput("pepSelect")
                                     )
                            )
                          )
                 ),
                 ######################################################################################################
                 tabPanel("Create Decision Rules", theme = "bootstrap.css",
                          fluidPage(
                            p(strong("Create your decision rule:")),
                            wellPanel(
                              fluidRow(
                                p(strong("RED FLAG"), style="color:black; background-color: red;",align = "center",style="font-size:125%;"),
                                p(strong("System performance is UNACCEPTABLE when:"),align = "center"),
                                p("1. greater than the selected % of peptides are", strong("out of control"),"and"),
                                p("2. greater than the selected # of QC metrics are", strong("out of control."))
                              ),
                              fluidRow(
                                column(2,
                                       br()
                                ),
                                column(5,
                                       p(strong("% out of control peptides: ")),
                                       numericInput('threshold_peptide_red', '', value = 70, min = 0, max = 100, step = 1)
                                ),
                                column(5,
                                       p(strong("# out of control QC metrics: ")),
                                       uiOutput("metricThresholdRed")
                                )
                              )
                            ),

                            wellPanel(
                              fluidRow(
                                p(strong("YELLOW FLAG"), style="color:black; background-color: yellow;",align = "center",style="font-size:125%;"),
                                p(strong("System performance is POOR when:"),align = "center"),
                                p("1. greater than the selected % of peptides are", strong("out of control"),"and"),
                                p("2. greater than the selected # of QC metrics are", strong("out of control.")),
                                p("Warning:The limits should be less than or equal to the the RED FLAG limits")
                              ),
                              fluidRow(
                                column(2,
                                       br()
                                ),
                                column(5,
                                       p(strong("% of out of control peptides: ")),
                                       uiOutput("peptideThresholdYellow")
                                ),
                                column(5,
                                       p(strong("# of out of control metrics: ")),
                                       uiOutput("metricThresholdYellow")
                                )
                              )
                            ),
                            wellPanel(
                              fluidRow(
                                p(strong("BLUE FLAG"), style="color:black; background-color: blue;",align = "center",style="font-size:125%;"),
                                p(strong("System performance is ACCEPTABLE when:"),align = "center"),
                                p("RED FLAG and YELLOW FLAG limits are not exceeded.")
                              )

                            )
                          )
                 ),
                 #####################################################################################################
                 tabPanel("Metric Summary",
                          tabsetPanel(

                            tabPanel("Descriptives: Boxplots for QC Metrics",
                                     tags$head(tags$style(type="text/css")),
                                     conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                      tags$div("It may take a while to load the plots, please wait...",
                                                               id="loadmessage")),
                                     plotlyOutput("box_plot", height = 2000)
                            ),

                            tabPanel("Overall Performance: Decision-maps",
                                     tags$head(tags$style(type="text/css")),
                                     conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                      tags$div("It may take a while to load the plots, please wait...",
                                                               id="loadmessage")),
                                     sidebarLayout(
                                       sidebarPanel(
                                         checkboxGroupInput("heatmap_controlChart_select", "Select your control chart",
                                                            choices = c("CUSUM Charts" = "CUSUM","XmR Chart" = "XmR"), selected = "XmR")
                                         #htmlOutput("heatmap_txt")
                                       ),
                                       mainPanel(plotOutput("heat_map")
                                       )
                                     )
                            ),

                            tabPanel("Detailed Performance: Plot summaries",
                                     tags$head(tags$style(type="text/css")),
                                     conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                      tags$div("It may take a while to load the plots, please wait...",
                                                               id="loadmessage")),
                                     sidebarLayout(
                                       sidebarPanel(
                                         checkboxGroupInput("summary_controlChart_select", "Select your control chart",
                                                            choices = c("CUSUM Charts" = "CUSUM","XmR Chart" = "XmR"), selected = "XmR")
                                         #htmlOutput("summary_decision_txt")
                                       ),
                                       mainPanel(
                                         plotOutput("plot_summary")
                                       )
                                     )
                            )
                          )
                 ),
                 ###################################################################################################
                 navbarMenu("Control Charts",
                            tabPanel("XmR Control Charts",
                                     uiOutput("XmR_tabset")
                            ),

                            tabPanel("CUSUMm and CUSUMv Control Charts",
                                     uiOutput("CUSUM_tabset")
                            ),
                            tabPanel("Change Point Analysis for Mean and Variability",
                                     uiOutput("CP_tabset")
                            )
                 ),
                 ###################################################################################################
                 tabPanel("Help",
                          tabsetPanel(
                            tabPanel("Metrics"
                                     ,h5(strong("Retention Time")),
                                     p("Retention time is the time it takes a solute to travel through the column. The retention time is assigned to
                                       the corresponding solute peak. The retention time is a measure of the amount of time a solute spends in a column.
                                       It is the sum of the time spent in the stationary phase and the mobile phase."
                                       ,a("visit for more info",href="http://www.britannica.com/science/retention-time")),
                                     br(),
                                     h5(strong("Total Peak Area")),
                                     p("Total Peak Area is the sum of all integrated signals for a certain peptide."),
                                     br(),
                                     h5(strong("Full Width at Half Maximum (FWHM)")),
                                     p("Full width at half maximum 'FWHM' is an expression of the extent of a
                                       function given by the difference between the two extreme values
                                       of the independent variable at which the dependent variable is equal
                                       to half of its maximum value."
                                       ,a("visit for more info",href="https://en.wikipedia.org/wiki/Full_width_at_half_maximum")),
                                     br(),
                                     h5(strong("Peak Assymetry")),
                                     p("Peak Assymetry is a measure of symetry for a peak. Calculated by taking 2*a/(a+b). Optimal value is around 1 for a Gaussian peak.")
                            ),

                            tabPanel("Plots"

                                     ,h5(strong("XmR control charts")),
                                     h5("Can detect large shifts and spikes in the mean and dispersion of suitability metric."),
                                     h5("The sequential differences between successive values as a measure
                                        of dispersion and individual observations are used to construct the plots."),
                                     p("A measure of dispersion can be estimated by computing the ranges of two consecutive observations. This
                                       approach is used to construct XmR charts. XmR charts plot
                                       original observations and moving ranges to investigate deviations from random process behaviour.
                                       XmR chart consists of 2 charts;
                                       Individuals (X) chart and Moving Range (mR) chart.
                                       ",
                                       a("visit for more info",href="https://en.wikipedia.org/wiki/Shewhart_individuals_control_chart")),

                                     br(),
                                     h5(strong("CUSUMm and CUSUMv control charts")),
                                     h5("Can detect small shifts and sustained drifts in the mean and dispersion of suitability metric."),
                                     h5("Time weighted cumulative sums are used to construct the plots."),
                                     p("A CUSUM chart is a time-weighted control chart that displays the cumulative sums
                                       'CUSUMs'. Because it is cumulative, even minor drifts in the process mean or dispersion and gradual deterioration of quality will cause steadily
                                       increasing or decreasing cumulative values. We introduce two CUSUM control charts: a mean CUSUM (CUSUMm) and
                                       a dispersion CUSUM (CUSUMv).",
                                       a("visit for more info",href="https://en.wikipedia.org/wiki/CUSUM")),

                                     br(),
                                     h5(strong("Change Point Analysis")),
                                     h5("Can identify the exact time of a change in the mean and dispersion of suitability metric. "),
                                     h5("Likelihood functions are plotted and the QC sample which maximizes the functions is considered as a candidate change point."),
                                     p("A change in the process parameters triggers a control chart to generate an out of
                                       control signal. The QC sample at which the signal is issued is considered as the
                                       stopping time and after the signal search for an assignable cause is recommended.
                                       However, the signal does not always designate that the special cause actually occurred
                                       at that certain time. A remedy to this problem is to use follow-up change
                                       point analysis along with control charts. Change point estimation procedures have a potential
                                       to save time by narrowing the search window. We introduce
                                       two change point models: step shift change model for mean and step shift change model for variance.",
                                   a("visit for more info",href="http://www.eng.fsu.edu/~pigna/pdf/"))
                                  ),

                         tabPanel("Documentation",
                                  h5(strong("MSstatsQC webpage")),
                                  p("Source codes, related documents and user manual can be found via our MSstats website"
                                    , a("visit for more info",href="http://www.msstats.org/msstatsqc")),
                                  br(),
                                  h5(strong("MSstatsQC Github")),
                                  p("Latest documantation is also available via our Github page"
                                    , a("visit for more info",href="https://github.com/srtaheri/msstats-qc"))
                                  )

                                  ))
#####################################################################################################
                       )
                            ),
    ###################################################################
    server = function(input, output,session) {
      #### Read data  ##################################################################################################
      data <- reactiveValues(df = NULL, metrics = NULL)

      observeEvent(custom_data, {
        data$df <- input_checking(custom_data)
        validate(
          need(!is.null(data$df), "Please upload your data"),
          need(is.data.frame(data$df), data$df)
        )
        data$metrics <- c(find_custom_metrics(data$df))
      })


      observeEvent(input$filein, {
        custom_data <- NULL
        file1 <- input$filein
        data$df <- input_checking(read.csv(file=file1$datapath, sep=",", header=TRUE, stringsAsFactors=TRUE))
        validate(
          need(!is.null(data$df), "Please upload your data"),
          need(is.data.frame(data$df), data$df)
        )
        #data$metrics <- c(COL.BEST.RET, COL.TOTAL.AREA, COL.FWHM, COL.PEAK.ASS, find_custom_metrics(data$df))
        data$metrics <- c(find_custom_metrics(data$df))
      }, priority = 20)

      observeEvent(input$sample_button, {
        custom_data <- NULL
        data$df <- input_checking(read.csv(paste0(appDir,"/Datasets/Sampledata_CPTAC_Study_9_1_Site54.csv")))
        validate(
          need(!is.null(data$df), "Please upload your data"),
          need(is.data.frame(data$df), data$df)
        )

        #data$metrics <- c(COL.BEST.RET, COL.TOTAL.AREA, COL.FWHM, COL.PEAK.ASS, find_custom_metrics(data$df))
        data$metrics <- c(find_custom_metrics(data$df))
      }, priority = 20)

      observeEvent(input$clear_button, {
        custom_data <- NULL
        data$df <- NULL
        data$metrics <- NULL
      }, priority = 20)
      ##### Precursor type selection #####################################################################################
      output$pepSelect <- renderUI({
        prodata <- data$df
        validate(
          need(!is.null(prodata), "Please upload your data"),
          need(is.data.frame(prodata), prodata)
        )
        selectInput("pepSelection","Choose peptide"
                    #,choices = c(levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET])),"all peptides")
                    ,choices = c(levels(prodata$Precursor),"all peptides")
        )
      })
      ######Show table of data #####################################################################################################
      output$prodata_table <- renderDataTable({
        validate(
          need(!is.null(data$df), "Please upload your data"),
          need(is.data.frame(data$df), data$df)
        )
        data$df
      }, options = list(pageLength = 25))
      ###### selection tab in Data Improt and selection #####################################################
      output$selectMeanSD <- renderUI({
        lapply(input$user_selected_metrics,
               function(x){
                 fluidRow(
                   column(4,paste(x,":")),
                   column(4,
                          numericInput(paste0("selectMean@",x),"mean",value = 1)
                   ),
                   column(4,
                          numericInput(paste0("selectSD@",x),"standard deviation",value = 1)
                   )
                 )
               })
      })


      output$selectGuideSet <- renderUI({
        fluidRow(
          column(6,
                 numericInput("L","Lower bound of guide set",value = 1, min = 1, step = 1)
          ),
          column(6,
                 numericInput("U","Upper bound of guide set", value = 5, min = 2, step = 1)
          )
        )


      })
      ###### Tab for selecting decision rule and metrics ###############################################
      output$metricThresholdRed <- renderUI({
        numOfMetrics <- length(input$user_selected_metrics)
        numericInput('threshold_metric_red', '', value = 2, min = 0, max = numOfMetrics, step = 1)
      })

      output$peptideThresholdYellow <- renderUI({
        threshold_peptide_red <- input$threshold_peptide_red
        numericInput('threshold_peptide_yellow', '', value = threshold_peptide_red - 1, min = 1, max = threshold_peptide_red, step = 1)
      })

      output$metricThresholdYellow <- renderUI({
        numOfMetrics <- length(input$user_selected_metrics)
        threshold_metric_red <- input$threshold_metric_red
        validate(
          need(!is.null(numOfMetrics),"loading..."),
          need(!is.null(threshold_metric_red),"loading...")
        )

        numericInput('threshold_metric_yellow', '', value = threshold_metric_red , min = 0, max = threshold_metric_red, step = 1)
      })

      output$metricSelection <- renderUI({
        checkboxGroupInput("user_selected_metrics","",
                           choices = c(data$metrics),
                           #selected = c(COL.PEAK.ASS,COL.BEST.RET,
                           #            COL.FWHM, COL.TOTAL.AREA),
                           inline = TRUE)
      })

      ################################################################# plots ###################################################
      #################################################################################################################
      output$XmR_tabset <- renderUI({

        validate(
          need(!is.null(data$df), "Please upload your data first"),
          need(is.data.frame(data$df), data$df),
          need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule")
        )
        is_guidset_selected <- FALSE
        if(input$selectGuideSetOrMeanSD == "Mean and standard deviation estimated from guide set") {
          is_guidset_selected <- TRUE
        }
        Tabs <- lapply(input$user_selected_metrics,
                       function(x) {
                         tabPanel(x,
                                  tags$head(tags$style(type="text/css")),
                                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                   tags$div("It may take a while to load the plots, please wait...",
                                                            id="loadmessage")),
                                  renderPlotly(render.QC.chart(data$df, input$pepSelection, input$L,
                                                               input$U, metric = x,
                                                               plot.method = "XmR", normalization = FALSE,
                                                               y.title1 = "Individual Value", y.title2 = "Moving Range",
                                                               selectMean = input[[paste0("selectMean@",x)]],selectSD = input[[paste0("selectSD@",x)]],
                                                               guidset_selected = is_guidset_selected)
                                  )

                         )
                       })
        do.call(tabsetPanel, Tabs)

      })
      ################################################################################################################
      #################################################################################################################
      output$CUSUM_tabset <- renderUI({
        validate(
          need(!is.null(data$df), "Please upload your data first"),
          need(is.data.frame(data$df), data$df),
          need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule")
        )
        is_guidset_selected <- FALSE
        if(input$selectGuideSetOrMeanSD == "Mean and standard deviation estimated from guide set") {
          is_guidset_selected <- TRUE
        }
        Tabs <- lapply(input$user_selected_metrics,
                       function(x) {

                         tabPanel(x,
                                  tags$head(tags$style(type="text/css")),
                                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                   tags$div("It may take a while to load the plots, please wait...",
                                                            id="loadmessage")),
                                  renderPlotly(render.QC.chart(data$df, input$pepSelection, input$L, input$U, metric = x, plot.method = "CUSUM", normalization = TRUE, y.title1 = "CUSUM mean", y.title2 = "CUSUM variation",selectMean = input[[paste0("selectMean@",x)]],selectSD = input[[paste0("selectSD@",x)]],guidset_selected = is_guidset_selected))
                         )
                       })

        do.call(tabsetPanel, Tabs)
      })
      ################################################################################################################
      #################################################################################################################
      output$CP_tabset <- renderUI({
        validate(
          need(!is.null(data$df), "Please upload your data first"),
          need(is.data.frame(data$df), data$df),
          need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule")
        )
        is_guidset_selected <- FALSE
        if(input$selectGuideSetOrMeanSD == "Mean and standard deviation estimated from guide set") {
          is_guidset_selected <- TRUE
        }
        Tabs <- lapply(input$user_selected_metrics,
                       function(x) {
                         tabPanel(x,
                                  tags$head(tags$style(type="text/css")),
                                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                   tags$div("It may take a while to load the plots, please wait...",
                                                            id="loadmessage")),
                                  renderPlotly(render.QC.chart(data$df, input$pepSelection, input$L, input$U, metric = x, plot.method = "CP", normalization = TRUE, y.title1 = "Change point for mean", y.title2 = "Change point for variation",selectMean = input[[paste0("selectMean@",x)]],selectSD = input[[paste0("selectSD@",x)]],guidset_selected = is_guidset_selected))
                         )
                       })

        do.call(tabsetPanel, Tabs)
      })
      ######################################################### height and width in Summary tab ########################################
      my_height <- reactive({
        if(length(input$user_selected_metrics) < 5) {
          my_height <- ceiling(length(input$summary_controlChart_select))*600
        }else if(length(input$user_selected_metrics) < 9) {
          my_height <- ceiling(length(input$summary_controlChart_select))*1200
        }else {
          my_height <- 1500
        }

      })
      # my_width <- reactive({
      #   if(length(input$user_selected_metrics) < 5) {
      #     my_height <- ceiling(ceiling(input$user_selected_metrics))*50
      #   }else if(length(input$user_selected_metrics) < 9) {
      #     my_height <- ceiling((input$user_selected_metrics)%%4)*50
      #   }else {
      #     my_height <- 1500
      #   }
      #
      # })
      heatmap_height <- reactive({
        heatmap_height <- ceiling(length(input$user_selected_metrics)*length(input$summary_controlChart_select))*230
      })

      heatmap_width <- reactive({
        prodata <- data$df
        heatmap_width <- nrow(prodata[prodata$Precursor == prodata$Precursor[1],])*20
        heatmap_width
      })
      ########################################################## box plot in Metric Summary tab ##########################################
      output$box_plot <- renderPlotly({
        prodata <- data$df
        validate(
          need(!is.null(prodata), "Please upload your data"),
          need(is.data.frame(prodata), prodata),
          need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule")
        )
        metrics_box.plot(prodata, data.metrics = input$user_selected_metrics)
      })

      ###############   summary plots and radar plots ############################################################################
      output$plot_summary <- renderPlot({

        prodata <- data$df
        validate(
          need(!is.null(prodata), "Please upload your data"),
          need(is.data.frame(prodata), prodata),
          need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule")
        )

        is_guidset_selected <- FALSE
        if(input$selectGuideSetOrMeanSD == "Mean and standard deviation estimated from guide set") {
          is_guidset_selected <- TRUE
        }
        listMean <- list()
        listSD <- list()
        for(metric in input$user_selected_metrics){
          listMean[[metric]] <- input[[paste0("selectMean@",metric)]]
          listSD[[metric]] <- input[[paste0("selectSD@",metric)]]
        }

        plots <- list()
        i <- 1
        for(method in input$summary_controlChart_select) {
          p1 <- NULL
          p2 <- NULL
          if(method == "XmR") {
            p1 <- XmR.Summary.plot(prodata, data.metrics = input$user_selected_metrics, input$L, input$U, listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected)
            p2 <- XmR.Radar.Plot(prodata, data.metrics = input$user_selected_metrics,input$L,input$U,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected)
          } else if(method == "CUSUM") {
            p1 <- CUSUM.Summary.plot(prodata, data.metrics = input$user_selected_metrics, input$L, input$U,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected)
            p2 <- CUSUM.Radar.Plot(prodata, data.metrics = input$user_selected_metrics, input$L,input$U,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected)
          }
          plots[[i]]   <- p1
          plots[[i+1]] <- p2

          i <- i+2
        }
        if(length(plots) > 0)
          do.call("grid.arrange", c(plots, ncol = 1))
      }, height = my_height)

      ################## decision message for XmR in summary tab #########
      #  output$summary_decision_txt <- renderUI({
      #    prodata <- data$df
      #    validate(
      #      need(!is.null(prodata), "Please upload your data"),
      #      need(is.data.frame(prodata), prodata),
      #      need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule")
      #    )
      #    peptideThresholdRed <- (as.numeric(input$threshold_peptide_red))/100 #this is the percentage of peptide user chooses for red flag
      #    metricThresholdRed <- as.numeric(input$threshold_metric_red) #this is the number of metric user chooses for red flag
      #    peptideThresholdYellow <- (as.numeric(input$threshold_peptide_yellow))/100 #this is the percentage of peptide user chooses for yellow flag
      #    metricThresholdYellow <- as.numeric(input$threshold_metric_yellow) #this is the number of metric user chooses for yellow flag
      #
      #    is_guidset_selected <- FALSE
      #    if(input$selectGuideSetOrMeanSD == "Mean and standard deviation estimated from guide set") {
      #      is_guidset_selected <- TRUE
      #    }
      #    listMean <- list()
      #    listSD <- list()
      #    for(metric in input$user_selected_metrics){
      #      listMean[[metric]] <- input[[paste0("selectMean@",metric)]]
      #      listSD[[metric]] <- input[[paste0("selectSD@",metric)]]
      #    }
      #    if("XmR" %in% input$summary_controlChart_select && "CUSUM" %in% input$summary_controlChart_select) {
      #      HTML(paste(
      #        decisionRule_warning_message_CUSUM(prodata,input$user_selected_metrics,method = "CUSUM", peptideThresholdRed,peptideThresholdYellow,metricThresholdRed,metricThresholdYellow,
      #                                           input$L, input$U, type = 2,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected),"\n\n",
      #        decisionRule_warning_message_XmR(prodata,input$user_selected_metrics,method = "XmR", peptideThresholdRed,peptideThresholdYellow,metricThresholdRed,metricThresholdYellow,
      #                                         input$L, input$U, type = 2,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected),
      #        sep = "<br/>"
      #      ))
      #    }else if("CUSUM" %in% input$summary_controlChart_select) {
      #      decisionRule_warning_message_CUSUM(prodata,input$user_selected_metrics,method = "CUSUM", peptideThresholdRed,peptideThresholdYellow,metricThresholdRed,metricThresholdYellow,
      #                                         input$L, input$U, type = 2,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected)
      #    }else if("XmR" %in% input$summary_controlChart_select) {
      #      decisionRule_warning_message_XmR(prodata,input$user_selected_metrics,method = "XmR", peptideThresholdRed,peptideThresholdYellow,metricThresholdRed,metricThresholdYellow,
      #                                       input$L, input$U, type = 2,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected)
      #    }else {
      #    }
      # })
      ################## decision message for XmR in heatmap tab #########
      # output$heatmap_txt <- renderUI({
      #   prodata <- data$df
      #   validate(
      #     need(!is.null(prodata), "Please upload your data"),
      #     need(is.data.frame(prodata), prodata),
      #     need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule")
      #   )
      #   peptideThresholdRed <- (as.numeric(input$threshold_peptide_red))/100 #this is the percentage of peptide user chooses for red flag
      #   metricThresholdRed <- as.numeric(input$threshold_metric_red) #this is the number of metric user chooses for red flag
      #   peptideThresholdYellow <- (as.numeric(input$threshold_peptide_yellow))/100 #this is the percentage of peptide user chooses for yellow flag
      #   metricThresholdYellow <- as.numeric(input$threshold_metric_yellow) #this is the number of metric user chooses for yellow flag
      #
      #   is_guidset_selected <- FALSE
      #   if(input$selectGuideSetOrMeanSD == "Mean and standard deviation estimated from guide set") {
      #     is_guidset_selected <- TRUE
      #   }
      #   listMean <- list()
      #   listSD <- list()
      #   for(metric in input$user_selected_metrics){
      #     listMean[[metric]] <- input[[paste0("selectMean@",metric)]]
      #     listSD[[metric]] <- input[[paste0("selectSD@",metric)]]
      #   }
      #   if("XmR" %in% input$heatmap_controlChart_select && "CUSUM" %in% input$heatmap_controlChart_select) {
      #     HTML(paste(
      #       decisionRule_warning_message_CUSUM(prodata,input$user_selected_metrics,method = "CUSUM", peptideThresholdRed,peptideThresholdYellow,metricThresholdRed,metricThresholdYellow,
      #                                        input$L, input$U, type = 2,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected),
      #       decisionRule_warning_message_XmR(prodata,input$user_selected_metrics,method = "XmR", peptideThresholdRed,peptideThresholdYellow,metricThresholdRed,metricThresholdYellow,
      #                                        input$L, input$U, type = 2,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected),
      #       sep = "<br/><br/>"
      #     ))
      #   }else if("CUSUM" %in% input$heatmap_controlChart_select) {
      #     decisionRule_warning_message_CUSUM(prodata,input$user_selected_metrics,method = "CUSUM", peptideThresholdRed,peptideThresholdYellow,metricThresholdRed,metricThresholdYellow,
      #                                      input$L, input$U, type = 2,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected)
      #   }else if("XmR" %in% input$heatmap_controlChart_select) {
      #     decisionRule_warning_message_XmR(prodata,input$user_selected_metrics,method = "XmR", peptideThresholdRed,peptideThresholdYellow,metricThresholdRed,metricThresholdYellow,
      #                                      input$L, input$U, type = 2,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected)
      #   }else {
      #   }
      #
      # })
      ############################# heat_map in Summary tab #############################################
      output$heat_map <- renderPlot({
        prodata <- data$df
        validate(
          need(!is.null(prodata), "Please upload your data"),
          need(is.data.frame(prodata), prodata),
          need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule"),
          need(!is.null(prodata$AcquiredTime),"To view heatmaps, the dataset should include Acquired Time column.")
        )

        peptideThresholdRed <- (as.numeric(input$threshold_peptide_red))/100
        peptideThresholdYellow <- (as.numeric(input$threshold_peptide_yellow))/100
        if(is.null(prodata$AcquiredTime)) return(NULL)

        is_guidset_selected <- FALSE
        if(input$selectGuideSetOrMeanSD == "Mean and standard deviation estimated from guide set") {
          is_guidset_selected <- TRUE
        }

        listMean <- list()
        listSD <- list()
        for(metric in input$user_selected_metrics){
          listMean[[metric]] <- input[[paste0("selectMean@",metric)]]
          listSD[[metric]] <- input[[paste0("selectSD@",metric)]]
        }

        plots <- list()
        i <- 1
        for(method in input$heatmap_controlChart_select) {
          p1 <- metrics_heat.map(prodata,
                                 data.metrics = input$user_selected_metrics, method = method,
                                 peptideThresholdRed, peptideThresholdYellow,input$L, input$U, type = 1,
                                 title = "Decision-map (Changes in mean of QC metric)",
                                 listMean = listMean, listSD = listSD, guidset_selected = is_guidset_selected)
          p2 <- metrics_heat.map(prodata,
                                 data.metrics = input$user_selected_metrics, method = method,
                                 peptideThresholdRed, peptideThresholdYellow,input$L, input$U, type = 2,
                                 title = "Decision-map (Changes in variability of QC metric)",
                                 listMean = listMean, listSD = listSD, guidset_selected = is_guidset_selected)
          plots[[i]]   <- p1
          plots[[i+1]] <- p2

          i <- i+2
        }
        if(length(plots) > 0)
          do.call("grid.arrange", c(plots, ncol = 1))

      }, height = heatmap_height, width = heatmap_width)
    }
)
}

