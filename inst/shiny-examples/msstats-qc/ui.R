library(shiny)
library(shinyBS)
library(plotly)
library(markdown)
library(shinythemes)
shinyUI(fluidPage(
  #theme = shinytheme("yeti"),
  shinyjs::useShinyjs(),
  titlePanel(title=p(strong("MSstatsQC"),align = "center",style="color:#0A4476;",style="font-size:170%;",
               style="font-family:arial;"),windowTitle = "MSstatsQC"),
  navbarPage(h4("Longitudinal system suitability monitoring and quality control for targeted proteomic experiments",style="color:darkblue;"),
             ##888888

#################################################################################################################
              tabPanel("Home",
                         tags$img(src='logo.png', height=220, width=220, style = "float: right"),
                         tags$img(src='home.png', height=200, width=500, style = "float: left"),
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
                         ("7)	Navigate results and download them for your reports"),
                         br(),
                         br(),

                         br(),
                         br(),
                         h5 ("Project team: "),
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
              tabPanel("Data import and selection",
                       tabsetPanel(
                         tabPanel("Data import",
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
                                  p(strong("Select metrics for all further analyses:")),

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
             tabPanel("Create decision rules", theme = "bootstrap.css",
                      fluidPage(
                        p(strong("Create your decision rule:")),
                        wellPanel(
                          fluidRow(
                            p(strong("RED FLAG"), style="color:black; background-color: red;",align = "center",style="font-size:125%;"),
                            p(strong("System performance is UNACCEPTABLE when:"),align = "center"),
                            p("1. greater than the selected % of peptides are", strong("out of control"),"and"),
                            p("2. greater than the selected # of metrics are", strong("out of control."))
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
                                 p(strong("# out of control metrics: ")),
                                 uiOutput("metricThresholdRed")
                            )
                          )
                        ),

                        wellPanel(
                          fluidRow(
                            p(strong("YELLOW FLAG"), style="color:black; background-color: yellow;",align = "center",style="font-size:125%;"),
                            p(strong("System performance is POOR when:"),align = "center"),
                            p("1. greater than the selected % of peptides are", strong("out of control"),"and"),
                            p("2. greater than the selected # of metrics are", strong("out of control.")),
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
              tabPanel("Metric summary",
                       tabsetPanel(

                         tabPanel("Descriptives : boxplots for metrics",
                                  tags$head(tags$style(type="text/css")),
                                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                   tags$div("It may take a while to load the plots, please wait...",
                                                            id="loadmessage")),
                                  plotlyOutput("box_plot", height = 2000)
                         ),

                         tabPanel("Overall performance : decision maps",
                                  tags$head(tags$style(type="text/css")),
                                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                   tags$div("It may take a while to load the plots, please wait...",
                                                            id="loadmessage")),
                                  sidebarLayout(
                                    sidebarPanel(
                                      checkboxGroupInput("heatmap_controlChart_select", "Select your control chart",
                                                         choices = c("CUSUM charts" = "CUSUM","XmR chart" = "XmR"), selected = "XmR")
                                      #htmlOutput("heatmap_txt")
                                    ),
                                    mainPanel(plotOutput("heat_map")
                                    )
                                  )
                         ),

                           tabPanel("Detailed performance: plot summaries",
                                    tags$head(tags$style(type="text/css")),
                                    conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                     tags$div("It may take a while to load the plots, please wait...",
                                                              id="loadmessage")),
                                    sidebarLayout(
                                      sidebarPanel(
                                        checkboxGroupInput("summary_controlChart_select", "Select your control chart",
                                                           choices = c("CUSUM charts" = "CUSUM","XmR chart" = "XmR"), selected = "XmR")
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
              navbarMenu("Control charts",
                         tabPanel("XmR control charts",
                                  uiOutput("XmR_tabset")
                                  ),

                         tabPanel("CUSUMm and CUSUMv control charts",
                                  uiOutput("CUSUM_tabset")
                                  ),
                         tabPanel("Change point analysis for mean and variability",
                                  uiOutput("CP_tabset")
                                  )
              ),
###################################################################################################
              tabPanel("Help",
                       tabsetPanel(
                         tabPanel("Metrics"
                                  ,h5(strong("Retention time")),
                                  p("Retention time is the time it takes a solute to travel through the column. The retention time is assigned to
                                    the corresponding solute peak. The retention time is a measure of the amount of time a solute spends in a column.
                                    It is the sum of the time spent in the stationary phase and the mobile phase."
                                    ,a("visit for more info",href="http://www.britannica.com/science/retention-time")),
                                  br(),
                                  h5(strong("Total peak area")),
                                  p("Total peak area is the sum of all integrated signals for a certain peptide."),
                                  br(),
                                  h5(strong("Full width at half maximum (FWHM)")),
                                  p("Full width at half maximum 'FWHM' is an expression of the extent of a
                                    function given by the difference between the two extreme values
                                    of the independent variable at which the dependent variable is equal
                                    to half of its maximum value."
                                    ,a("visit for more info",href="https://en.wikipedia.org/wiki/Full_width_at_half_maximum")),
                                  br(),
                                  h5(strong("Peak assymetry")),
                                  p("Peak assymetry is a measure of symetry for a peak. Calculated by taking 2*a/(a+b). Optimal value is around 1 for a Gaussian peak.")
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
                                  h5(strong("Change point analysis")),
                                  h5("Can identify the exact time of a change in the mean and dispersion of suitability metric. "),
                                  h5("Likelihood functions are plotted and the sample which maximizes the functions is considered as a candidate change point."),
                                  p("A change in the process parameters triggers a control chart to generate an out of
                                    control signal. The sample at which the signal is issued is considered as the
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
  ))
