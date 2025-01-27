library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyFiles)
library(DT)
library(ggplot2)
library(plotly)
library(magrittr)
library(purrr)
library(RColorBrewer)
library(imputeLCMD)
library(impute)
library(pcaMethods)
library(ggalt)
library(ggrepel)
library(ComplexHeatmap)
library(svglite)
library(ragg)
library(stringr)
library(rmarkdown)
library(data.table)
library(gtools)
library(dplyr)
library(tidyr)
library(scales)
library(rstudioapi)
library(reshape2)

####UI----
options(shiny.maxRequestSize=1000*1024^2)

ui <- dashboardPage(
  dashboardHeader(title = "Proteom'IC"),
  dashboardSidebar( width = 130,
                    sidebarMenu(id = "tabs",
                                menuItem("Data", tabName = "ImportData"),
                                menuItem("Descriptives Plots", tabName = "DescPlot"),
                                menuItem("Statistical Plots", tabName = "StatPlot"),
                                menuItem("Proteins Plots",tabName="ProtPlot"),
                                downloadButton("MyKnit", "Report",style="color: #000")
                    )
  ),
  dashboardBody(
    tabItems(
      ##### Read data ----
      tabItem(tabName = "ImportData",
              tabBox( width = 0,
                      tabPanel("Import",
                               fluidRow(
                                 box(title="1. Required files",
                                     width=4,
                                     h5(HTML("<b>For Maxquant analysis, folder must contain :</b>")),
                                     h6("   - proteinGroups.txt,"),
                                     h6("   - peptides.txt,"),
                                     h6("   - mqpar.xml,"),
                                     h6("   - evidence.txt,"),
                                     h6("   - #runningTimes.txt "),
                                     h5(HTML("<b>For DIANN analysis, folder must contain :</b>")),
                                     h6("   - report.pg_matrix.tsv.txt,"),
                                     h6("   - report.log.txt,"),
                                     h6("   - report.pr_matrix.tsv,"),
                                     h6("   - report.tsv ")),
                                 box(title="2. Select proteinGroup.txt or report.pg_matrix.tsv.",
                                     width=4,
                                     actionButton("MyFile","Browse..."),
                                     textInput("MyPath",label="",value=),
                                     textOutput("MySoft")),
                                 box(title="3. Expression columns",
                                     width=4,
                                     textInput("MyRegEx",label="Enter the regular expression to detect expression columns",value="(.*)"),
                                 )
                               ),
                               dataTableOutput("Mydf")
                      ),
                      
                      tabPanel("Experimental Design",
                               fluidRow(
                                 box(title = "1. Enter the list of condition separate by ';'",
                                     width = 3,
                                     textInput("MyCond",label=""),
                                     h5(HTML("<font color='red'>Warning :</font>")),
                                     h5(HTML("<font color='red'>If ExperimentalDesign.txt already exists and if you want to change the conditions list, please delete or change the name of this file.</font>"))),
                                 box(title = "2. Choose color for each condition",
                                     width = 3,
                                     uiOutput("colors")),
                                 box(title = "3. Select condition for each sample",
                                     width = 3,
                                     uiOutput("SampleCond")),
                                 box(title="4. Download ExperimentalDesign",
                                     width=3,
                                     downloadButton("ExportPhenodata",label="Download"))
                               ),
                               dataTableOutput("phenodata")
                      ),
                      
                      tabPanel("Statistics",
                               selectInput(
                                 inputId = "MyStat",
                                 label = "Make T-test ?",
                                 choices  = c("YES",'NO'),
                                 selected="YES"
                               ),
                               conditionalPanel(condition = "input.MyStat == 'YES'",
                                                fluidRow(
                                                  box(title="1. Comparisons",
                                                      width=4,
                                                      selectizeInput(inputId="MyComp",
                                                                     label="Select comparisons to test :",
                                                                     choices = "",
                                                                     selected = NULL,
                                                                     multiple = TRUE)
                                                  ),
                                                  box(title="2. Type of T-test",
                                                      width=4,
                                                      selectInput(inputId = "MyType",
                                                                  label = "",
                                                                  choices = paste(rep(c("two.sided", "less", "greater"),Times=2),rep(c("Unpaired","Paired"),each=3)),
                                                                  selected = "two.sided Unpaired")
                                                  ),
                                                  box(title="3. Threasholds",
                                                      width=4,
                                                      selectInput(inputId = "ThreasholdType",
                                                                  label = "Select the type of threashold :",
                                                                  choices = c("Pvalue","Qvalue"),
                                                                  selected = "Pvalue"),
                                                      numericInput(inputId="SignifThreashold",
                                                                   label="Choose the significant P/Q value threashold :",
                                                                   value=0.05),
                                                      numericInput(inputId="FCThreashold",
                                                                   label="Choose the Fold Change threashold (1 means no difference) :",
                                                                   value=1))
                                                )),
                               actionButton("Calculatefdata","OK",width="100%"),
                               progressBar(id = "pb1", value = 0,title = "Process 0/3"),
                               dataTableOutput("TableNbProtSignif"),
                               dataTableOutput("TableAppearedDsiappeared"),
                               dataTableOutput("fdata")
                      )
              )
      ),
      #####Descriptiveplot ----
      tabItem(tabName = "DescPlot",
              
              tabBox( width = 0,
                      ###### NbProt ----
                      tabPanel("Nbprot",
                               plotlyOutput('nbprot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "NbProtformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthNbProt", "Width of downloaded plot (cm)", 20),
                               numericInput("HeightNbProt", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadNbProtplot", "Download")
                      ),
                      ###### ProfilPlot ----
                      tabPanel("ProfilPlots",
                               plotlyOutput('ProfilConta',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "ProfilContaformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthProfilConta", "Width of downloaded plot (cm)", 20),
                               numericInput("HeightProfilConta", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadProfilConta", "Download"),
                               
                               plotlyOutput('ProfilTop5',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "ProfilTop5format",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthProfilTop5", "Width of downloaded plot (cm)", 20),
                               numericInput("HeightProfilTop5", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadProfilTop5", "Download")
                      ),
                      ###### PCA ----
                      tabPanel("PCA",  
                               selectInput(
                                 inputId = "TypeDataPCA",
                                 label = "Type of data used to generate PCA",
                                 choices = c("Z-score(intensity)", "Log2(intensity)"),
                                 selected = "Z-score(intensity)"
                               ),
                               selectInput(
                                 inputId = "FilterType",
                                 label = "Type of filter used to generate PCA",
                                 choices = c("% of Valid value in at least one group", "% of Valid value in each group","% of Valid value in total"),
                                 selected = "% of Valid value in at least on group"
                               ),
                               numericInput("PcFilter","% of valid values used to generate PCA",70),
                               actionButton("DrawPCA","Draw PCA"),
                               helpText("WARRNING : PCA take a long time to be displayed."),
                               textOutput("NbProtPCA"),
                               plotOutput('PCA',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "PCAformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthPCA", "Width of downloaded plot (cm)", 26),
                               numericInput("HeightPCA", "Height of downloaded plot (cm)", 17),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadPCAplot", "Download"),
                               
                               plotlyOutput('PCAprotInteractive',width = "100%", height = "400"),
                               selectizeInput(
                                 inputId = "GenePCA",
                                 label = "Select which genes to show",
                                 choices = NULL,
                                 multiple = TRUE
                               ),
                               numericInput("ForcePCA", "Force to avoid overlap of gene names", 1),
                               plotOutput('PCAprot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "PCAprotformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthPCAprot", "Width of downloaded plot (cm)", 26),
                               numericInput("HeightPCAprot", "Height of downloaded plot (cm)", 17),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadPCAprotplot", "Download")
                      ),
                      
                      ###### Pearson Correlations ----
                      tabPanel("Pearson correlations", 
                               checkboxInput(
                                 inputId = "PearsonCluster",
                                 label = "Add Sample clustering",
                                 value = FALSE
                               ),
                               plotOutput('Pearson',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "Pearsonformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthPearson", "Width of downloaded plot (cm)", 20),
                               numericInput("HeightPearson", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadPearsonplot", "Download"),
                               br(),
                               plotOutput('PearsonBox',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "PearsonBoxformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthPearsonBox", "Width of downloaded plot (cm)", 20),
                               numericInput("HeightPearsonBox", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadPearsonBoxplot", "Download")
                      ),
                      
                      ###### Euclidean distances ----
                      tabPanel("Euclidean distances", 
                               checkboxInput(
                                 inputId = "EuclidCluster",
                                 label = "Add Sample clustering",
                                 value = FALSE
                               ),
                               plotOutput('Euclid',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "Euclidformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthEuclid", "Width of downloaded plot (cm)", 20),
                               numericInput("HeightEuclid", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadEuclidplot", "Download"),
                               br(),
                               plotOutput('EuclidBox',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "EuclidBoxformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthEuclidBox", "Width of downloaded plot (cm)", 20),
                               numericInput("HeightEuclidBox", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadEuclidBoxplot", "Download")
                      ),
                      ###### Heatmap Log2(intensity) ----
                      tabPanel("Heatmap Log2(intensity)", 
                               colourpicker::colourInput(
                                 inputId = "LowColor",
                                 label = "Choose the color of low expression proteins",
                                 value = "#7AFFC3"
                               ),
                               colourpicker::colourInput(
                                 inputId = "MiddleColor",
                                 label = "Choose the color of middle expression proteins",
                                 value = "#056660"
                               ),
                               colourpicker::colourInput(
                                 inputId = "HighColor",
                                 label = "Choose the color of high expression proteins",
                                 value = "#000000"
                               ),
                               checkboxInput(
                                 inputId = "SampleCluster",
                                 label = "Add Sample clustering",
                                 value = FALSE
                               ),
                               plotOutput('HeatmapLog2',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "HeatmapLog2format",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthHeatmapLog2", "Width of downloaded plot (cm)", 15),
                               numericInput("HeightHeatmapLog2", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadHeatmapLog2plot", "Download")
                      ),
                      
                      ###### Scaterplot ----
                      tabPanel("Scater plot", 
                               h5(HTML("<font color='red'>Warning :</font>")),
                               h5(HTML("<font color='red'>This plot is not exported into the report.</font>")),
                               selectInput(
                                 inputId = "Sample1",
                                 label = "Sample on X axis",
                                 choices = NULL,
                                 selected = NULL
                               ),
                               selectInput(
                                 inputId = "Sample2",
                                 label = "Sample on Y axis",
                                 choices = NULL,
                                 selected = NULL
                               ),
                               textOutput('CorScaterplot'),
                               actionButton("DrawScaterPlot","Draw Scaterplot"),
                               plotlyOutput('Scaterplotinteractive',width = "100%", height = "400"),
                               selectizeInput(
                                 inputId = "GenesScaterplot",
                                 label = "Select genes to show",
                                 choices = NULL,
                                 multiple = TRUE
                               ),
                               numericInput("ForceScaterplot", "Force to avoid overlap of gene names", 1),
                               plotOutput('Scaterplot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "ScaterPlotformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthScaterPlot", "Width of downloaded plot (cm)", 15),
                               numericInput("HeightScaterPlot", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadScaterPlot", "Download")
                      ),
                      ###### MeanScaterplot ----
                      tabPanel("Mean Scater plot",
                               h5(HTML("<font color='red'>Warning :</font>")),
                               h5(HTML("<font color='red'>This plot is not exported into the report.</font>")),
                               selectInput(
                                 inputId = "Cond1",
                                 label = "Condition on X axis",
                                 choices = NULL,
                                 selected = NULL
                               ),
                               selectInput(
                                 inputId = "Cond2",
                                 label = "Condition on y axis",
                                 choices = NULL,
                                 selected = NULL
                               ),
                               textOutput('CorMeanScaterplot'),
                               actionButton("DrawMeanScaterPlot","Draw Scaterplot"),
                               plotlyOutput('MeanScaterplotinteractive',width = "100%", height = "400"),
                               selectizeInput(
                                 inputId = "GenesMeanScaterplot",
                                 label = "Select genes to show",
                                 choices = NULL,
                                 multiple = TRUE
                               ),
                               numericInput("ForceMeanScaterplot", "Force to avoid overlap of gene names", 1),
                               plotOutput('MeanScaterplot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "MeanScaterPlotformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthMeanScaterPlot", "Width of downloaded plot (cm)", 15),
                               numericInput("HeightMeanScaterPlot", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadMeanScaterPlot", "Download")
                      )
              )
      ),
      # #####Statisticalplot ----
      tabItem(tabName = "StatPlot",
              selectInput(
                inputId = "ComparisonSelect",
                label = "Select the comparison",
                choices = NULL,
                selected = NULL
              ),
              textOutput("NbProtSignif"),
              tabBox( width = 0,
                      
                      ###### PvaluePlot ----
                      tabPanel("Pvalues",
                               plotOutput("PvalPlot",width = "100%", height = "400"),
                               selectInput(
                                 inputId = "PvalPlotformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthPvalPlot", "Width of downloaded plot (cm)", 20),
                               numericInput("HeightPvalPlot", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadPvalPlot", "Download")),
                      ###### FCPlot ----
                      tabPanel("FC",
                               plotOutput("FCplot",width = "100%", height = "400"),
                               selectInput(
                                 inputId = "FCPlotformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthFCPlot", "Width of downloaded plot (cm)", 20),
                               numericInput("HeightFCPlot", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadFCPlot", "Download")),
                      ###### VolcanoPlot ----
                      tabPanel("Volcanoplot",
                               fluidRow(
                                 box(title = "Colors",
                                     width = 4,
                                     colourpicker::colourInput(
                                       inputId = "UpColor",
                                       label = "Choose the color of up regulated proteins",
                                       value = "#FF4040"
                                     ),
                                     colourpicker::colourInput(
                                       inputId = "DownColor",
                                       label = "Choose the color of down regulated proteins",
                                       value = "#4C4CFF"
                                     ),
                                     colourpicker::colourInput(
                                       inputId = "NSColor",
                                       label = "Choose the color of non significant proteins",
                                       value = "#C4C4C4"
                                     )    
                                 ),
                                 box(title="Lines of threashold",
                                     width=4,
                                     checkboxInput(
                                       inputId = "ShowHline",
                                       label = "Add P/Qvalue threashhold line",
                                       value = FALSE
                                     ),
                                     checkboxInput(
                                       inputId = "ShowVline",
                                       label = "Add Log2(FC) threashhold line",
                                       value = FALSE
                                     ),
                                     selectInput(
                                       inputId = "LineType",
                                       label = "Select the type of the line",
                                       choices = c( "solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678"),
                                       selected = "solid"
                                     ),
                                     colourpicker::colourInput(
                                       inputId = "LineColor",
                                       label = "Choose the color of the line",
                                       value = "#000000"
                                     )
                                 )
                               ),
                               plotlyOutput('IntercativeVolcanoPlot',width = "100%", height = "400"),
                               fluidRow(
                                 box(title="Gene names",
                                     width=4,
                                     selectizeInput(
                                       inputId = "sel_gene_nm",
                                       label = "Select which gene names to show",
                                       choices = NULL,
                                       multiple = TRUE
                                     ),
                                     numericInput("Force", "Force to avoid overlap of gene names", 1),
                                 ),
                                 box(title="Gene points",
                                     width=4,
                                     selectizeInput(
                                       inputId = "GenePoints",
                                       label = "Select which gene points to show",
                                       choices = NULL,
                                       multiple = TRUE
                                     ),
                                     numericInput("PointSize", "Points size", 1),
                                     colourpicker::colourInput(
                                       inputId = "PointColor",
                                       label = "Points Color",
                                       value = "#000000"
                                     )
                                     
                                 )
                               ),
                               plotOutput('VolcanoPlot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "Volcanoformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthVolcano", "Width of downloaded plot (cm)", 15),
                               numericInput("HeightVolcano", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadVolcanoplot", "Download"),
                               
                               
                      ),
                      ###### Heatmap Significant----
                      tabPanel("Heatmap",
                               colourpicker::colourInput(
                                 inputId = "LowColorSignif",
                                 label = "Choose the color of low expression proteins",
                                 value = "#4C4CFF"
                               ),
                               colourpicker::colourInput(
                                 inputId = "MiddleColorSignif",
                                 label = "Choose the color of middle expression proteins",
                                 value = "#FFFFFF"
                               ),
                               colourpicker::colourInput(
                                 inputId = "HighColorSignif",
                                 label = "Choose the color of high expression proteins",
                                 value = "#FF4040"
                               ),
                               checkboxInput(
                                 inputId = "ShowGeneName",
                                 label = "Show gene names",
                                 value = FALSE
                               ),
                               plotOutput('HMSignif',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "HMSignifformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthHMSignif", "Width of downloaded plot (cm)", 15),
                               numericInput("HeightHMSignif", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadHMSignif", "Download")
                      ),
                      ###### Heatmap Anova----
                      tabPanel("Anova",
                               textOutput("NbProtSignifAnova"),
                               colourpicker::colourInput(
                                 inputId = "LowColorSignifAnova",
                                 label = "Choose the color of low expression proteins",
                                 value = "#4C4CFF"
                               ),
                               colourpicker::colourInput(
                                 inputId = "MiddleColorSignifAnova",
                                 label = "Choose the color of middle expression proteins",
                                 value = "#FFFFFF"
                               ),
                               colourpicker::colourInput(
                                 inputId = "HighColorSignifAnova",
                                 label = "Choose the color of high expression proteins",
                                 value = "#FF4040"
                               ),
                               checkboxInput(
                                 inputId = "ShowGeneNameAnova",
                                 label = "Show gene names",
                                 value = FALSE
                               ),
                               plotOutput('HMAnova',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "HMAnovaformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthHMAnova", "Width of downloaded plot (cm)", 15),
                               numericInput("HeightHMAnova", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadHMAnova", "Download")
                      )
              )
      ),
      tabItem("ProtPlot",
              selectInput(
                inputId = "ProtPlot_GeneName",
                label = "Select the Gene to show",
                choices = NULL,
                selected = NULL
              ),
              tabBox( width = 0,
                      
                      ###### ProtPlot ----
                      tabPanel("Protein level",
                               h5(HTML("<font color='red'>Warning :</font>")),
                               h5(HTML("<font color='red'>This plot is not exported into the report.</font>")),
                               textInput("ytitle", "y axis title", "Intensity"),
                               selectInput(
                                 inputId = "BarOrLine",
                                 label = "Select the type of plot display",
                                 choices = c("Bar plot", "Line plot"),
                                 selected = "Bar plot"
                               ),
                               plotlyOutput('ProtBarplot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "ProtBarPlotformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthProtBarPlot", "Width of downloaded plot (cm)", 10),
                               numericInput("HeightProtBarPlot", "Height of downloaded plot (cm)", 10),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadProtBarPlot", "Download")
                      ),
                      ###### PepPlot ----
                      tabPanel("Peptide level",
                               h5(HTML("<font color='red'>Warning :</font>")),
                               h5(HTML("<font color='red'>This plot is not exported into the report.</font>")),
                               plotlyOutput('PepLinePlot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "PepLinePlotformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               numericInput("WidthPepLinePlot", "Width of downloaded plot (cm)", 30),
                               numericInput("HeightPepLinePlot", "Height of downloaded plot (cm)", 15),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadPepLinePlot", "Download")
                      )
                      
              )
      )
    )
  )
)


####Server ----
server <- function(input, output, session) {
  
  Pepdata<-eventReactive(input$dataFile,{
    out <- readRDS(input$dataFile$datapath)
    out<-out[["Pepdata"]]
    out[,3:ncol(out)]<-log2(out[,3:ncol(out)])
    out
  })
  
  
  #####Functions  ----  
  
  Exportggplot<-function(graph,basename,format,width,height){
    downloadHandler(
      filename = function() {
        paste0(basename,".", req(format))
      },
      content = function(file) {
        ggsave(file,
               device=req(format),
               plot=graph,
               dpi=600,
               width=width,
               height=height,
               units="cm")
      })
  }
  
  Exportheatmap<-function(graph,basename,format,width,height){
    downloadHandler(
      filename = function() {
        paste0(basename,".", req(format))
      },
      content = function(file) {
        if (format == "png") {
          ragg::agg_png(file,
                        res = 600,
                        width=width,
                        height=height, 
                        units = "cm"
          )
        } else if (format == "pdf") {
          pdf(file,
              width=width*0.393701,
              height=height*0.393701)
        } else if (format == "svg") {
          svglite::svglite(file,
                           width=width/2,
                           height =height/2)
        }
        ComplexHeatmap::draw(graph)
        dev.off()
      }
    )}
  
  Density<-function(x,y,points){
    #CalcDensityOnGrid
    MyData<-data.frame(x=x,y=y)
    MyData<-MyData[!is.na(MyData[,1])&!is.na(MyData[,2]),]
    
    xmin<-min(MyData[,1])
    xmax<-max(MyData[,1])
    ymin<-min(MyData[,2])
    ymax<-max(MyData[,2])
    
    #GetValuesOnGrid
    xvals<-MyData[,1]
    xStep<-(xmax - xmin) / points
    yvals<-MyData[,2]
    yStep<-(ymax - ymin) / points
    
    n = length(xvals)
    
    #CalcCovariance
    CalcCovariance<-function(data){
      n = nrow(data)
      p = ncol(data)
      
      means<-apply(data,2,mean)
      
      cov = matrix(0,nrow=p,ncol=p)
      for (i in 1:p){
        for (j in 1:i){
          cov[i, j] = sum((data[, i] - means[i]) * (data[, j] - means[j]))
          cov[i, j] = cov[i, j]/ n
          cov[j, i] = cov[i, j]
        }
      }
      return(cov)
    }
    
    cov<-CalcCovariance(MyData)#OK
    
    fact = n^(1 / 6)#OK
    hinv = fact/sqrt(abs(cov))
    hinv[c(1,2), 1] <-hinv[c(1,2),1]*xStep#OK
    hinv[c(1,2), 2] <-hinv[c(1,2),2]*yStep#OK
    
    dx = as.integer(1.0 / hinv[1, 1] * 5)#OK
    dy = as.integer(1.0 / hinv[2, 2] * 5)#OK
    
    values = matrix(0,ncol=points,nrow = points)
    
    xind = as.integer(floor((xvals - xmin) / xStep))#OK
    yind = as.integer(floor((yvals - ymin) / yStep))#OK
    
    for (i in 1:n){
      
      ii <-max(xind[i] - dx, 1): min(xind[i] + dx, points)
      jj <-max(yind[i] - dy, 1): min(yind[i] + dy, points)
      
      #MatrixTimesVector & StandardGaussian
      a1<-((hinv[1,1]*(ii - xind[i]))+(hinv[2,1]*(ii - xind[i])))^2
      a1<-matrix(rep(a1,length(jj)),nrow=length(ii))
      
      a2<-((hinv[1,2]*(jj - yind[i]))+(hinv[2,2]*(jj - yind[i])))^2
      a2<-matrix(rep(a2,each=length(ii)),ncol=length(jj))
      
      MySum<-a1+a2
      StandardGaussian<-exp(-0.5*MySum)/(2*pi)#OK
      
      values[ii, jj] <-values[ii, jj]+StandardGaussian
      
    }
    values<-as.data.frame(values)
    values[values == "NaN"]<-0 
    
    values<-values/max(values)#OK
    xmat <- seq(xmin, xmax, length.out = points)
    ymat <- seq(ymin, ymax, length.out = points)
    
    dvals = rep(0,length(x))
    for (i in 1:length(x)){
      xx = x[i]
      yy = y[i]
      if (!is.na(xx) && !is.na(yy)){
        xind = length(xmat[xmat<=xx])
        yind = length(ymat[ymat<=yy])
        dvals[i] = values[xind, yind]
      } else{
        dvals[i] = NA
      }
    }
    return(dvals)
  }
  
  DataScaterPlot<-function(x,y,edata,fdata){
    MyData<-data.frame(x=x,y=y)
    MyData$Density<-Density(x,y,300)
    MyData$edataRowNames<-row.names(edata)
    MyData<-merge(MyData,fdata,by.x="edataRowNames",by.y="RowNamesfdata",all.x=TRUE)
    return(MyData)
  }
  
  DrawScaterplot<-function(MyData,TitleX,TitleY){
    
    hmcol<- c(colorRampPalette(c("green","yellow"))(5),
              colorRampPalette(c("yellow","orange"))(5),
              colorRampPalette(c("orange","red"))(15),
              colorRampPalette(c("red","blue"))(25),
              colorRampPalette(c("blue","cyan"))(50))
    
    ggplot(data=MyData,aes_string(x=colnames(MyData)[2],y=colnames(MyData)[3],color="Density",label="Gene_ProteinGroup"))+
      geom_point(size=0.5)+
      scale_color_gradientn(colours = hmcol)+
      geom_abline(slope=1,intercept = 0,colour ="white",linewidth=1)+
      theme(panel.background = element_rect(fill = "black"),
            panel.grid.major=element_line(colour="grey40"),
            panel.grid.minor=element_line(colour="grey40"))+
      xlab(TitleX) +
      ylab(TitleY)
  }
  
  #####Import  ----

  MyFile <-eventReactive(input$MyFile,{
    MyPath<-selectFile(caption="Select report.pg_matrix.tsv or proteinGroup.txt file.",label="OK")

    df<-as.data.frame(fread(MyPath,header = T, sep = "\t"))
    MySoft<-if(str_detect(string=MyPath,pattern="proteinGroups.txt")==TRUE){"MaxQ"}else{"DIANN"}
    if(MySoft=="MaxQ"){
      df <- df[df$Reverse == ''| is.na(df$Reverse),]
      df <- df[df$'Only identified by site' == '' | is.na(df$'Only identified by site'),] 
      colnames(df)[1]<-"Protein.Group"
      colnames(df)[which(colnames(df)=="Gene names")]<-"Genes"
    }else{
      df<-df[,c(1:5)]
    } 
    
    out<-paste(df$Genes,df$Protein.Group)
    
    updateSelectizeInput(
      inputId = "sel_gene_nm",
      choices = out,
      server = TRUE,
      selected = NULL
    )
    updateSelectizeInput(
      inputId = "GenePoints",
      choices = out,
      server = TRUE,
      selected = NULL
    )
    updateSelectizeInput(
      inputId = "ProtPlot_GeneName",
      choices = out,
      server = TRUE,
      selected = NULL
    )
    updateSelectizeInput(
      inputId = "GenePCA",
      choices = out,
      server = TRUE,
      selected = NULL
    )
    updateSelectizeInput(
      inputId = "GenesScaterplot",
      choices = out,
      server = TRUE,
      selected = NULL
    )
    updateSelectizeInput(
      inputId = "GenesMeanScaterplot",
      choices = out,
      server = TRUE,
      selected = NULL
    )
    
    MyPath
  })
  
  observeEvent(input$MyFile,{
    updateTextInput(session, "MyPath", label="",
                    value = MyFile())
  })
  

  
  MySoft<-reactive({
    if(str_detect(string=MyFile(),pattern="proteinGroups.txt")==TRUE){"MaxQ"}else{"DIANN"}
  })
  
  output$MySoft<-renderText({paste("Software :",MySoft())}) 
  
  observe({
    if (MySoft()=="MaxQ"){
      Label<-"3. Enter the generic expression to detect columns to work with"
      Value<-"LFQ intensity "
    }else{
      Label<-"3. Enter the regular expression to extract Sample names"
      Value<-"(.*)"
    }
    updateTextInput(session, "MyRegEx", 
                    label=Label, 
                    value = Value)
  })
  
  generatedf<-function(){
    MyPath<-input$MyPath
    MyRegEx<-input$MyRegEx
    df<-as.data.frame(fread(MyPath,header = T, sep = "\t"))
    MySoft<-if(str_detect(string=MyPath,pattern="proteinGroups.txt")==TRUE){"MaxQ"}else{"DIANN"}
    
    if(MySoft=="MaxQ"){
      df <- df[df$Reverse == ''| is.na(df$Reverse),]
      df <- df[df$'Only identified by site' == '' | is.na(df$'Only identified by site'),] 
      x<-df[,str_detect(colnames(df),MyRegEx)==TRUE]
      df<-df[,c(1:2,6:8,which(colnames(df)=="Potential contaminant"))]
      df<-cbind(df,x)
      colnames(df)[1]<-"Protein.Group"
      colnames(df)[4]<-"Genes"
      colnames(df)[6]<-"Potential.contaminant"
      
      x<-colnames(df)[7:ncol(df)]
      x<-gsub(x=x,pattern=MyRegEx,replacement="")
      colnames(df)[7:ncol(df)]<-x
    }else{
      if(length(which(colnames(df)=="Protein.Ids"))==0){
        df<-df[,c(1,1,2:ncol(df))]
        colnames(df)[2]<-"Protein.Ids"
      }
      df$Potential.contaminant<-str_detect(df$Protein.Group,"Conta_")
      df[df$Potential.contaminant==TRUE,"Potential.contaminant"]<-"+"
      df[df$Potential.contaminant==FALSE,"Potential.contaminant"]<-""
      df<-df[,c(1:5,ncol(df),6:(ncol(df)-1))]
      
      x<-colnames(df)[7:ncol(df)]
      x<-str_extract(string=x,pattern=".*\\\\(.*)\\.d$",group=1)
      colnames(df)<-c(colnames(df)[1:6],x)
      
      x<-colnames(df)[7:ncol(df)]
      x<-str_extract(string=x,pattern=MyRegEx,group=1)
      colnames(df)[7:ncol(df)]<-x
    }  
    return(df)
  }
  
  df<-reactive({generatedf()})
  
  output$Mydf <- renderDataTable({
    datatable(generatedf(),
              options = list(
                pageLength=5, scrollX='400px'))
    
  })
  
  
  
  edata<-reactive({
    df<-generatedf()
    df<-df[df$Potential.contaminant=="",]
    edata <- df[,7:length(colnames(df))]
    
    edata <- log2(edata)
    edata[edata == -Inf]<-NA
    
    edata <- edata[rowSums(is.na(edata)) != ncol(edata), ]
    edata
  })
  
  
  
  
  
  #####Experimental Design  ----
  MyCond<-reactive({
    unlist(str_split(input$MyCond,";"))
  })
  
  output$SampleCond <- renderUI({
    df<-generatedf()
    MySample<-colnames(df)[7:ncol(df)]
    
    x<-as.list(MySample)
    names(x)<-MySample
    updateSelectInput(inputId="Sample1",choices=x,selected=x[[1]])
    updateSelectInput(inputId="Sample2",choices=x,selected=x[[2]])

    if(file.exists(paste0(dirname(input$MyPath),"/ExperimentalDesign.txt"))==TRUE){
      phenodata<-as.data.frame(fread(paste0(dirname(input$MyPath),"/ExperimentalDesign.txt"),
                                     header = T, 
                                     sep = "\t"))
      
      x<-paste(unique(phenodata$Condition),collapse=";")
      updateTextInput(inputId="MyCond",value=x)
      
      MyCond<-as.list(unique(phenodata$Condition))
      names(MyCond)<-unique(phenodata$Condition)
      
      updateSelectInput(inputId = "Cond1",choices = MyCond,selected=MyCond[[1]])
      updateSelectInput(inputId = "Cond2",choices = MyCond,selected=MyCond[[2]])
      
      MyComparison<-unlist(str_split(input$MyCond,";"))
      MyComparison<-permutations(n=length(MyComparison),r=2,v=MyComparison)
      MyComparison<-paste(MyComparison[,1],MyComparison[,2],sep="/")
      MyComparison<-as.list(MyComparison)
      names(MyComparison)<-unlist(MyComparison)
      updateSelectizeInput(inputId="MyComp",
                           choices=MyComparison)
      updateSelectInput(inputId = "ComparisonSelect",
                        choices = MyComparison)
      
      lapply(1:nrow(phenodata),function(x){
        selectInput(inputId=phenodata$ShortName[x],
                    label = phenodata$ShortName[x],
                    choices = MyCond,
                    selected=phenodata$Condition[x])
      })
      
    }else{
      MyCond<-as.list(unlist(str_split(input$MyCond,";")))
      names(MyCond)<-unlist(str_split(input$MyCond,";"))
      
      updateSelectInput(inputId = "Cond1",choices = MyCond,selected=MyCond[[1]])
      updateSelectInput(inputId = "Cond2",choices = MyCond,selected=MyCond[[2]])
      
      MyComparison<-unlist(str_split(input$MyCond,";"))
      MyComparison<-permutations(n=length(MyComparison),r=2,v=MyComparison)
      MyComparison<-paste(MyComparison[,1],MyComparison[,2],sep="/")
      MyComparison<-as.list(MyComparison)
      names(MyComparison)<-unlist(MyComparison)
      updateSelectizeInput(inputId="MyComp",
                           choices=MyComparison)
      updateSelectInput(inputId = "ComparisonSelect",
                        choices = MyComparison)
      
      lapply(MySample, function(x){
        selectInput(x,
                    label = x,
                    choices = MyCond,
                    selected=NULL)
      })
    }
    
    
  })
  
  output$colors <- renderUI({
    
    MyCond<-unlist(str_split(input$MyCond,";"))

    if(file.exists(paste0(dirname(input$MyPath),"/ExperimentalDesign.txt"))==TRUE){
      phenodata<-as.data.frame(fread(paste0(dirname(input$MyPath),"/ExperimentalDesign.txt"),
                                     header = T, 
                                     sep = "\t"))
      
      MyColors<-c()
      for(i in MyCond){
        x<-phenodata[phenodata$Condition==i,"CondColor"]
        x<-unique(as.vector(x))
        MyColors<-c(MyColors,x)
      }
      
      purrr::map2(
        MyCond,
        MyColors,
        ~ colourpicker::colourInput(
          inputId = session$ns(.x),
          paste("Color", .x),
          value = .y
        )
      )
      
    }else{
      purrr::map2(
        MyCond,
        colorRampPalette(brewer.pal(9, "Set1"))(if(length(MyCond())>9){length(MyCond())}else{9})[1:length(MyCond())],
        ~ colourpicker::colourInput(
          inputId = session$ns(.x),
          paste("Color", .x),
          value = .y
        )
      )
    }
  })
  
  generatephenodata<-function(){
    df<-generatedf()
    MySample<-colnames(df)[7:ncol(df)]
    
    phenodata<-data.frame(ShortName=MySample,Condition=MySample,CondColor=MySample)
    for(i in MySample){
      phenodata[phenodata$ShortName==i,"Condition"]<- input[[i]]
    }
    
    MyCond<-unlist(str_split(input$MyCond,";"))
    
    for(i in MyCond){
      phenodata[phenodata$Condition==i,"CondColor"]<- input[[i]]
    }
    
    phenodata
    
  }
  
  output$phenodata<- renderDataTable({
    datatable(generatephenodata(),
              options=list(paging=FALSE))
  })
  
  phenodata<-reactive({
    generatephenodata()
  })
  
  output$ExportPhenodata<-downloadHandler(
    filename = function(){"ExperimentalDesign.txt"},
    content=function(file){
      write.table(phenodata(),
                  file,
                  row.names = FALSE,
                  sep = "\t")
    })
  
  ######Comparison############### 
  
  generatefdata<-function(){
    df<-generatedf()
    fdata <- df[,1:6]
    select_index <- match(row.names(edata()),row.names(fdata))
    fdata <- fdata[select_index,]
    fdata<-cbind(fdata,
                 RowNamesfdata=row.names(fdata),
                 Gene_ProteinGroup=paste(fdata$Genes,fdata$Protein.Group))
    
    if(input$MyStat=="NO"){
      fdata<-fdata[,c(ncol(fdata)-1,1:(ncol(fdata)-2),ncol(fdata))]
      fdata
    }else{
      
      MyComparison<-input$MyComp
      MyAlternative<-str_split(input$MyType," ")[[1]][1]
      MyPaired<-str_split(input$MyType," ")[[1]][2]=="Paired"
      MySignifType<-input$ThreasholdType
      MySignifVal<-input$SignifThreashold
      MySignifFC<-log2(input$FCThreashold)
      
      x<-as.list(input$MyComp)
      names(x)<-input$MyComp
      updateSelectInput(inputId = "ComparisonSelect",choices=x)
      
      
      #T-test
      TtestApply<-function(rank,x,MyListX=MyListX,MyListY=MyListY){
        
        MyTtest<-t.test(x=as.numeric(x[rank,MyListX]),
                        y=as.numeric(x[rank,MyListY]),
                        alternative = MyAlternative,
                        paired=MyPaired, 
                        var.equal=FALSE,
                        na.action="na.omit")
        MyPvalue<-MyTtest$p.value
        if(MyPaired==FALSE){
          MyFC<-MyTtest$estimate['mean of x']-MyTtest$estimate['mean of y']
        }else{
          MyFC<-MyTtest$estimate['mean difference']#'mean of the differences'
        }
        
        res<-c(MyPvalue,MyFC)
        return(res)
      }
      
      
      FilterPreTtest <- function(data, Comparaison, MinPc1Condition, MinNbAllCondition) {
        
        MyListX<-which(phenodata()$Condition == str_split(Comp,"/")[[1]][1])
        MyListY<-which(phenodata()$Condition == str_split(Comp,"/")[[1]][2])
        
        Mydata<-data[,c(MyListX,MyListY)]
        
        ValidProt<-rep(NA, nrow(Mydata))
        for(i in 1:nrow(Mydata)){
          
          nbVV <- tapply(as.numeric(Mydata[i,]), INDEX=phenodata()[c(MyListX,MyListY),"Condition"], FUN=complete.cases)#return TRUE FALSE Valid Values for each sample per group
          nbVV <- lapply(nbVV,sum)#nb VV per group
          
          PcVV<-nbVV
          
          for(condition in str_split(Comp,"/")[[1]]) {
            PcVV[[condition]] <- PcVV[[condition]]/table(phenodata()[c(MyListX,MyListY),"Condition"])[[condition]]*100#%VV per condition
          }
          
          ValidProt[i] <- any(PcVV >= MinPc1Condition) & nbVV[1]>=MinNbAllCondition & nbVV[2]>=MinNbAllCondition#return TRUE there is at least 1 group with 70%VV and at least 3VV in the other
        }
        
        return(ValidProt)
      }
      
      FilterPreTtestPaired <- function(data, Comparaison, MinPc1Condition) {
        
        MyListX<-which(phenodata()$Condition == str_split(Comp,"/")[[1]][1])
        MyListY<-which(phenodata()$Condition == str_split(Comp,"/")[[1]][2])
        
        Mydata<-data[,c(MyListX,MyListY)]
        
        ValidProt<-rep(NA, nrow(Mydata))
        for(i in 1:nrow(Mydata)){
          
          nbVV <- tapply(as.numeric(Mydata[i,]),
                         INDEX=phenodata()[c(MyListX,MyListY),"Condition"], 
                         FUN=complete.cases)#return TRUE FALSE Valid Values for each sample per group
          
          nbVV<-paste(nbVV[[1]],nbVV[[2]],sep = "_")
          nbVV<-length(nbVV[nbVV=="TRUE_TRUE"])
          PcVV<-100*nbVV/length(MyListX)#%Valid Values
          
          ValidProt[i] <- PcVV >= MinPc1Condition #return TRUE if there is at least 70% Valid values 
        }
        
        return(ValidProt)
      }
      
      for(Comp in MyComparison){
        #position of columns of the comparison  
        MyListX<-which(phenodata()$Condition == str_split(Comp,"/")[[1]][1])
        MyListY<-which(phenodata()$Condition == str_split(Comp,"/")[[1]][2])
        
        
        #Filter
        if(MyPaired==FALSE){
          DataToTest <- edata()[FilterPreTtest(data=edata(), Comparaison=Comp, MinPc1Condition=70, MinNbAllCondition=3),]
        }else{
          DataToTest <- edata()[FilterPreTtestPaired(data=edata(), Comparaison=Comp, MinPc1Condition=60),]
        }
        
        
        
        #T-test
        Myres<-mapply(TtestApply,1:nrow(DataToTest),MoreArgs=list(DataToTest,MyListX,MyListY))
        Myres<-as.data.frame(t(Myres))
        
        #add Qvalues
        Myres$Qvalue<-p.adjust(p=Myres[,1], method = "BH")
        
        #presentation of results
        Myres<-as.data.frame(cbind(Myres[,1],Myres[,3],Myres[,2]))
        colnames(Myres)<-paste(c("Pvalue","Qvalue","Log2"),Comp)
        Myres$RowNameDataToTest<-row.names(DataToTest)
        
        fdata<-merge(fdata,Myres,all.x=TRUE, all.y=FALSE,by.x="RowNamesfdata",by.y="RowNameDataToTest",sort=FALSE)
      }
      
      
      #VV
      VV<-as.data.frame(row.names(edata()))
      for(Comp in MyComparison){
        MyListX<-which(phenodata()$Condition == str_split(Comp,"/")[[1]][1])
        MyListY<-which(phenodata()$Condition == str_split(Comp,"/")[[1]][2])
        
        #Counting per conditions
        x<-apply(edata()[,MyListX],
                 1,
                 function(x){ length(x)-sum(is.na(x))}
        )
        y<-apply(edata()[,MyListY],
                 1,
                 function(x){ length(x)-sum(is.na(x))}
        )
        
        VV$MyCombin<-paste0(x,rep("_",length(x)),y)
        colnames(VV)[ncol(VV)]<-paste0("Nb Valid values ",str_split(Comp,"/")[[1]][1],"_",str_split(Comp,"/")[[1]][2])
        
      }
      
      #Add to fdata
      fdata<-merge(fdata,VV,by.x="RowNamesfdata",by.y="row.names(edata())",all.x=TRUE)
      
      #reorder columns fdata
      nColVV<-ncol(VV)-1
      startVV<-ncol(fdata)-nColVV+1
      
      Newfdata<-as.data.frame(fdata[,1:8])
      
      for(i in 1:nColVV){
        Newfdata<-cbind(Newfdata,fdata[,(9+((i-1)*3)):(11+((i-1)*3))])
        
        Newfdata<-cbind(Newfdata,fdata[,(startVV+i-1)])
        colnames(Newfdata)[ncol(Newfdata)]<-colnames(fdata)[(startVV+i-1)]
        
      }
      fdata<-Newfdata
      
      
      #Significant
      SeuilPQ<-1+if(MySignifType=="Pvalue"){0}else{1}
      
      for(i in 1:length(MyComparison)){
        
        Comp<-MyComparison[i]
        CompX<- str_split(Comp,"/")[[1]][1]
        CompY<- str_split(Comp,"/")[[1]][2]
        
        MaxX<-length(which(phenodata()$Condition == CompX))
        MaxY<-length(which(phenodata()$Condition == CompY))
        
        Signif<-function(x){
          #select columns of the comparison
          x<-x[(9+((i-1)*4)):(12+((i-1)*4))]
          
          MyRes<-""
          if(!is.na(x[SeuilPQ])){
            if(as.numeric(x[SeuilPQ])<MySignifVal && as.numeric(x[3])>log2(input$FCThreashold)){
              MyRes<-paste("Increased in",CompX)
            }else if(as.numeric(x[SeuilPQ])<MySignifVal && as.numeric(x[3])<(-log2(input$FCThreashold))){
              MyRes<-paste("Decreased in",CompX)
            }
            
          }
          
          VVX<-str_split(x[4],"_")[[1]][1]
          VVY<-str_split(x[4],"_")[[1]][2]
          
          if(VVX==0 && VVY==MaxY){
            MyRes<-paste("Disappeared in",CompX)
          }else if(VVX==MaxX && VVY==0){
            MyRes<-paste("Appeared in",CompX)
          }
          
          return(MyRes)
        }
        
        
        res<-apply(fdata,1,Signif)
        fdata$NewCol<-res
        colnames(fdata)[ncol(fdata)]<-paste("Significant",Comp)
      }
      
      #reorder columns of fdata
      nColS<-length(MyComparison)
      startS<-ncol(fdata)-nColS+1
      
      Newfdata<-as.data.frame(fdata[,1:8])
      
      for(i in 1:nColS){
        Newfdata<-cbind(Newfdata,fdata[,(startS+i-1)])
        colnames(Newfdata)[ncol(Newfdata)]<-colnames(fdata)[(startS+i-1)]
        
        Newfdata<-cbind(Newfdata,fdata[,(9+((i-1)*4)):(12+((i-1)*4))])
        
        
      }
      fdata<-Newfdata
      
      #ANOVA
      y<-edata()
      res<-as.data.frame(row.names(y))
      colnames(res)<-"Prot"

      for(i in 1:nrow(y)){

        x<-y[i,1:(ncol(y))]
        #filter
        x<-as.data.frame(t(x))
        x$Sample<-row.names(x)
        colnames(x)[1]<-c("Intensity")
        x<-merge(x,phenodata(),by.x="Sample",by.y="ShortName", all.x = TRUE)

        nbVV <- tapply(x, INDEX=x$Condition, FUN=complete.cases)#return TRUE FALSE Valid Values for each sample per group
        nbVV <- lapply(nbVV,sum)#nb VV per group

        PcVV<-nbVV

        for(condition in unique(x$Condition)) {
          PcVV[[condition]] <- PcVV[[condition]]/table(phenodata()$Condition)[[condition]]*100#%VV per condition
        }

        MyFilter <- any(PcVV >= 70) & sum(nbVV>=3)>=2#return TRUE there is at least 1 group with 70%VV and at least 3VV in an other one


        if(MyFilter==TRUE){

          #calcul anova
          x<-x[!is.na(x$Intensity),]
          x <- aov(Intensity ~ Condition, data = x)

          #extract Pvalue
          x <- summary(x)
          x<-x[[1]]$`Pr(>F)`[1]
          res$`Anova Pval`[i]<-x

        }else{res$`Anova Pval`[i]<-NA}
      }

      #Anova Qvalues
      res$`Anova Qval`<-p.adjust(p=res$`Anova Pval`, method = "BH")

      #Anova Significant
      if (input$ThreasholdType=="Pvalue"){n=0}else{n=1}
      res$`Anova significant`[res[,2+n]<input$SignifThreashold]<-"Anova +"
      res<-res[,c("Prot","Anova significant","Anova Pval","Anova Qval")]

      fdata<-merge(fdata,res,by.x="RowNamesfdata",by.y="Prot",all.x=TRUE)

      
      #Anova VV
      VVCalculator<-function(x){
        x<-cbind(Prot=1,x)
        x<-reshape2::melt(x,id.vars="Prot")
        x<-merge(x,phenodata(), by.x="variable",by.y="ShortName",all.x=TRUE)
        x<-tapply(x,INDEX=x$Condition,FUN=complete.cases)
        x<-lapply(x,sum)
        res<-paste(x,collapse = "_")
        y<-paste(names(x),collapse = "_")
        return(c(res,y))
      }
      
      VVdata<-c()
      for(i in 1:nrow(edata())){
        VVdata<-c(VVdata,VVCalculator(edata()[i,])[1])
      }  
      VVdata<- data.frame(id=row.names(edata()), Anova.Valid.Values=VVdata)
      fdata<-merge(fdata,VVdata, by.x="RowNamesfdata",by.y="id")
      colnames(fdata)[ncol(fdata)]<-paste("Anova Valid Values", VVCalculator(edata()[i,])[2])
      
      #Anova apparition/Disparition
      AppDispDetector<-function(x){
        x<-cbind(Prot=1,x)
        x<-reshape2::melt(x,id.vars="Prot")
        x<-merge(x,phenodata(), by.x="variable",by.y="ShortName",all.x=TRUE)
        x<-tapply(x,INDEX=x$Condition,FUN=complete.cases)
        VV<-unlist(lapply(x,sum))
        Tot<-unlist(lapply(x,length))
        x<-VV/Tot
        x<-any(x==1) & any(x==0)
        return(x)
      }
      
      AppDisp<-c()
      for(i in 1:nrow(edata())){
        AppDisp<-c(AppDisp,AppDispDetector(edata()[i,]))
      } 
      AppDisp<- data.frame(id=row.names(edata()), Anova.App.Disp=AppDisp)
      fdata<-merge(fdata,AppDisp, by.x="RowNamesfdata",by.y="id")
      fdata[fdata$Anova.App.Disp==TRUE,"Anova significant"]<-paste(fdata[fdata$Anova.App.Disp==TRUE,"Anova significant"],"Disappeared in at least one group",sep="/")
      fdata$`Anova significant`<-gsub(pattern="NA/", replacement="",x=fdata$`Anova significant`)
      
      fdata<-fdata[,-ncol(fdata)]
      
      
      
      fdata
    }
  }
  
  
  fdata<-eventReactive(input$Calculatefdata,{
    updateProgressBar(
      id = "pb1",
      value = 1, total = 3,
      title = "Process1/3"
    )
    generatefdata()
    
  })
  
  output$fdata<-renderDataTable({
    datatable(fdata(),
              options = list(pageLength=5, scrollX='400px'))
    })
  
  generateTableNbProtSignif<-function(){
    if(input$MyStat=="NO"){matrix()}else{
    
    fdata<-generatefdata()
    MyComparison<-input$MyComp
    MySignifFC<-log2(input$FCThreashold)
    NbComp<-length(MyComparison)
    
    x<-which(str_detect(string=colnames(fdata),pattern="Pvalue"))
    x<-unlist(lapply(X=x,FUN=function(x) x+0:2))
    MyTable<-fdata[,x]
    
    # MyTable to long
    MyTable$Prot<-row.names(MyTable)
    MyTable<-  reshape(MyTable,
                       idvar="Prot",
                       varying = colnames(MyTable)[1:(ncol(MyTable)-1)],
                       times= colnames(MyTable)[1:(ncol(MyTable)-1)],
                       direction = "long",
                       v.name="Value")
    
    #Add Type
    x<-str_split(MyTable$time," ")
    
    x<-unlist(x)
    MyTable$Type<-x[seq(from=1,to=length(x)-1,by=2)]
    MyTable$Comp<-x[seq(from=2,to=length(x),by=2)]
    MyTable<-MyTable %>% select(-'time')
    
    MyTable<-  reshape(MyTable,
                       idvar=c("Prot","Comp"),#Names of lines
                       timevar = "Type",#Names of columns
                       direction = "wide",
                       v.name="Value")#Values
    
    #FDR pval0.05 and 0.01
    xPval5pc<-MyTable[MyTable$Value.Pvalue<0.05,]
    xPval5pc<-tapply(xPval5pc$Value.Qvalue,INDEX=xPval5pc$Comp,max)
    xPval5pc<-round(xPval5pc*100,1)
    xPval5pc<-data.frame(Comp=names(xPval5pc),FDR=xPval5pc)
    
    xPval1pc<-MyTable[MyTable$Value.Pvalue<0.01,]
    xPval1pc<-tapply(xPval1pc$Value.Qvalue,INDEX=xPval1pc$Comp,max)
    xPval1pc<-round(xPval1pc*100,1)
    xPval1pc<-data.frame(Comp=names(xPval1pc),FDR=xPval1pc)
    
    xQval5pc<-MyTable[MyTable$Value.Qvalue<0.05,]
    xQval5pc<-tapply(xQval5pc$Value.Pvalue,INDEX=xQval5pc$Comp,max)
    xQval5pc<-round(xQval5pc,3)
    xQval5pc<-data.frame(Comp=names(xQval5pc),SeuilPval=xQval5pc)
    
    xQval1pc<-MyTable[MyTable$Value.Qvalue<0.01,]
    xQval1pc<-tapply(xQval1pc$Value.Pvalue,INDEX=xQval1pc$Comp,max)
    xQval1pc<-round(xQval1pc,3)
    xQval1pc<-data.frame(Comp=names(xQval1pc),SeuilPval=xQval1pc)
    
    
    #Add Signif
    MyTable$SignifPval5pcFC<-(MyTable[,"Value.Pvalue"]<0.05
                              & (MyTable[,"Value.Log2"]<(-log2(input$FCThreashold) )
                                 | MyTable[,"Value.Log2"]>log2(input$FCThreashold)))
    MyTable$SignifPval1pcFC<-(MyTable[,"Value.Pvalue"]<0.01
                              & (MyTable[,"Value.Log2"]<(-log2(input$FCThreashold) )
                                 | MyTable[,"Value.Log2"]>log2(input$FCThreashold)))
    MyTable$SignifQval5pcFC<-(MyTable[,"Value.Qvalue"]<0.05
                              & (MyTable[,"Value.Log2"]<(-log2(input$FCThreashold) )
                                 | MyTable[,"Value.Log2"]>log2(input$FCThreashold)))
    MyTable$SignifQval1pcFC<-(MyTable[,"Value.Qvalue"]<0.01
                              & (MyTable[,"Value.Log2"]<(-log2(input$FCThreashold) )
                                 | MyTable[,"Value.Log2"]>log2(input$FCThreashold)))
    
    MyTable$Pval<-!is.na(MyTable[,"Value.Pvalue"])
    
    #delete unused columns
    MyTable<-MyTable %>% select(-c('Value.Pvalue','Value.Qvalue','Value.Log2'))
    
    
    #Counting
    Counting<-matrix(MyComparison,nrow=length(MyComparison),ncol=1)
    
    #Pval<0.05
    x<-as.data.frame(table(MyTable$Comp,MyTable$SignifPval5pcFC))
    x<-x[x$Var2==TRUE,]
    x<-x[,c(1,3)]
    colnames(x)[2]<-"Pvalue<0.05"
    Counting<-merge(Counting,x,by.y="Var1",by.x="V1",all.x=TRUE)
    
    #Pval<0.01
    x<-as.data.frame(table(MyTable$Comp,MyTable$SignifPval1pcFC))
    x<-x[x$Var2==TRUE,]
    x<-x[,c(1,3)]
    colnames(x)[2]<-"Pvalue<0.01"
    Counting<-merge(Counting,x,by.y="Var1",by.x="V1",all.x=TRUE)
    
    #Qval<0.05
    x<-as.data.frame(table(MyTable$Comp,MyTable$SignifQval5pcFC))
    x<-x[x$Var2==TRUE,]
    x<-x[,c(1,3)]
    colnames(x)[2]<-"Qvalue<0.05"
    Counting<-merge(Counting,x,by.y="Var1",by.x="V1",all.x=TRUE)
    
    #Qval<0.01
    x<-as.data.frame(table(MyTable$Comp,MyTable$SignifQval1pcFC))
    x<-x[x$Var2==TRUE,]
    x<-x[,c(1,3)]
    colnames(x)[2]<-"Qvalue<0.01"
    Counting<-merge(Counting,x,by.y="Var1",by.x="V1",all.x=TRUE)
    
    #Nombre de protéines testées
    x<-as.data.frame(table(MyTable$Comp,MyTable$Pval))
    x<-x[x$Var2==TRUE,]
    x<-x[,c(1,3)]
    colnames(x)[2]<-"Nb tested proteins"
    Counting<-merge(Counting,x,by.y="Var1",by.x="V1",all.x=TRUE)
    
    Counting[is.na(Counting)]<-0
    colnames(Counting)[1]<-"Comparison"
    
    #Add FDR for seuil P-val
    for(Comp in MyComparison){
      i<-which(Counting$Comparison==Comp)
      FDR5pc<-xPval5pc[which(xPval5pc$Comp==Comp),"FDR"]
      Counting[i,2]<-paste0(Counting[i,2]," (FDR=",FDR5pc,"%)")
      
      FDR1pc<-xPval1pc[which(xPval1pc$Comp==Comp),"FDR"]
      Counting[i,3]<-paste0(Counting[i,3]," (FDR=",FDR1pc,"%)")
      
      Pval5pc<-xQval5pc[which(xQval5pc$Comp==Comp),"SeuilPval"]
      Counting[i,4]<-paste0(Counting[i,4]," (Pvalue=",Pval5pc,")")
      
      Pval1pc<-xQval1pc[which(xQval1pc$Comp==Comp),"SeuilPval"]
      Counting[i,5]<-paste0(Counting[i,5]," (Pvalue=",Pval1pc,")")
      
    }
    Counting
    }
  }
  
  MyTableSignif<-eventReactive(input$Calculatefdata,{
    updateProgressBar(
      id = "pb1",
      value = 2, total = 3,
      title = "Process2/3"
    )
    generateTableNbProtSignif()
  })
  
  output$TableNbProtSignif<-renderDataTable({
    datatable(MyTableSignif(),
              options = list(pageLength=10, scrollX='400px'))
    
    })

  
  generateTableApDis<-function(){
    if(input$MyStat=="NO"){matrix()}else{
    fdata<-generatefdata()
    MyComparison<-input$MyComp
    NbComp<-length(MyComparison)
    x<-which(str_detect(string=colnames(fdata),pattern="Significant"))
    MyTable<-as.data.frame(fdata[,x] )
    
    #Counting
    Counting<-data.frame(Comparison=MyComparison,Numerator.Appeared=rep(NA,NbComp),Numerator.Disappeared=rep(NA,NbComp))
    
    for( i in 1:NbComp){
      MyTitle<-MyComparison[i]
      NbAppeared<-length(MyTable[str_detect(MyTable[,i],"Appeared")==TRUE,i])
      NbDisappeared<-length(MyTable[str_detect(MyTable[,i],"Disappeared")==TRUE,i])
      
      Counting[i,2]<-NbAppeared
      Counting[i,3]<-NbDisappeared
    }
    
    Counting
    }
  }
  
  MyTableApDis<-eventReactive(input$Calculatefdata,{
    updateProgressBar(
      id = "pb1",
      value = 3, total = 3,
      title = "Process 3/3")
    
    generateTableApDis()
    
    
  })
  
  output$TableAppearedDsiappeared<-renderDataTable({
      datatable(MyTableApDis(),
              options = list(pageLength=10, scrollX='400px'))
    
    })
  
  #####Pepdata ----
  generatePepdata<-function(){
    if(MySoft()=="DIANN"){
      #import peptides data
      Pepdata<-fread(paste0(dirname(input$MyPath),"/report.pr_matrix.tsv"),sep="\t",header=TRUE)
      Pepdata<-as.data.frame(Pepdata)
      
      
      #extract sample name
      x<-colnames(Pepdata)[11:ncol(Pepdata)]
      x<-str_extract(string=x,pattern=".*\\\\(.*)\\.d$",group=1)
      x<-str_extract(string=x,pattern=input$MyRegEx,group=1)
      colnames(Pepdata)[11:ncol(Pepdata)]<-x
      
      #Prepare Pepdata
      Pepdata$Gene_ProteinGroup<-paste(Pepdata$Genes,Pepdata$Protein.Group)
      Pepdata<-Pepdata[,-c(1:9)]
      Pepdata<-Pepdata[,c(1,ncol(Pepdata),2:(ncol(Pepdata)-1))]
      colnames(Pepdata)[1]<-"Peptide"
      
    }else{
      #import peptides data
      Pepdata<-fread(paste0(dirname(input$MyPath),"/peptides.txt"),sep="\t",header=TRUE)
      Pepdata<-as.data.frame(Pepdata)
      x<-Pepdata[,c("Sequence","Protein group IDs")]
      Pepdata<-cbind(x,Pepdata[,which(str_detect(colnames(Pepdata),"Intensity "))])
      Pepdata<-separate_rows(Pepdata,'Protein group IDs',sep=";")
      x<-fread(MyFile(),select=c("Protein IDs","Gene names","id"))
      x$Gene_ProteinGroup<-paste(x$`Gene names`,x$`Protein IDs`)
      x<-x[,3:4]
      x$id<-as.character(x$id)
      Pepdata<-left_join(Pepdata, x, by = join_by('Protein group IDs' == id))
      Pepdata<-Pepdata[,-2]
      Pepdata<-Pepdata[,c(1,ncol(Pepdata),2:(ncol(Pepdata)-1))]
      colnames(Pepdata)[1]<-"Peptide"
      colnames(Pepdata)[3:ncol(Pepdata)]<-str_extract(colnames(Pepdata)[3:ncol(Pepdata)],"Intensity (.*)",group=1)
      Pepdata[Pepdata==0]<-NA
    }
    
    Pepdata[,3:ncol(Pepdata)]<-log2(Pepdata[,3:ncol(Pepdata)])
    Pepdata
  }
  
  Pepdata<-reactive({
    generatePepdata()
  })
  
  
  
  
  #####Descriptive plots  ----
  ######ProfilPlotConta----
  
  GraphConta<-reactive({
    df<-generatedf()
    phenodata<-generatephenodata()
    
    x<-df[,7:ncol(df)]
    x$Prot<-row.names(x)
    x<-reshape(x,
               idvar="Prot",
               varying = colnames(x)[1:(ncol(x)-1)],
               times= colnames(x)[1:(ncol(x)-1)], 
               direction = "long",
               v.name="Intensity")
    colnames(x)[2]<-"Sample"
    
    #annotation of contaminant proteins, BSA and Trypsin
    Myfdata<-data.frame(RowNamesfdata=row.names(df),
                        ProteinGroup=df$Protein.Group,
                        Genes=df$Genes,
                        Contaminant=df$Potential.contaminant)
    Myfdata$Type<-"Other"
    Myfdata[Myfdata$Contaminant=="+","Type"]<-"Contaminant"
    Myfdata[str_detect(Myfdata$ProteinGroup,"P00761")==TRUE,"Type"]<-"Trypsin" 
    Myfdata[str_detect(Myfdata$ProteinGroup,"P02769")==TRUE,"Type"]<-"BSA"
    
    x<-merge(x,Myfdata,by.x="Prot",by.y="RowNamesfdata",all.x=TRUE)
    x<-x[order(x$Type,decreasing=FALSE),]
    
    if(any(str_detect(x$Type,"Trypsin"))==TRUE){
      MyTrypsin<-x[str_detect(x$Type,"Trypsin")==TRUE,]
      MyAllOther<-x[str_detect(x$Type,"Trypsin")==FALSE,]
      x<-rbind(MyTrypsin,MyAllOther)
    }
    
    
    
    level_orderX <- factor(x$Sample, level = phenodata$ShortName)
    level_orderProt<-factor(x$Prot, level = unique(x$Prot))
    
    OrderType<-rev(unique(x$Type))
    TypeProt<-factor(x$Type,level=OrderType)
    
    MyColors<-OrderType
    MyColors[MyColors=="BSA"]<-"#008000"
    MyColors[MyColors=="Contaminant"]<-"black"
    MyColors[MyColors=="Other"]<-"grey"
    MyColors[MyColors=="Trypsin"]<-"#FF0080"
    
    MySize<-OrderType
    MySize[MySize!="Other"]<-1
    MySize[MySize=="Other"]<-0.5
    MySize<-as.numeric(MySize)
    
    
    ggplot(x, aes(x =  level_orderX, y = Intensity, group = level_orderProt,size = TypeProt,color=TypeProt, label=Genes)) + 
      geom_line()+
      scale_size_manual(values = MySize) +
      scale_colour_manual(values = MyColors)+
      theme_light()+
      theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.25),
            axis.title.x = element_blank())+
      ggtitle("Profilplot of Intensity per proteins (Contaminants)") 
    
  })
  
  output$ProfilConta<-renderPlotly({ GraphConta() })
  
  output$downloadProfilConta <-Exportggplot(graph=GraphConta(),
                                            basename="ProfilPlotConta",
                                            format=input$ProfilContaformat,
                                            width=input$WidthProfilConta,
                                            height=input$HeightProfilConta)
  ######ProfilPlotTop5----
  
  GraphTop5<-reactive({
    df<-generatedf()
    phenodata<-generatephenodata()
    
    x<-df[df$Potential.contaminant!="+",7:ncol(df)]
    x$Prot<-row.names(x)
    x<-reshape(x,
               idvar="Prot",
               varying = colnames(x)[1:(ncol(x)-1)],
               times= colnames(x)[1:(ncol(x)-1)], 
               direction = "long",
               v.name="Intensity")
    colnames(x)[2]<-"Sample"
    
    #annotation of contaminant proteins, BSA and Trypsin
    Myfdata<-data.frame(RowNamesfdata=row.names(df),
                        ProteinGroup=df$Protein.Group,
                        Genes=df$Genes,
                        Contaminant=df$Potential.contaminant)
    Myfdata$Type<-"Other"
    Myfdata<-Myfdata[Myfdata$Contaminant!="+",]
    
    x<-merge(x,Myfdata,by.x="Prot",by.y="RowNamesfdata",all.x=TRUE)
    
    #Find Top5 of proteins
    MyMaxLFQ<-tapply(x$Intensity, x$Prot, function(x){max(x,na.rm=TRUE)})
    MyMaxLFQ<-as.data.frame(MyMaxLFQ)
    MyTOP5<-row.names(top_n(MyMaxLFQ,n=5))
    
    ExtractFirstAccession <- function(x){
      FirstAccession<-strsplit(x, split = ";") [[1]][1]
      return(FirstAccession)
    }
    
    MyAnnotation<-distinct(x[which(x$Prot==MyTOP5),c("Prot","ProteinGroup","Genes")])
    MyAnnotation$Accession <-as.character(lapply(MyAnnotation$ProteinGroup,ExtractFirstAccession))
    MyAnnotation$Gene<-as.character(lapply(MyAnnotation$Genes,ExtractFirstAccession))
    
    for(i in 1:5){
      x[x$Prot==MyTOP5[i],"Type"]<-paste(MyAnnotation$Accession[i],MyAnnotation$Gene[i])
    }
    
    #plot
    x<-x[order(x$Type),]
    
    MyOtherrows<-x[str_detect(x$Type,"Other")==TRUE,]
    MyTOP5rows<-x[str_detect(x$Type,"Other")==FALSE,]
    x<-rbind(MyTOP5rows,MyOtherrows)
    
    
    level_orderX <- factor(x$Sample, level = phenodata$ShortName)
    level_orderProt<-factor(x$Prot, level = unique(x$Prot))
    
    OrderType<-rev(unique(x$Type))
    TypeProt<-factor(x$Type,level=OrderType)
    
    MyColors<-OrderType
    for(i in 2:6){
      MyColors[MyColors==OrderType[i]]<-hue_pal()(5)[i-1]
    }
    MyColors[MyColors=="Other"]<-"grey"
    
    
    MySize<-OrderType
    MySize[MySize!="Other"]<-1
    MySize[MySize=="Other"]<-0.5
    MySize<-as.numeric(MySize)
    
    
    ggplot(x, aes(x =  level_orderX, y = Intensity, group = level_orderProt, size = TypeProt, color=TypeProt,label=Genes)) + 
      geom_line() +
      scale_size_manual(values = MySize) +
      scale_colour_manual(values = MyColors) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.25),
            axis.title.x = element_blank()) +
      ggtitle("Profilplot of LFQ intensity per proteins (TOP5 without Contaminants)") 
    
  })
  
  output$ProfilTop5<-renderPlotly({ GraphTop5() })
  
  output$downloadProfilTop5 <-Exportggplot(graph=GraphTop5(),
                                           basename="ProfilPlotTop5",
                                           format=input$ProfilTop5format,
                                           width=input$WidthProfilTop5,
                                           height=input$HeightProfilTop5)
  
  ######NbProt----
  
  GraphNbProt<-reactive({
    req(edata(),phenodata())
    
    NbProt <- colSums(!is.na(edata()))
    
    x <- data.frame(Sample=colnames(edata()),
                    NbProt=NbProt,
                    Condition=phenodata()$Condition,
                    CondColor=phenodata()$CondColor)
    
    Conditions<-factor(x$Condition,levels=unique(x$Condition))
    MyColor<-unique(x$CondColor)
    MySample<-factor(x$Sample,levels=x$Sample)
    
    ggplot(x, aes(x=MySample, y=NbProt, color=Conditions,fill=Conditions)) +
      geom_bar(stat="identity")+
      ylim(0,NA)+
      scale_colour_manual(values = MyColor)+
      scale_fill_manual(values = MyColor)+
      theme_light()+
      theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.25),
            axis.title.x = element_blank())+
      ggtitle("Number of Protein per sample")
    
  }) 
  
  output$nbprot<-renderPlotly({ GraphNbProt() })
  
  output$downloadNbProtplot <-Exportggplot(graph=GraphNbProt(),
                                           basename="NbProt",
                                           format=input$NbProtformat,
                                           width=input$WidthNbProt,
                                           height=input$HeightNbProt)
  
  
  ######PCA----
  
  DataToShow <-eventReactive(input$DrawPCA,{
    req(edata(),phenodata())
    
    FilterGroup <- function(data, MinPc,type=c("1Group","EachGroup","Total")) {
      ValidProt<-rep(NA,nrow(data))
      
      
      if(type=="Total"){
        ValidProt<-100*rowSums(!is.na(edata()))/ncol(edata())
        ValidProt<-ValidProt>=MinPc
      }else{
        for(i in 1:nrow(data)){
          nbVV <- tapply(as.numeric(data[i,]),
                         INDEX=phenodata()[,"Condition"],
                         FUN=complete.cases)#show TRUE FALSE Valid Values for each sample per group
          nbVV <- lapply(nbVV,sum)#count number of Valid Values per group
          
          for(condition in levels(factor(phenodata()[,"Condition"]))) {
            nbVV[[condition]] <- nbVV[[condition]]/table(phenodata()[,"Condition"])[[condition]]*100#%Valid Values per condition.
          }
          
          if (type=="1Group"){
            ValidProt[i] <- any(nbVV >= MinPc)#return TRUE if there is at least one group with 70% of valid values.
          }else if(type=="EachGroup"){
            ValidProt[i] <- !any(nbVV < MinPc)
          }
        }
      }
      return(ValidProt)
    }
    
    if(input$FilterType=="% of Valid value in at least one group"){
      MyType<-"1Group"
    } else if (input$FilterType=="% of Valid value in each group"){
      MyType<-"EachGroup"
    } else if (input$FilterType=="% of Valid value in total"){
      MyType<-"Total"
    }
    
    edata()[FilterGroup(edata(), MinPc = input$PcFilter,type=MyType),]
    
  })
  
  output$NbProtPCA<-renderText({
    paste(as.character(nrow(DataToShow())),"/",as.character(nrow(edata())),"proteins used for PCA.")
  })
  
  my.prc <-eventReactive(input$DrawPCA,{
    req(edata(),phenodata(),DataToShow())
    
    DataToShow <- impute.MinProb(DataToShow(),q = 0.01,tune.sigma =1)
    
    if(input$TypeDataPCA=="Z-score(intensity)"){
      x<-colnames(DataToShow)
      DataToShow <-t(apply(DataToShow,1,scale))
      colnames(DataToShow)<-x  
      prcomp(t(DataToShow) , center=F, scale=F)
    }else if (input$TypeDataPCA=="Log2(intensity)"){
      prcomp(t(DataToShow) , center=T, scale=F)
    }
    
    
  })
  
  PCATitle <-eventReactive(input$DrawPCA,{
    req(edata(),phenodata(),DataToShow())
    paste("PCA", input$TypeDataPCA ,"of proteins with", input$PcFilter,input$FilterType)
  })
  
  GraphPCA<-reactive({
    req(edata(),phenodata(),MyCond())
    
    #PC1 versus PC2
    var_explained <- round(my.prc()$sdev^2/sum(my.prc()$sdev^2)*100,1)
    
    Conditions<-as.factor(phenodata()$Condition)
    MyColors<-levels(Conditions)
    MyColors<-match(MyColors,phenodata()$Condition)
    MyColors<-phenodata()$CondColor[MyColors]
    
    x<-as.data.frame(my.prc()$x ) 
    x$Sample<-rownames(x)
    
    ggplot(x,aes(x=PC1,y=PC2,color=Conditions)) +
      geom_point(size=4) +
      scale_colour_manual(values = MyColors)+
      scale_fill_manual(values = MyColors)+
      geom_text_repel(data = x,
                      aes(label = Sample),
                      nudge_y = 2,size=5,
                      segment.size  = 0.5,
                      force=1,
                      max.overlaps = Inf)+
      geom_encircle(
        aes(group = Conditions,
            fill = Conditions),
        alpha = 0.1,
        size = 2,
        show.legend = FALSE,
        na.rm = TRUE)+
      theme_bw(base_size=17) + 
      labs(x=paste0("PC1: ",var_explained[1],"%"),
           y=paste0("PC2: ",var_explained[2],"%")) +
      theme(legend.position="right")+
      ggtitle(PCATitle())
  }) 
  
  output$PCA<-renderPlot({ GraphPCA() })
  
  output$downloadPCAplot <-Exportggplot(graph=GraphPCA(),
                                        basename="PCA",
                                        format=input$PCAformat,
                                        width=input$WidthPCA,
                                        height=input$HeightPCA)
  
  
  GraphPCAprot<-reactive({
    req(my.prc(), fdata())
    
    x<-as.data.frame(my.prc()$rotation ) 
    x$ProtNumber<-rownames(x)
    y<-fdata()[,c(1,2,5,8)]
    x<-merge(x,y,by.x="ProtNumber",by.y="RowNamesfdata",all.x=TRUE, all.y=FALSE)
    
    ggplot(x,aes(x=PC1,y=PC2,label=Gene_ProteinGroup)) +
      geom_point(color="grey")+
      geom_text_repel(data = x[match(input$GenePCA,x$Gene_ProteinGroup),],
                      aes(label = Genes),
                      segment.size  = 0.5,
                      force=input$ForcePCA,
                      color = "black",
                      max.overlaps = Inf)+
      theme_bw()
  }) 
  
  output$PCAprotInteractive<-renderPlotly({ GraphPCAprot() })
  
  output$PCAprot<-renderPlot({ GraphPCAprot() })
  
  output$downloadPCAprotplot <-Exportggplot(graph=GraphPCAprot(),
                                            basename="PCAprot",
                                            format=input$PCAprotformat,
                                            width=input$WidthPCAprot,
                                            height=input$HeightPCAprot)
  
  ######Pearson correlations----
  ####### Heatmap ----
  GraphPearson<-reactive({
    req(edata(),phenodata(),MyCond())
    
    MyConditions<-data.frame(Conditions=phenodata()$Condition,row.names=colnames(edata()))
    
    MycolorCond<-phenodata()[,c(2,3)]
    MycolorCond<-distinct(MycolorCond)
    x<-MycolorCond[,2]
    names(x) <- MycolorCond[,1]
    MycolorCond<-list(x)
    names(MycolorCond)<-"Conditions"
    
    matDist <- cor(edata(),use = "pairwise.complete.obs")
    matDist[matDist == 1]<-NA 
    
    hmcol<- c(colorRampPalette(c("red","sienna1"))(250),
              colorRampPalette(c("sienna1","white"))(150),
              colorRampPalette(brewer.pal(9, 'GnBu'))(100))
    
    MyBreaks<-c(seq(from=0,to=0.4999,length.out=250),
                seq(from=0.5,to=0.7999,length.out=150),
                seq(from=0.8,to=1,length.out=100))
    
    ComplexHeatmap::pheatmap(matDist,  
             color=hmcol,
             breaks=MyBreaks,
             border_color = NA,
             cluster_rows=input$PearsonCluster,
             cluster_cols = input$PearsonCluster,
             legend_breaks=c(0, 0.5, 0.8, 0.9, 1),
             annotation_col = MyConditions,
             annotation_row = MyConditions,
             annotation_colors=MycolorCond,
             annotation_names_col = FALSE,
             annotation_names_row =  FALSE,
             show_rownames=F,
             show_colnames=T,
             main = "Pearson's correlation between log2(Intensity)of proteins",
             na_col = "#000000",
             # fontsize = 13,
             # fontsize_col = 13,
             angle_col = "90"
    )
  }) 
  
  
  output$Pearson<-renderPlot({ 
    GraphPearson() 
  })
  
  output$downloadPearsonplot <-Exportheatmap(graph=GraphPearson(),
                                             basename="PearsonHeatmap",
                                             format=input$Pearsonformat,
                                             width=input$WidthPearson,
                                             height=input$HeightPearson)
  
  ####### Boxplot ----
  GraphPearsonBox<-reactive({
    req(edata(),phenodata(),MyCond())
    matDist <- cor(edata(),use = "pairwise.complete.obs")
    matDist[matDist == 1]<-NA
    matDist<-as.data.frame(matDist)
    matDist$Sample1<-row.names(matDist)
    matDist<-reshape(matDist,
                     idvar="Sample1",
                     varying = colnames(matDist)[1:(ncol(matDist)-1)],
                     times= colnames(matDist)[1:(ncol(matDist)-1)], 
                     direction = "long",
                     v.name="Pearson")
    matDist<-matDist[is.na(matDist$Pearson)==FALSE,]
    colnames(matDist)[2]<-"Sample2"
    matDist<-merge(x=matDist,
                   y=phenodata(),
                   by.x="Sample1",
                   by.y="ShortName",
                   all.x=TRUE,
                   all.y=FALSE)
    colnames(matDist)[4:5]<-paste0(colnames(matDist)[4:5],"1")
    
    matDist<-merge(x=matDist,
                   y=phenodata(),
                   by.x="Sample2",
                   by.y="ShortName",
                   all.x=TRUE,
                   all.y=FALSE)
    colnames(matDist)[6:7]<-paste0(colnames(matDist)[6:7],"2")
    
    for(i in 1:nrow(matDist)){
      matDist$Cond1_Cond2[i]<-paste(sort(c(matDist$Condition1[i],matDist$Condition2[i])),collapse = "/")
    }
    MyColor<-rep("lightgray",length(levels(factor(matDist$Cond1_Cond2))))
    
    MyNames<-boxplot(data=matDist,Pearson~Cond1_Cond2)$names
    MyNames<-str_split(MyNames,pattern="/")
    for(i in 1:length(MyNames)){
      if(MyNames[[i]][1]==MyNames[[i]][2]){
        #Find color condition
        x<-matDist[matDist$Condition1==MyNames[[i]][1],]$CondColor1
        MyColor[i]<-x
      }
    }
    
    matDist$Cond1_Cond2 <- as.factor(matDist$Cond1_Cond2)
    
    ggplot(matDist, aes(x=Cond1_Cond2, y=Pearson,fill=Cond1_Cond2)) + 
      geom_boxplot(outlier.shape = NA)+
      scale_fill_manual(values=MyColor)+
      geom_jitter(shape=16, position=position_jitter(0.2),size=0.5)+
      theme_light()+
      theme(axis.text.x = element_text(angle = 90,hjust=1),
            axis.title.x = element_blank())+
      theme(legend.position="none")+
      ggtitle("Boxplot of Pearson's correlation between samples") 
    
  }) 
  
  output$PearsonBox<-renderPlot({ GraphPearsonBox() })
  
  output$downloadPearsonBoxplot <-Exportggplot(graph=GraphPearsonBox(),
                                               basename="PearsonBox",
                                               format=input$PearsonBoxformat,
                                               width=input$WidthPearsonBox,
                                               height=input$HeightPearsonBox)
  
  
  ######Euclidean distances ----
  #######Heatmap ----
  GraphEuclid<-reactive({
    req(edata(),phenodata(),MyCond())
    
    MyConditions<-data.frame(Conditions=phenodata()$Condition,row.names=colnames(edata()))
    
    MycolorCond<-phenodata()[,c(2,3)]
    MycolorCond<-distinct(MycolorCond)
    x<-MycolorCond[,2]
    names(x) <- MycolorCond[,1]
    MycolorCond<-list(x)
    names(MycolorCond)<-"Conditions"
    
    matDist <- as.matrix(dist(t(edata())))
    matDist[matDist == 0]<-NA 
    
    
    
    hmcol<- rev(colorRampPalette(brewer.pal(9, 'GnBu'))(100))
    
    if(ncol(edata())<3){
      NULL
    }else{
      ComplexHeatmap::pheatmap(matDist,  
               color=hmcol,
               border_color = NA,
               cluster_rows=input$EuclidCluster,
               cluster_cols = input$EuclidCluster,
               annotation_col = MyConditions,
               annotation_row = MyConditions,
               annotation_colors=MycolorCond,
               annotation_names_col = FALSE,
               annotation_names_row =  FALSE,
               show_rownames=F,
               show_colnames=T,
               main = "Euclidean distance between log2(Intensity)of proteins",
               na_col = "#000000",
               # fontsize = 13,
               # fontsize_col = 13,
               angle_col = "90"
      )
    }
  }) 
  
  
  output$Euclid<-renderPlot({ GraphEuclid() })
  
  output$downloadEuclidplot <-Exportheatmap(graph=GraphEuclid(),
                                            basename="EuclideanHeatmap",
                                            format=input$Euclidformat,
                                            width=input$WidthEuclid,
                                            height=input$HeightEuclid)
  
  #######BoxPlot ----
  GraphEuclidBox<-reactive({
    req(edata(),phenodata(),MyCond())
    matDist <- as.matrix(dist(t(edata())))
    matDist[matDist == 0]<-NA 
    matDist<-as.data.frame(matDist)
    matDist$Sample1<-row.names(matDist)
    matDist<-reshape(matDist,
                     idvar="Sample1",
                     varying = colnames(matDist)[1:(ncol(matDist)-1)],
                     times= colnames(matDist)[1:(ncol(matDist)-1)], 
                     direction = "long",
                     v.name="Euclidean.distance")
    matDist<-matDist[is.na(matDist$Euclidean.distance)==FALSE,]
    colnames(matDist)[2]<-"Sample2"
    matDist<-merge(x=matDist,
                   y=phenodata(),
                   by.x="Sample1",
                   by.y="ShortName",
                   all.x=TRUE,
                   all.y=FALSE)
    colnames(matDist)[4:5]<-paste0(colnames(matDist)[4:5],"1")
    
    matDist<-merge(x=matDist,
                   y=phenodata(),
                   by.x="Sample2",
                   by.y="ShortName",
                   all.x=TRUE,
                   all.y=FALSE)
    colnames(matDist)[6:7]<-paste0(colnames(matDist)[6:7],"2")
    
    for(i in 1:nrow(matDist)){
      matDist$Cond1_Cond2[i]<-paste(sort(c(matDist$Condition1[i],matDist$Condition2[i])),collapse = "/")
    }
    
    MyColor<-rep("lightgray",length(levels(factor(matDist$Cond1_Cond2))))
    
    MyNames<-boxplot(data=matDist,Euclidean.distance~Cond1_Cond2)$names
    MyNames<-str_split(MyNames,pattern="/")
    for(i in 1:length(MyNames)){
      if(MyNames[[i]][1]==MyNames[[i]][2]){
        #Find color of the condition
        MyColor[i]<-matDist[matDist$Condition1==MyNames[[i]][1],]$CondColor1
      }
    }
    
    matDist$Cond1_Cond2 <- as.factor(matDist$Cond1_Cond2)
    
    ggplot(matDist, aes(x=Cond1_Cond2, y=Euclidean.distance,fill=Cond1_Cond2)) + 
      geom_boxplot(outlier.shape = NA)+
      scale_fill_manual(values=MyColor)+
      geom_jitter(shape=16, position=position_jitter(0.2),size=0.5)+
      theme_light()+
      theme(axis.text.x = element_text(angle = 90,hjust=1),
            axis.title.x = element_blank())+
      theme(legend.position="none")+
      ggtitle("Boxplot of Euclidean distance between samples") 
    
    
  }) 
  
  output$EuclidBox<-renderPlot({ 
    GraphEuclidBox() 
  })
  
  output$downloadEuclidBoxplot <-Exportggplot(graph=GraphEuclidBox(),
                                              basename="EuclideanBox",
                                              format=input$EuclidBoxformat,
                                              width=input$WidthEuclidBox,
                                              height=input$HeightEuclidBox)
  
  ######Heatmap Log2 ----
  GraphLog2<-reactive({
    req(edata(),phenodata(),MyCond())
    
    # Select proteins with at least 1 valid value in total
    x <- edata()[rowSums(is.na(edata())) != ncol(edata()), ]
    
    # Mean intensity
    x$MeanLFQ<-apply(x,1,function(x){mean(x,na.rm = TRUE)})
    
    # Sort proteins by mean intensity
    x<-x[order(x$MeanLFQ,decreasing = TRUE),]
    
    # Delete mean columns
    x<-x %>% select(-'MeanLFQ')
    
    #select 2000 prot if nbprot>2000
    nbprot<-nrow(x)
    if (nbprot>2000){
      step<-ceiling(nbprot/2000)
      x<-x[seq(from=1, to=nbprot,by=step ),]
    }
    
    hmcol<- colorRampPalette(c(input$LowColor,input$MiddleColor,input$HighColor))(100)
    
    MyConditions<-data.frame(Conditions=phenodata()$Condition,row.names=colnames(edata()))
    
    MycolorCond<-phenodata()[,c(2,3)]
    MycolorCond<-distinct(MycolorCond)
    y<-MycolorCond[,2]
    names(y) <- MycolorCond[,1]
    MycolorCond<-list(y)
    names(MycolorCond)<-"Conditions"
    
    
    #heatmap with sample clustering
    ComplexHeatmap::pheatmap(x,
             cluster_rows=FALSE,
             cluster_cols = input$SampleCluster,
             clustering_distance_cols = "euclidean",
             clustering_method = "average",
             col=hmcol,
             annotation_col = MyConditions,
             annotation_colors = MycolorCond,
             na_col = "grey",
             show_rownames=F,
             show_colnames=T,
             annotation_names_col = FALSE,
             annotation_names_row =  FALSE,
             main = "Heatmap of Log2(Intensity)",
             fontsize = 8,
             fontsize_col = 8,
             border_color = NA,
             angle_col="90")
    
  }) 
  
  
  output$HeatmapLog2<-renderPlot({ GraphLog2() })
  
  output$downloadHeatmapLog2plot <-Exportheatmap(graph=GraphLog2(),
                                                 basename="HeatmapLog2",
                                                 format=input$HeatmapLog2format,
                                                 width=input$WidthHeatmapLog2,
                                                 height=input$HeightHeatmapLog2)
  
  
  ######Scaterplot ----
  
  output$CorScaterplot<-renderText({
    x<-edata()[,which(colnames(edata())==input$Sample1)]
    y<-edata()[,which(colnames(edata())==input$Sample2)]
    paste("r =",round(cor(x,y,use="pairwise.complete.obs"),5))
  })
  
  
  DensityScater <-eventReactive(input$DrawScaterPlot,{
    req(edata(),fdata())
    
    x<-edata()[,which(colnames(edata())==input$Sample1)]
    y<-edata()[,which(colnames(edata())==input$Sample2)]
    
    DataScaterPlot(x,y,edata(),fdata())
    
  })
  
  DrawScater<-eventReactive(input$DrawScaterPlot,{
    req(DensityScater())
    DrawScaterplot(DensityScater(),TitleX=input$Sample1,TitleY=input$Sample2)
    
  })
  
  ScaterPlot<-reactive({
    req(DrawScater())
    
    g<-DrawScater()
    g+geom_text_repel(data = DensityScater()[match(input$GenesScaterplot,DensityScater()$Gene_ProteinGroup),],
                      aes(label = Genes),
                      segment.size  = 0.5,
                      force=input$ForceScaterplot,
                      color = "white",
                      max.overlaps = Inf)
    
    
  })
  
  output$Scaterplotinteractive<-renderPlotly({ ScaterPlot() })
  
  output$Scaterplot<-renderPlot({ ScaterPlot() })
  
  output$downloadScaterPlot <-Exportggplot(graph=ScaterPlot(),
                                           basename=paste0("ScaterPlot_",input$Sample1,"_",input$Sample2),
                                           format=input$ScaterPlotformat,
                                           width=input$WidthScaterPlot,
                                           height=input$HeightScaterPlot)
  
  ######MeanScaterplot ----
  
  output$CorMeanScaterplot<-renderText({
    x<-edata()[,which(phenodata()$Condition == input$Cond1)]
    x<-apply(x,MARGIN=1,mean, na.rm = TRUE)
    y<-edata()[,which(phenodata()$Condition == input$Cond2)]
    y<-apply(y,MARGIN=1,mean, na.rm = TRUE)
    paste("r =",round(cor(x,y,use="pairwise.complete.obs"),5))
  })
  
  
  DensityMeanScater  <-eventReactive(input$DrawMeanScaterPlot,{
    req(edata(),fdata())
    
    x<-edata()[,which(phenodata()$Condition == input$Cond1)]
    x<-apply(x,MARGIN=1,mean, na.rm = TRUE)
    y<-edata()[,which(phenodata()$Condition == input$Cond2)]
    y<-apply(y,MARGIN=1,mean, na.rm = TRUE)
    
    DataScaterPlot(x,y,edata(),fdata())
    
  })
  
  DrawMeanScater<-eventReactive(input$DrawMeanScaterPlot,{
    req(DensityMeanScater())
    DrawScaterplot(DensityMeanScater(),TitleX=input$Cond1,TitleY=input$Cond2)
    
  })
  
  MeanScaterPlot<-reactive({
    req(DrawMeanScater())
    
    g<-DrawMeanScater()
    g+geom_text_repel(data = DensityMeanScater()[match(input$GenesMeanScaterplot,DensityMeanScater()$Gene_ProteinGroup),],
                      aes(label = Genes),
                      segment.size  = 0.5,
                      force=input$ForceMeanScaterplot,
                      color = "white",
                      max.overlaps = Inf)
    
    
  })
  
  output$MeanScaterplotinteractive<-renderPlotly({ MeanScaterPlot() })
  
  output$MeanScaterplot<-renderPlot({ MeanScaterPlot() })
  
  output$downloadMeanScaterPlot <-Exportggplot(graph=MeanScaterPlot(),
                                               basename=paste0("MeanScaterPlot_",input$Cond1,"_",input$Cond2),
                                               format=input$MeanScaterPlotformat,
                                               width=input$WidthMeanScaterPlot,
                                               height=input$HeightMeanScaterPlot)
  
  #####Statistical plots  ----
  
  output$NbProtSignif<-renderText({
    req(MyListeProt(),MyColNumber())
    
    Pval<-max(fdata()[which(MyListeProt()==TRUE),MyColNumber()])
    Qval<-max(fdata()[which(MyListeProt()==TRUE),(MyColNumber()+1)])
    paste(as.character(length(which(MyListeProt()==TRUE))),"/",as.character(length(MyListeProt())),"significant proteins.","\nPvalue=",Pval,",FDR=",Qval,".")
  })
  
  ######Pvalplot ----
  GraphPval<-reactive({
    MyComparison<-input$MyComp
    nbComp<-length(MyComparison)
    
    MyColNumber<-which(str_detect(string=colnames(fdata()),pattern="Pvalue"))
    DataToShow<-fdata()[,c(2,MyColNumber)]
    colnames(DataToShow)<-c("Prot",MyComparison)
    DataToShow<-  reshape(DataToShow,
                          idvar="Prot",
                          varying = colnames(DataToShow)[2:ncol(DataToShow)],
                          times= colnames(DataToShow)[2:ncol(DataToShow)],
                          direction = "long",
                          v.name="Value")
    DataToShow<-DataToShow[is.na(DataToShow$Value)==FALSE,]
    colnames(DataToShow)[2:3]<-c("Comp", "Pvalue")
    
    
    Mycondition<-levels(factor(DataToShow$Comp))
    
    MyColor<-colorRampPalette(brewer.pal(8, "Set2"))(if(nbComp>8){nbComp}else{8})
    MyColor<-MyColor[1:nbComp]
    
    ggplot(DataToShow,aes(x=Pvalue,color=Comp)) +
      geom_histogram(aes(y = after_stat(density)),alpha=0,position = "identity")+
      geom_density(alpha=0, linewidth=1.1)+
      scale_color_manual(values=MyColor)+
      guides(color = guide_legend(title = "Comparison"))+
      theme_light()
    
  })
  
  output$PvalPlot<-renderPlot({ GraphPval() })
  
  output$downloadPvalPlot <-Exportggplot(graph=GraphPval(),
                                         basename="Pvalues",
                                         format=input$PvalPlotformat,
                                         width=input$WidthPvalPlot,
                                         height=input$HeightPvalPlot)
  
  ######FCplot ----
  GraphFC<-reactive({
    MyComparison<-input$MyComp
    nbComp<-length(MyComparison)
    
    MyColNumber<-which(str_detect(string=colnames(fdata()),pattern="Log2"))
    DataToShow<-fdata()[,c(2,MyColNumber)]
    colnames(DataToShow)<-c("Prot",MyComparison)
    DataToShow<-  reshape(DataToShow,
                          idvar="Prot",
                          varying = colnames(DataToShow)[2:ncol(DataToShow)],
                          times= colnames(DataToShow)[2:ncol(DataToShow)],
                          direction = "long",
                          v.name="Value")
    DataToShow<-DataToShow[is.na(DataToShow$Value)==FALSE,]
    colnames(DataToShow)[2:3]<-c("Comp", "Log2FC")
    
    Mycondition<-levels(factor(DataToShow$Comp))
    
    MyColor<-colorRampPalette(brewer.pal(8, "Set2"))(if(nbComp>8){nbComp}else{8})
    MyColor<-MyColor[1:nbComp]
    
    mu<-tapply(X=DataToShow$Log2FC,INDEX=DataToShow$Comp,FUN=median)
    mu2<-tapply(X=DataToShow$Log2FC,INDEX=DataToShow$Comp,FUN=sd)
    mu<-data.frame(Comp=names(mu),Mediane=mu, Ecart.Type=mu2)
    
    MyAnnotation<-c()
    for(i in 1:nbComp){
      MyAnnotation<-c(MyAnnotation,mu$Comp[i],
                      paste(" Median = ",round(mu$Mediane[i],2)),
                      paste(" sd = ",round(mu$Ecart.Type[i],2)),
                      " ")
    }
    
    MyTable<-fdata()[,MyColNumber] 
    MyTable<-as.data.frame(MyTable)
    MyPositionAnnotation<-min(MyTable,na.rm = TRUE)  
    
    d<-c()
    for(i in 1:ncol(MyTable)){
      x<-MyTable[!is.na(MyTable[,i]),i]
      d <- max(d,max(density(x)$y))
    }
    
    ggplot(DataToShow,aes(x=Log2FC,color=Comp)) +
      geom_density(alpha=0, linewidth=1.1)+
      scale_color_manual(values=MyColor)+
      geom_vline(data=mu, 
                 aes(xintercept=Mediane, color=Comp),
                 linetype="dashed")+
      annotate("text",
               x=rep(MyPositionAnnotation,4*nbComp),
               y=seq(from= d, to= d-(d*0.05*4*nbComp), length.out= 4*nbComp),
               size=3.5,
               col = rep(MyColor,each=4),
               label = MyAnnotation,
               hjust = 0)+
      labs(y = "Density")+
      guides(color = guide_legend(title = "Comparison"))+
      theme_light()
  })
  
  output$FCplot<-renderPlot({GraphFC()})
  output$downloadFCPlot <-Exportggplot(graph=GraphFC(),
                                       basename="FCplot",
                                       format=input$FCPlotformat,
                                       width=input$WidthFCPlot,
                                       height=input$HeightFCPlot)
  
  ######Volcanoplot ----
  
  GraphVolcano<-reactive({
    req(fdata())
    
    n<-if(input$ThreasholdType=="Pvalue"){0}else{1}
    
    x<-fdata()[is.na(fdata()[,MyColNumber()])==FALSE,]
    x<-x[,c(MyColNumber(),MyColNumber()+n,MyColNumber()+2,8,5)]
    x$Type<-input$NSColor
    x[x[,2]<input$SignifThreashold & x[,3]>log2(input$FCThreashold),"Type"]<-input$UpColor
    x[x[,2]<input$SignifThreashold & x[,3]<(-log2(input$FCThreashold)),"Type"]<-input$DownColor
    x[,1]<--log10(x[,1])
    colnames(x)<-gsub(pattern="\\s",replacement="_",x=colnames(x))
    colnames(x)<-gsub(pattern="\\/",replacement="_",x=colnames(x))
    
    g<-ggplot(x,aes_string(x=colnames(x)[3],y=colnames(x)[1],color="Type",label="Gene_ProteinGroup")) +
      geom_point(size=1) +
      scale_colour_manual(values = levels(as.factor(x$Type)))+
      geom_point(data=x[match(input$GenePoints,x$Gene_ProteinGroup),],
                 color=input$PointColor,
                 size=input$PointSize)+
      geom_text_repel(data = x[match(input$sel_gene_nm,x$Gene_ProteinGroup),],
                      aes_string(label = "Genes"),
                      segment.size  = 0.5,
                      force=input$Force,
                      color = "black",
                      max.overlaps = Inf)+
      theme_bw(base_size=10) + 
      labs(x=paste0("Log2(",input$ComparisonSelect,")"),
           y="-log10(Pvalue)") +
      theme(legend.position="none")+
      ggtitle(paste("Volcano plot",input$ComparisonSelect))
    
    if(input$ShowHline==TRUE){
      
      MyThreashold<-min(x[x[,2]<input$SignifThreashold,1])
      
      g<-g+geom_hline(linetype=input$LineType,
                      color=input$LineColor,
                      yintercept=MyThreashold
      )
    }
    
    if(input$ShowVline==TRUE){
      g<-g+geom_vline(linetype=input$LineType,
                      color=input$LineColor,
                      xintercept=log2(input$FCThreashold)
      )+
        geom_vline(linetype=input$LineType,
                   color=input$LineColor,
                   xintercept=-log2(input$FCThreashold)
        )
    }
    g
  })
  
  
  
  output$VolcanoPlot<-renderPlot({ GraphVolcano() })
  
  output$IntercativeVolcanoPlot<-renderPlotly({ GraphVolcano() })
  
  output$downloadVolcanoplot <-Exportggplot(graph=GraphVolcano(),
                                            basename=paste0("Volcanoplot_",input$ComparisonSelect),
                                            format=input$Volcanoformat,
                                            width=input$WidthVolcano,
                                            height=input$HeightVolcano)
  
  
  
  ######Heatmap significant ----
  MyColNumber <-reactive({
    req(fdata())
    MyColNumber<-which(colnames(fdata())==paste("Pvalue",input$ComparisonSelect))
    MyColNumber
  })
  
  MyListeProt <-reactive({
    req(fdata(),MyColNumber())
    n<-if(input$ThreasholdType=="Pvalue"){0}else{1}
    
    #filter significative proteins
    MyListeProt<-(fdata()[,(MyColNumber()+n)]<input$SignifThreashold
                  & (fdata()[,(MyColNumber()+2)]<(-log2(input$FCThreashold ))
                     |fdata()[,(MyColNumber()+2)]>log2(input$FCThreashold)))
    MyListeProt
  })
  
  
  GraphHMSignif<-reactive({
    req(fdata(),phenodata(),MyListeProt(),edata())
    
    #select columns of the comparison
    MyListX<-which(phenodata()$Condition == str_split(input$ComparisonSelect,"/")[[1]][1])
    MyListY<-which(phenodata()$Condition == str_split(input$ComparisonSelect,"/")[[1]][2])
    
    #filter significative proteins
    Genes<-fdata()[which(MyListeProt()==TRUE),c(5,2)]#"Genes","Protein.Group"
    Genes[is.na(Genes$Genes),"Genes"]<-Genes[is.na(Genes$Genes),2]#"Protein.Group"
    Genes<-Genes$Genes
    
    MyListeProt<-fdata()[which(MyListeProt()==TRUE),"RowNamesfdata"]
    MyListeProt<-data.frame(x1=MyListeProt,x2=MyListeProt)
    
    DataToShow<-cbind(edata(),RowNamesedata=row.names(edata()))
    DataToShow<-merge(DataToShow,MyListeProt,by.x="RowNamesedata",by.y="x1",all=FALSE)
    DataToShow<-DataToShow %>% select(-c('RowNamesedata','x2'))
    
    DataToShow<-DataToShow[,c(MyListX,MyListY)]
    
    if(nrow(MyListeProt)>1){  
      #Z-score
      MyZscore<-t(apply(DataToShow,1,scale))
      colnames(MyZscore)<-colnames(DataToShow)
      rownames(MyZscore)<-Genes
      
      hmcol<- colorRampPalette(c(input$LowColorSignif,input$MiddleColorSignif,input$HighColorSignif))(100)
      hmcol<-c(input$LowColorSignif,hmcol,input$HighColorSignif)
      
      MyBreaks<-seq(from=-2,to=2,length.out=100)
      MinBreaks<-if(min(MyZscore,na.rm=TRUE)<(-2)){min(MyZscore,na.rm=TRUE)}else{-2.0001}
      MaxBreaks<-if(max(MyZscore,na.rm=TRUE)>2){max(MyZscore,na.rm=TRUE)}else{2.0001}
      MyBreaks<-c(MinBreaks,MyBreaks,MaxBreaks)
      
      MyConditions<-data.frame(Conditions=phenodata()$Condition[c(MyListX,MyListY)],row.names=colnames(MyZscore))
      
      MycolorCond<-phenodata()[c(MyListX[1],MyListY[1]),c(2,3)]
      x<-MycolorCond[,2]
      names(x) <- MycolorCond[,1]
      MycolorCond<-list(x)
      names(MycolorCond)<-"Conditions"
      
      cor_matrix <- cor(MyZscore, use = "pairwise.complete.obs", method = "pearson")
      dist_matrix <- as.dist(1 - cor_matrix)
      dist_matrix[is.na(dist_matrix)]<-2
      hc_cols <- hclust(dist_matrix,method = "complete")
      
      
      cor_matrix <- cor(t(MyZscore), use = "pairwise.complete.obs", method = "pearson")
      dist_matrix <- as.dist(1 - cor_matrix)
      dist_matrix[is.na(dist_matrix)]<-2
      hc_rows <- hclust(dist_matrix,method = "complete")
      
      ComplexHeatmap::pheatmap(MyZscore,
               cluster_rows = hc_rows,
               cluster_cols = hc_cols,
               color=hmcol,
               breaks=MyBreaks,
               legend_breaks=c(if(MinBreaks==-2.0001){-2}else{c(round(MinBreaks,1),-2)},
                               if(MaxBreaks==2.0001){2}else{c(2,round(MaxBreaks,1))}),
               annotation_col = MyConditions,
               annotation_colors=MycolorCond,
               na_col = "grey",
               show_rownames=input$ShowGeneName,
               show_colnames=TRUE,
               annotation_names_col = FALSE,
               annotation_names_row =  FALSE,
               main = paste0("Heatmap of Zscore of differential proteins ",
                             input$ComparisonSelect,
                             " (",input$ThreasholdType,input$SignifThreashold,if(input$FCThreashold!=1){paste0(", FC>",input$FCThreashold)}else{""},
                             ")"
               ),
               border_color = NA,
               angle_col="90"
      )
    }
  })
  
  output$HMSignif<-renderPlot({ GraphHMSignif() })
  
  output$downloadHMSignif <-Exportheatmap(graph=GraphHMSignif(),
                                          basename=paste0("HeatmapSignificant_",input$ComparisonSelect),
                                          format=input$HMSignifformat,
                                          width=input$WidthHMSignif,
                                          height=input$HeightHMSignif)
  
  ######Heatmap Anova ----
  MyColNumberAnova<-reactive({
    req(fdata())
    MyColNumber<-which(colnames(fdata())=="Anova Pval")
    MyColNumber
  })
  
  MyListeProtAnova <-reactive({
    req(fdata(),MyColNumberAnova())
    n<-if(input$ThreasholdType=="Pvalue"){0}else{1}
    
    #filter significative proteins
    MyListeProt<-fdata()[,(MyColNumberAnova()+n)]<input$SignifThreashold
    MyListeProt
  })
  
  output$NbProtSignifAnova<-renderText({
    req(MyListeProtAnova(),MyColNumberAnova())
    
    Pval<-max(fdata()[which(MyListeProtAnova()==TRUE),MyColNumberAnova()])
    Qval<-max(fdata()[which(MyListeProtAnova()==TRUE),(MyColNumberAnova()+1)])
    paste(as.character(length(which(MyListeProtAnova()==TRUE))),"/",as.character(length(MyListeProtAnova())),"significant proteins.","\nPvalue=",Pval,",FDR=",Qval,".")
  })
  
  GraphHMAnova<-reactive({
    req(fdata(),phenodata(),edata())

    #filter significative proteins
    Genes<-fdata()[which(MyListeProtAnova()==TRUE),c(5,2)]#"Genes","Protein.Group"
    Genes[is.na(Genes$Genes),"Genes"]<-Genes[is.na(Genes$Genes),2]#"Protein.Group"
    Genes<-Genes$Genes
    
    MyListeProt<-fdata()[which(MyListeProtAnova()==TRUE),"RowNamesfdata"]
    MyListeProt<-data.frame(x1=MyListeProt,x2=MyListeProt)
    
    DataToShow<-cbind(edata(),RowNamesedata=row.names(edata()))
    DataToShow<-merge(DataToShow,MyListeProt,by.x="RowNamesedata",by.y="x1",all=FALSE)
    DataToShow<-DataToShow %>% select(-c('RowNamesedata','x2'))

    
    if(nrow(MyListeProt)>1){  
      #Z-score
      MyZscore<-t(apply(DataToShow,1,scale))

      colnames(MyZscore)<-colnames(DataToShow)
      rownames(MyZscore)<-Genes
      
      hmcol<- colorRampPalette(c(input$LowColorSignifAnova,input$MiddleColorSignifAnova,input$HighColorSignifAnova))(100)
      hmcol<-c(input$LowColorSignifAnova,hmcol,input$HighColorSignifAnova)
      
      MyBreaks<-seq(from=-2,to=2,length.out=100)
      MinBreaks<-if(min(MyZscore,na.rm=TRUE)<(-2)){min(MyZscore,na.rm=TRUE)}else{-2.0001}
      MaxBreaks<-if(max(MyZscore,na.rm=TRUE)>2){max(MyZscore,na.rm=TRUE)}else{2.0001}
      MyBreaks<-c(MinBreaks,MyBreaks,MaxBreaks)
      
      MyConditions<-data.frame(Conditions=phenodata()$Condition,row.names=colnames(MyZscore))
      
      MycolorCond<-phenodata()[,c(2,3)]
      x<-MycolorCond[,2]
      names(x) <- MycolorCond[,1]
      MycolorCond<-list(x)
      names(MycolorCond)<-"Conditions"
      
      cor_matrix <- cor(MyZscore, use = "pairwise.complete.obs", method = "pearson")
      dist_matrix <- as.dist(1 - cor_matrix)
      dist_matrix[is.na(dist_matrix)]<-2
      hc_cols <- hclust(dist_matrix,method = "complete")
      
      
      cor_matrix <- cor(t(MyZscore), use = "pairwise.complete.obs", method = "pearson")
      dist_matrix <- as.dist(1 - cor_matrix)
      dist_matrix[is.na(dist_matrix)]<-2
      hc_rows <- hclust(dist_matrix,method = "complete")
      
      
      ComplexHeatmap::pheatmap(MyZscore,
                               cluster_rows = hc_rows,
                               cluster_cols = hc_cols,
                               color=hmcol,
                               breaks=MyBreaks,
                               legend_breaks=c(if(MinBreaks==-2.0001){-2}else{c(round(MinBreaks,1),-2)},
                                               if(MaxBreaks==2.0001){2}else{c(2,round(MaxBreaks,1))}),
                               annotation_col = MyConditions,
                               annotation_colors=MycolorCond,
                               na_col = "grey",
                               show_rownames=input$ShowGeneNameAnova,
                               show_colnames=TRUE,
                               annotation_names_col = FALSE,
                               annotation_names_row =  FALSE,
                               main = paste0("Heatmap of Zscore of Anova differential proteins ",
                                             " (",input$ThreasholdType,input$SignifThreashold,")"
                               ),
                               border_color = NA,
                               angle_col="90"
      )
    }
  })
  
  output$HMAnova<-renderPlot({ GraphHMAnova() })
  
  output$downloadHMAnova <-Exportheatmap(graph=GraphHMAnova(),
                                          basename="HeatmapAnova_",
                                          format=input$HMAnovaformat,
                                          width=input$WidthHMAnova,
                                          height=input$HeightHMAnova)
  
  
  
  
  #####Proteins plots  ----
  ###### ProtPlot ----
  df2<-reactive({
    Myedata<-2^edata()
    Myedata<-cbind(Myedata,RowNamesedata=row.names(edata()))
    # Myfdata<-fdata()[,c(1,ncol(fdata()))]
    Myfdata<-fdata()[,c("RowNamesfdata","Gene_ProteinGroup")]
    df<-merge(Myfdata,Myedata,by.x="RowNamesfdata",by.y="RowNamesedata")
    df<-df[,-1]
    df<-reshape(data=df,
                idvar="Gene_ProteinGroup",
                varying = colnames(df)[2:ncol(df)],
                times=colnames(df)[2:ncol(df)],
                v.names="Intensity",
                direction="long")
    
    #Ajoute conditions
    Myphenodata<-phenodata()[,1:2]
    df<-merge(df,Myphenodata,by.x="time",by.y="ShortName",all.x=TRUE)
    colnames(df)[1]<-"Sample"
    df<-df[!is.na(df$Intensity),]
    df
  })
  
  GraphProt<-reactive({
    req(df2(),MyCond())
    i<-input$ProtPlot_GeneName
    x<-df2()[df2()$Gene_ProteinGroup==i,]
    
    for(j in MyCond()){
      if(length(which(x$Condition==j))==0){
        y<-data.frame(Sample=NA,
                      Gene_ProteinGroup=i,
                      Intensity=0,
                      Condition=j)
        x<-rbind(x,y)
      }
    }
    
    #Calcul SD et Mean
    x.summary <- x %>%
      group_by(Condition) %>%
      summarise(
        sd = sd(Intensity, na.rm = TRUE),
        Intensity = mean(Intensity)
      )
    x.summary$Sample<-NA
    
    x<-x[x$Intensity>0,]
    
    g<-ggplot(x, aes(x=Condition, y=Intensity, label=Sample)) +
      geom_errorbar( aes(ymin = Intensity-sd, ymax = Intensity+sd), 
                     data = x.summary, width = 0.3) +
      theme_light()+
      ylab(input$ytitle)+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=90, hjust=1),
            plot.title = element_text(size=8))+
      ggtitle(i)
    
    if(input$BarOrLine=="Bar plot"){
      g<-g+geom_col(data = x.summary, fill = "grey", color = "grey50") +
        geom_point(pch=16,size=4,alpha=0.4, color = "black")
    }else{
      g<-g+geom_line(data=x.summary,aes(x=Condition,y=Intensity,group=1))+
        geom_point(data=x.summary,pch=16,size=4)
    }
    g
  })
  
  output$ProtBarplot<-renderPlotly({ GraphProt() })
  
  output$downloadProtBarPlot <-Exportggplot(graph=GraphProt(),
                                            basename=input$ProtPlot_GeneName,
                                            format=input$ProtBarPlotformat,
                                            width=input$WidthProtBarPlot,
                                            height=input$HeightProtBarPlot)
  
  
  ###### PeptidePlot ----
  GraphPep<-reactive({
    req(Pepdata())
    i<-input$ProtPlot_GeneName
    x<-Pepdata()[Pepdata()$Gene_ProteinGroup==i,]
    
    x<-reshape(data=x,
               idvar=c("Gene_ProteinGroup","Peptide"),
               varying=colnames(x)[3:ncol(x)],
               times=colnames(x)[3:ncol(x)],
               v.names="Log2Intensity",
               direction="long")
    colnames(x)[3]<-"Sample"
    x$OxyM<-" Unmodified"
    x$OxyM[str_detect(string=x$Precursor.Id,pattern="(UniMod:35)")]<-"OxyM"
    
    #Pepide plot
    g<-ggplot(x,aes(x=Sample,y=Log2Intensity,group=Peptide,color=Peptide))+
      geom_line(aes(linetype=OxyM),linewidth=1.1)+
      geom_point(size=2)+
      theme_light()+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=90, hjust=1),
            plot.title = element_text(size=8),
            legend.position = "none")+
      ggtitle(i)
    
    #Add prot data
    Myfdata<-fdata()
    y<-Myfdata[Myfdata$Gene_ProteinGroup==i,"RowNamesfdata"]
    y<-edata()[row.names(edata())==y,]
    y<-as.data.frame(t(y))
    y$Sample<-row.names(y)
    colnames(y)[1]<-"Log2Intensity"
    y$Peptide<-"Protein"
    
    g<-g+geom_line(data=y,aes(x=Sample,y=Log2Intensity),color="black",linewidth=1.2,linetype="dotted")+
      geom_point(data=y,pch=1,size=3,color="black")
    g
    
  })
  
  output$PepLinePlot<-renderPlotly({ GraphPep() })
  
  output$downloadPepLinePlot <-Exportggplot(graph=GraphPep(),
                                            basename=paste0(input$ProtPlot_GeneName,"_peptides"),
                                            format=input$PepLinePlotformat,
                                            width=input$WidthPepLinePlot,
                                            height=input$HeightPepLinePlot)
  
  ######KNIT###############

  output$MyKnit<-downloadHandler(
    filename = function() { 
        x<-dirname(input$MyPath)
        x<-gsub(pattern="/txt",replacement="",x)
        x<-gsub(pattern="/combined",replacement="",x)
        x<-unlist(strsplit(x, split = "/",fixed=TRUE))
        ProjectName<-x[length(x)]
        
        
        if(input$MyStat=="YES"){
          
          MyPaired<-if(str_split(input$MyType," ")[[1]][2]=="Paired"){".Paired"}else{""}
          PQ<-str_sub(input$ThreasholdType,1,4)
          MySignifVal<-paste0(input$SignifThreashold*100,"pc")
          MySignifFC<-if(input$FCThreashold==1){""}else{paste0("_FC",input$FCThreashold)}
          MyAlternative<-str_split(input$MyType," ")[[1]][1]
        }else{
          MyPaired<-""
          PQ<-""
          MySignifVal<-""
          MySignifFC<-""
          MyAlternative<-""
        }
        
        MyFileBase<-paste0(ProjectName,"_",MyAlternative,MyPaired,"_",PQ,MySignifVal,MySignifFC,".docx")
        MyFileBase
    },
    content = function(file) {
      rmarkdown::render("RML_DIANN_MaxQprocess3.2.Rmd", 
                        output_file = file,
                        output_format = 'word_document',
                        params = list(MyFile=MyFile(),
                                      MySoft=MySoft(),
                                      df=df(),
                                      GraphConta=GraphConta(),
                                      GraphTop5=GraphTop5(),
                                      GraphNbProt=GraphNbProt(),
                                      GraphEuclid=GraphEuclid(),
                                      GraphEuclidBox=GraphEuclidBox(),
                                      GraphPearson=GraphPearson(),
                                      GraphPearsonBox=GraphPearsonBox(),
                                      GraphLog2=GraphLog2(),
                                      MyTableSignif=if(input$MyStat=="NO"){NULL}else{MyTableSignif()},
                                      MyTableApDis=if(input$MyStat=="NO"){NULL}else{MyTableApDis()},
                                      GraphPval=if(input$MyStat=="NO"){NULL}else{GraphPval()},
                                      GraphFC=if(input$MyStat=="NO"){NULL}else{GraphFC()},
                                      edata=edata(),
                                      fdata=fdata(),
                                      phenodata=phenodata(),
                                      MyRegEx=input$MyRegEx,
                                      Pepdata=Pepdata(),
                                      PCAparams=list(TypeDataPCA=input$TypeDataPCA,
                                                     FilterType=input$FilterType,
                                                     PcFilter=input$PcFilter),
                                      MyStat=if(input$MyStat=="YES"){TRUE}else{FALSE},
                                      Ttestparams=if(input$MyStat=="NO"){NULL}else{
                                                  list(MyComp=input$MyComp,
                                                       TypeTtest=input$MyType,
                                                       ThreasholdType=input$ThreasholdType,
                                                       SignifThreashold=input$SignifThreashold,
                                                       FCThreashold=input$FCThreashold)},
                                      Volcanoparams=if(input$MyStat=="NO"){NULL}else{
                                                  list(UpColor=input$UpColor,
                                                         DownColor=input$DownColor,
                                                         NSColor=input$NSColor,
                                                         GenePoints=input$GenePoints,
                                                         PointColor=input$PointColor,
                                                         PointSize=input$PointSize,
                                                         sel_gene_nm=input$sel_gene_nm,
                                                         Force=input$Force,
                                                         ShowHline=input$ShowHline,
                                                         ShowVline=input$ShowVline,
                                                         LineType=input$LineType,
                                                         LineColor=input$LineColor)},
                                      HMsignifparams=if(input$MyStat=="NO"){NULL}else{
                                                  list(LowColorSignif=input$LowColorSignif,
                                                          MiddleColorSignif=input$MiddleColorSignif,
                                                          HighColorSignif=input$HighColorSignif,
                                                          ShowGeneName=input$ShowGeneName)},
                                      HMAnova=if(input$MyStat=="NO"){NULL}else{GraphHMAnova()}
                        )
      )
    })
  
}


####Run----
shinyApp(ui, server)
