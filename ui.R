mypackages <- c("shiny", "shinyhelper", "magrittr", "pvclust", "rhandsontable", "shinyFiles", "MASS")
checkpkg <- mypackages[!(mypackages %in% installed.packages()[,"Package"])]
if(length(checkpkg)) install.packages(checkpkg, dependencies = TRUE)


library(shiny)
library(readr)
library(shinyhelper)
library(magrittr) # allows you to use %>%
library(pvclust)
library(MASS)
library(rhandsontable)
library(shinyFiles)


ui <- fluidPage(
  
  titlePanel(strong("BINMAT: FOR FRAGMENT ANALYSIS DATA")),
  img(src='clevercow.png', height = '150px', width = '450px'),
  titlePanel(h4("Created by Clarke van Steenderen")),
  br(),
  sidebarLayout(
    sidebarPanel( strong("Click the ? to view the help file:", style = "color:blue")
                  %>% helper(type = "markdown", content = "BINMAT_HelpFile", colour = "blue", icon = "question-circle"),
                  br(),
                  #radioButtons("format", "Select data format:", choices = c(".csv", ".xlsx"), selected = ".csv"),
                  #radioButtons("choice", "What type of data do you have?", choices = c("The whole dataset are replicates for one species/location" =1, "Each pair in the dataset are replicates" =2), selected = 2),
                  fileInput("inFile", "Select a .csv file", accept = c(".csv")),
                 # actionButton("preview", "PREVIEW", style="color: #fff; background-color: black; border-color: white", icon("thumbs-o-up")),
           
                  actionButton("check", "Check my data for unwanted values", style="color: #fff; background-color: black; border-color: white"),
                  br(), br(),
                  actionButton("act", "Consolidate matrix", style="color: #fff; background-color: black; border-color: white", icon("magic")),
                  br(), br(),
                  htmlOutput("consod_done"), tags$head(tags$style("#consod_done{color: darkred; font-size: 20px;}")),
                  br(), br(),
                  downloadButton("downloadData", "Download Matrix", style="color: #fff; background-color: darkblue; border-color: white"),
                  br(), 
                  tableOutput("text1") # check for unwanted values
                  
    ),
    
    mainPanel(
    # these are separate 'warnings' merely so that the output are on separate lines
      tabsetPanel(tabPanel(strong("MATRIX"),
                           br(),
                           
                           htmlOutput("warning5"), tags$head(tags$style("#warning5{color: blue; font-size: 15px;}")),
                           br(),
                           htmlOutput("warning4"), tags$head(tags$style("#warning4{color: blue; font-size: 15px;}")),
                           htmlOutput("warning"), tags$head(tags$style("#warning{color: blue; font-size: 15px;}")),
                           htmlOutput("warning2"), tags$head(tags$style("#warning2{color: blue; font-size: 15px;}")),
                           htmlOutput("warning3"), tags$head(tags$style("#warning3{color: blue; font-size: 15px;}")),
                           tableOutput("table")
                           
                           #this changes the text font and colour of a particular output text. In this case,
                           # the 'warning' outputs 
                           
      ),
      
      tabPanel(strong("SUMMARY"),
               br(),
               actionButton("summary", "Summary info", style="color: #fff; background-color: black; border-color: white"),
               br(), br(),
               tableOutput("text3"), # summary info
               downloadButton("download_summary", "Download Summary", style="color: #fff; background-color: darkblue; border-color: white")
      ),
      
      tabPanel(strong("ERROR RATES"),
               br(),
               actionButton("repro", "Check Error rates", style="color: #fff; background-color: black; border-color: white"),
               br(), br(),
               tableOutput("text2"), # error rates
               downloadButton("download_errors", "Download Errors", style="color: #fff; background-color: darkblue; border-color: white"),
               br(), br(),
               numericInput("jacc_remove", "Remove samples with a jaccard error greater than:", 0, step = 0.1), 
               br(),
               actionButton("remove", "Remove", style="color: #fff; background-color: black; border-color: white", icon("times")),
               br()
      ),
      
      tabPanel(strong("UPGMA TREE"), 
               br(),
               fileInput("inFile_upgma", "Select a consolidated .csv file", accept = ".csv"),
               numericInput("boot", "Number of bootstrap replicates for clustering tree:", 10, step = 1, min = 1),
               actionButton("drawTree", "Create Clustering Tree", style="color: #fff; background-color: black; border-color: white", icon("pencil")),
               downloadButton("downloadTree", "Download Clustering Tree", style="color: #fff; background-color: darkblue; border-color: white"),
               br(), br(),
               plotOutput("upgmaTree", height = "800px", width= "1500px")
               
      ),
      
      tabPanel(strong("nMDS PLOT"),
               br(),
               strong("Upload an already consolidated binary matrix."),
               br(), br(),
               fileInput("inFile_mds", "Select a .csv file", accept = c(".csv")),
               selectInput("dist_method", "Select a distance method", choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "binary"),
               htmlOutput("warn1"), tags$head(tags$style("#warn1{color: blue; font-size: 15px;}")),
               htmlOutput("warn2"), tags$head(tags$style("#warn2{color: blue; font-size: 15px;}")),
               htmlOutput("warn3"), tags$head(tags$style("#warn3{color: blue; font-size: 15px;}")),
               rHandsontableOutput("grps"),
               br(),
               br(),
               numericInput("k_value", "No. of dimensions (k): ", 2, step = 1, min = 2, max = 3),
               br(),
               actionButton("plotMDS", "PLOT nMDS", style="color: #fff; background-color: darkred; border-color: white; font-size:130%", icon("pencil")),
               br(), br(),
              # downloadButton("download_subsetted_dat", "Download filtered data set", style="color: #fff; background-color: darkblue; border-color: white"),
              # downloadButton("download_removed_data", "Download removed samples", style="color: #fff; background-color: darkblue; border-color: white"),
               downloadButton("download_dist_matrix","Download distance matrix", style="color: #fff; background-color: darkblue; border-color: white"),
               br(), br(), br(),
               sliderInput("cexSize", "Point size", min = 0.1, max = 5, step = 0.1, value = 1),
               radioButtons("display_labs", "Show sample labels?", choices = c("No", "Yes"), selected = "No"),
               plotOutput("mdsPlot", height = "600px", width = "700px"),
               downloadButton("downloadMDS", "Download MDS Plot", style="color: #fff; background-color: darkblue; border-color: white"),
               br(), br()
               
               
      ),
      
      tabPanel(strong("nMDS Validation"),
               
               br(),
               actionButton("scree", "Scree Plot", style="color: #fff; background-color: darkblue; border-color: white", icon("pencil")),
               actionButton("shep", "Shepard Plot", style="color: #fff; background-color: darkblue; border-color: white", icon("pencil")),
               plotOutput("screePlot", height = "600px", width = "700px"),
               br(),
               downloadButton("download_scree", "Download Scree Plot", style="color: #fff; background-color: darkblue; border-color: white"),
               br(),
               plotOutput("shepPlot", height = "600px", width = "700px"),
               br(),
               downloadButton("download_shep", "Download Shepard Plot", style="color: #fff; background-color: darkblue; border-color: white"),
               br(), br()
               
      ),
      
      tabPanel(strong("Filter data"),
               br(), br(),
               numericInput("peak_thresh", "Remove samples with peaks less than:", 0, step = 1, min = 0),
               htmlOutput("msg"), tags$head(tags$style("#warn1{color: blue; font-size: 15px;}")),
               br(), br(),
               actionButton("filter_samples", "CHECK", style="color: #fff; background-color: darkred; border-color: white; font-size:130%", icon("thumbs-o-up")),
               br(), br(),
               downloadButton("download_filtered", "Download filtered data", style="color: #fff; background-color: darkblue; border-color: white"),
               downloadButton("download_removed", "Download removed samples",  style="color: #fff; background-color: darkblue; border-color: white"),
               br(), br()
      )
      
      
      )
      
    )       
    
  )
)

