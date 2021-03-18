
mypackages <- c("ape", "data.table", "dplyr", "shiny", "Rmisc", "shinyhelper", "gtools", "magrittr", "shinyFiles", "shinythemes", "shinyalert", "splits", "phytools", "reshape2", "devtools", "ggplot2")
checkpkg <- mypackages[!(mypackages %in% installed.packages()[,"Package"])]
if(length(checkpkg)) install.packages(checkpkg, dependencies = TRUE)

library(shiny)
library(ape)
library(splits)
library(shinyhelper)
library(magrittr) # allows you to use %>%
library(shinythemes)
library(shinyFiles)
library(reshape2)
library(shinyalert)
library(phytools)
library(ggplot2)
library(gtools)
library(Rmisc)
library(dplyr)
library(data.table)

ggthemes = list("Classic" = theme_classic(),
                "Dark" = theme_dark(),
                "Minimal" = theme_minimal(),
                "Grey" = theme_grey(),
                "Light" = theme_light(),
                "Black/White" = theme_bw(),
                "Void" = theme_void())

ui <- fluidPage(theme = shinytheme("flatly"),
                
                useShinyalert(), # set up the shinyalert package
                
                titlePanel(strong(h1("SPEDE-SAMPLER-GMYC"))),
                div(img(src='spede_sampler_R_logo.png',  height = '150px', width = '200px'), style="display: block; margin-left: auto; margin-right: auto;"),
                titlePanel(h3("Run GMYC analyses on resampled ML trees for SPEcies DElimitation", )),
                
                titlePanel(h4("Created by Clarke van Steenderen")),
                br(),
                
                tabsetPanel(
                          tabPanel(strong("Home: multiple ML trees"),
                          br(), br(),
                    #sidebarLayout(
                        sidebarPanel(width = 10,
                            h5(strong('Select the folder containing your tree files:')),
                            shinyDirButton('directory', 'Folder select', 'Please select a folder containing your tree files', style="color: black; background-color: white; border-color: black"),
                            br(), br(),
                            h5("You have selected the folder path: "),
                            br(),
                            textOutput('folder_path'),
                            br(), br(),
                            h3('OR'),
                            br(),
                            textInput(inputId = 'raw_file_path', label = 'Manually insert a file path: '),
                            br(),
                            radioButtons("data_type", label="Maximum Likelihood program used to produce tree files:", choices = c("FastTree", "RAxML")),
                          tags$div(title="Check if you want the results to be reproducible if you run this again with the same data",
                            checkboxInput("set_seed", label = "Set a seed?", value = FALSE),
                          ),
                            #checkboxInput("group_info", label="Check this box if you wish to upload predefined grouping information. If yes, upload a .csv file, and select the columns containing your grouping information and sample names from the dropdown menu.", value = FALSE),
                            
                            fileInput("predefined_groups", label="Upload a .csv file containing predefined groups for your samples:", accept = ".csv"),
                            selectInput("col.group", "Select Group Column:", choices=NULL),
                            selectInput("sample_names", "Select Sample Name Column:", choices=NULL),
                            br(),
                            actionButton("run_gmyc", label = strong("RUN"), style="color: black; background-color: lightgreen; border-color: green", icon("check-square")),
                            br(),br()
                                                           ),
                    
                            mainPanel() ),
                    
                            tabPanel(strong("View Data"),
                                     br(), br(),
                                     actionButton('all_data', label=strong('Show all data'), style="color: black; background-color: lightgreen; border-color: green", icon("edit")),
                                     downloadButton("download_clust_ent_data", label = strong("Download all data"), style="color: black; background-color: lightgreen; border-color: green"),
                                     br(), br(),
                                     actionButton('summary_data', label = strong('Show summary table'), style="color: black; background-color: lightblue; border-color: darkblue", icon("edit")),
                                     downloadButton("download_stat_summary", label = strong("Download summary table"), style="color: black; background-color: lightblue; border-color: darkblue"),
                                     br(), br(),
                                     tableOutput("data_table"),
                                     br()
                               
                            ),
         
                            tabPanel(strong("Plot Results"), 
                                     br(), br(),
                                     selectInput("ggtheme_plots", "Select ggplot Theme:", choices = names(ggthemes), selected = ggthemes["Classic"]),
                                     actionButton("plot_clusts", label = strong("Plot clusters vs entities"), style="color: black; background-color: lightgreen; border-color: green", icon("drafting-compass")), 
                                     downloadButton("download_clust_plot", label = strong("Download"), style="color: black; background-color: lightgreen; border-color: green"),
                                     br(), br(),
                                     selectInput("clust_vs_ent_plot_point_colours", "Point colour: ", choices = c("black", "blue", "red", "darkgreen")),
                                     br(), 
                                     actionButton("plot_boxplot", label = strong("Plot boxplot"), style="color: black; background-color: lightblue; border-color: darkblue", icon("drafting-compass")), 
                                     downloadButton("download_boxplot", label = strong("Download"), style="color: black; background-color: lightblue; border-color: darkblue"),
                                     br(), br(), br(),
                                     actionButton("plot_clusts_vs_iterations", label = strong("Plot clusters vs iteration file"), style="color: black; background-color: lightyellow; border-color: black", icon("drafting-compass")), 
                                     downloadButton("download_clusts_vs_iterations", label = strong("Download"), style="color: black; background-color: lightyellow; border-color: black"),
                                     br(), br(),
                                     selectInput("plot_clusts_vs_iterations_point_colours", "Point colour: ", choices = c("black", "blue", "red", "darkgreen")),
                                     selectInput("plot_clusts_vs_iterations_line_colour", "Line colour: ", choices = c("black", "grey", "lightblue", "salmon", "lightgreen", "white")),
                                     br(), br(),
                                     actionButton("plot_ents_vs_iterations", label = strong("Plot entities vs iteration file"), style="color: black; background-color: lightpink; border-color: black", icon("drafting-compass")), 
                                     downloadButton("download_ents_vs_iterations", label = strong("Download"), style="color: black; background-color: lightpink; border-color: black"),
                                     br(), br(),
                                     selectInput("plot_ents_vs_iterations_point_colours", "Point colour: ", choices = c("black", "blue", "red", "darkgreen")),
                                     selectInput("plot_ents_vs_iterations_line_colour", "Line colour: ", choices = c("black", "grey", "lightblue", "salmon", "lightgreen")),
                                     br(), br(),
                                     plotOutput("clust_ent_plot", height = "600px"),
                                     br(), br(),
                                     
                                     ),
                    
                          tabPanel(strong("Plot Trees"),
                                    br(), br(),
                                    selectInput("select_tree", label = "Select a tree to plot", choices = NULL),
                                    sliderInput("tip_label_size", label = "Tip label size: ", value = 0.8, min = 0.1, max = 2),
                                    sliderInput("support_value_size", label = "Support value size: ", value = 1, min = 0.1, max = 5),
                                    sliderInput("line_width", label = "Branch line width: ", value = 1, min = 0.1, max =5),
                                    selectInput("support_value_type", label = "Select which support values to display: ", choices = c("Original ML", "GMYC estimate"), selected = "GMYC estimate"),
                                    selectInput("support_value_col", label = "Support value colour: ", choices = c("grey", "lightblue", "salmon", "lightgreen", "lightyellow", "white"), selected = "lightgreen"),
                                    selectInput("support_value_frame", label = "Support value frame:", choices = c("none", "circle", "rect"), selected = "rect"),
                                    selectInput("branch_col", label = "Branch colour: ", choices = c("black", "blue", "lightblue", "red", "green", "orange"), selected = "blue"),
                                    br(),
                                    actionButton("plot_gmyc_tree", label = strong("Plot GMYC tree result"), style="color: black; background-color: lightpink; border-color: black", icon("drafting-compass")),
                                    downloadButton("download_gmyc_tree", label = strong("Download"), style="color: black; background-color: lightblue; border-color: darkblue"),
                                    br(), br(),
                                    plotOutput("gmyc_tree", height = "1250px"),
                                    br(), br()
                          ),
                      
                          tabPanel(strong("Percentage Matches"),
                                    br(), br(),
                                    selectInput("select_tree_speclist", label = "Select a tree file", choices = NULL),
                                   
                                  tags$div(title="View the GMYC analysis output for each tree file uploaded",
                                    actionButton("view_gmyc_spec", label = strong("View GMYC species list"), style="color: black; background-color: lightyellow; border-color: black", icon("edit")), 
                                    br(), br(),
                                  ),
                                    
                                    downloadButton("download_gmyc_spec", label = strong("Download"), style="color: black; background-color: lightyellow; border-color: black"),
                                    br(), br(),
                                  tags$div(title="View the % match between what the GMYC analysis considers a species, to the predefined groups you uploaded",
                                    actionButton("view_match_data", label = strong("View Matches"), style="color: black; background-color: lightgreen; border-color: green", icon("edit")),
                                    br(), br(),
                                  ),
                                    
                                    downloadButton("download_match_data", label = strong("Download"), style="color: black; background-color: lightgreen; border-color: green"),
                                    br(), br(),
                                  tags$div(title="View the average, standard deviation, minimum and maximum values for the GMYC data and % matches",
                                    actionButton("view_summary_match_data", label = strong("View Matches Summary"), style="color: black; background-color: lightblue; border-color: darkblue", icon("edit")),
                                    br(), br(),
                                  ),
                                    downloadButton("download_match_data_summary", label = strong("Download"), style="color: black; background-color: lightblue; border-color: darkblue"), 
                                    br(),br(),
                                    tableOutput("matches"),
                                    br(), br()
                          
                        ),
                    
                          tabPanel(strong("Plot Percentage matches"),
                                   br(), br(),
                                   actionButton("plot_matches", label = strong("Plot"), style="color: black; background-color: lightyellow; border-color: black", icon("drafting-compass")),
                                   downloadButton("download_match_plot", label = strong("Download"), style="color: black; background-color: lightyellow; border-color: black"),
                                   br(), br(),
                                   selectInput("plot_matches_point_colours", "Point colour: ", choices = c("black", "blue", "red", "darkgreen")),
                                   selectInput("plot_matches_line_colour", "Line colour: ", choices = c("black", "grey", "lightblue", "salmon", "lightgreen", "white")),
                                   selectInput("ggtheme_matches", "Select ggplot Theme:", choices = names(ggthemes), selected = ggthemes["Classic"]),
                                   br(), br(),
                                   plotOutput("match_plot", height = "600px"),
                                   br(), br()
                                   
                        ),
                    
                        tabPanel(strong("GMYC Oversplitting"),
                                 br(), br(),
                                 actionButton("GMYC_oversplit_table_view", label = strong("View Summary Table"), style="color: black; background-color: lightblue; border-color: darkblue"),
                                 downloadButton("GMYC_oversplit_table_download", label = strong("Download"), style="color: black; background-color: lightblue; border-color: darkblue"),
                                 br(), br(),
                                 actionButton("GMYC_oversplit_boxplot", label = strong("Boxplot"), style="color: black; background-color: lightyellow; border-color: black"),
                                 downloadButton("GMYC_oversplit_boxplot_download", label = strong("Download"), style="color: black; background-color: lightyellow; border-color: black"),
                                 br(), br(),
                                 actionButton("GMYC_oversplit_barplot", label = strong("Bar chart"), style="color: black; background-color: lightgreen; border-color: darkgreen"),
                                 downloadButton("GMYC_oversplit_barplot_download", label = strong("Download"), style="color: black; background-color: lightgreen; border-color: darkgreen"),
                                 br(), br(),
                                 selectInput("GMYC_oversplit_ggtheme", "Select ggplot Theme:", choices = names(ggthemes), selected = ggthemes["Classic"]),
                                 selectInput("GMYC_barchart_fill", "Barchart fill: ", choices = c("black", "lightgrey", "white", "lightblue", "lightgreen"), selected = "white"),
                                 selectInput("GMYC_barchart_outline", "Barchart outline: ", choices = c("black", "darkblue", "white", "darkgreen"), selected = "black"),
                                 br(), br(),
                                 tableOutput("GMYC_oversplit_table"),
                                 br(), br(),
                                 plotOutput("GMYC_oversplit_plot"),
                                 br(), br()
                                 
                        ),
                    
                        tabPanel(strong("Plot for multiple-column data"),
                                 br(), br(),
                                 fileInput("multiple_input", label = "Upload a .csv file with multiple columns of output data:", accept = ".csv", multiple = TRUE),
                                 
                                 textInput("title_multiple_input", label = "Title: ", value = "Title", width = "600px"),
                                 textInput("x_lab_multiple_input", label = "X-axis label: ", value = "Resampled data (%)", width = "600px"),
                                 textInput("y_lab_multiple_input", label = "Y-axis label: ", value = "Measure variable", width = "600px"),
                                 #textInput("y_interval_multiple_input", label = "Y-axis tick-mark interval: ", value = "10", width = "200px"),
                                 numericInput("y_interval_multiple_input", label = "Y-axis tick-mark interval: ", value = 1, min = 1),
                                 selectInput("ggtheme_multiple", "Select ggplot Theme:", choices = names(ggthemes), selected = ggthemes["Classic"]),
                                 hr(),
                                 actionButton("multiple_input_boxplot", label = strong("Boxplot"), style="color: black; background-color: lightyellow; border-color: black", icon("drafting-compass")),
                                 downloadButton("download_multiple_input_boxplot", label = strong("Download"), style="color: black; background-color: lightyellow; border-color: black", icon("drafting-compass")),
                                 br(), br(),
                                 actionButton("plot_multiple_input", label = strong("Line Plot"), style="color: black; background-color: lightblue; border-color: darkblue", icon("drafting-compass")),
                                 downloadButton("download_multiple_input_plot", label = strong("Download"), style="color: black; background-color: lightblue; border-color: darkblue"),
                                 br(), br(),
                                 
                                 hr(),
                                 h4(strong("LINE PLOT TWEAKS:")),
                                 
                                 br(), 
                                 selectInput("error_bar_type", label = "Select what error bars should represent: ", choices = c("ci", "sd", "se")),
                                 selectInput("multiple_input_line_type", label = 'Select a line type: ', choices = c("blank" = 0, "solid" = 1, "dashed" = 2, "dotted" = 3), selected = 1),
                                 selectInput('multiple_input_line_col', label = "Select a line colour: ", choices = c("black", "grey", "lightblue", "salmon", "lightgreen")),
                                 sliderInput("multiple_input_line_width", label = "Line width: ", value = 1, min = 0.5, max = 5, step = 0.5),
                                 selectInput("multiple_input_error_bar_color", label = "Select a colour for the error bars: ", choices = c("grey", "lightblue", "black", "white")),
                                 sliderInput("multiple_input_point_size", label = "Point size: ", value = 1, min = 1, max = 10, step = 0.5),
                                 selectInput("multiple_input_point_shape", label = "Point shape: ", choices = c("Round filled" = 16, "Round open" = 1, "+" = 3, "X" = 4, "Square" = 15, "Triangle" = 17, "Diamond" = 18)),
                                 selectInput("multiple_input_point_colour", label = "Point colour: ", choices = c("blue", "lightblue", "red", "green", "forestgreen", "black", "grey"), selected = "black"),
                                 
                                 br(), br(),
                                 plotOutput("multiple_input_plot", height = "600px"),
                                 br(), br()
                                
                                 
                        )
                            
                )
)

#)