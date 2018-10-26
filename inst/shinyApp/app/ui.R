library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)

ui <- dashboardPage(
    
    dashboardHeader
    (
        title="Analysis Tool"
    ),
    
    dashboardSidebar
    (
        useShinyjs(),
        sidebarMenu
        (
            id="tabs",
            menuItem("Files Selection", tabName="t_1"),
            menuItem("Scoring", tabName="t_3")
        )
    ),
    
    dashboardBody
    (
        useShinyjs(),
        tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
        tabItems
        (
            tabItem
            (
                tabName="t_1",
                h2("Files Selection"),
                fluidRow
                (
                    id="t_1_fr",
                    box
                    (
                        title="Project Selection", id="t_1_1",width=9,
                        selectInput("t_1_1_select_project", "Select Project", choices=NULL),
                        actionButton("t_1_1_save_project", "Save Project"),
                        actionButton("t_1_1_load_project", "Load Project"),
                        actionButton("t_1_1_remove_project", "Remove Project", style="float:right;background-color:gray;color:orange"),
                        textInput("t_1_1_name_project","Project Name (creation only)",value = "[project name]"),
                        actionButton("t_1_1_create_project", "Create New Project")
                    ),
					box
					(
						id="t_1_2",width=3,
						actionButton("t_1_2_add", "Add files",style="width:90%;margin-left:4.8%"),
						downloadButton("t_1_2_dl", "Download Selection",style="width:90%;margin-left:4.8%"),
						actionButton("t_1_2_rm", "Remove Selection",style="width:90%;margin-left:4.8%")
					),
                    fluidRow
                    (
                        id="t_1_3"
                    )
                )
            ),
            
            tabItem
            (
                tabName="t_3",
                h2("Scoring"),
                fluidRow
                (
                    id="t_3_fr",
                    box
					(
						height = "8vh", width=12, id="t_3_2",
						actionButton("t_3_2_run", "Run Analysis",width="100%")
					),
                    shinyjs::hidden(fluidRow
                    (
                        id="t_3_3_fr",
                        tabBox
                        (
                            id="t_3_3",width=12,
                            tabPanel
                            (
                                "Annotations Visualization",id="t_3_3_1",
                                box
                                (
                                    width=6,
                                    selectInput("t_3_3_1_fileSel", "Select file", choices = NULL, selected = NULL),
                                    actionButton("t_3_3_1_plotButton", "Plot")
                                ),
                                box
                                (
                                    width=6,
                                    selectInput("t_3_3_1_methodSel", "Select Algorithm", choices = NULL, selected = NULL),
                                    selectInput("t_3_3_1_runSel", "Select Run", choices = NULL, selected = NULL)
                                ),
                                fluidRow
                                (
                                    id="t_3_3_1_fr",
                                    box
                                    (
                                        id="t_3_3_1_ref",
                                        imageOutput("t_3_3_1_refPlot"),
                                        div
                                        (
                                            selectInput("t_3_3_1_refPopSel", "Select Populations", choices=NULL, multiple = T),
                                            style="margin-top:10%"
                                        ),
                                        width=5
                                    ),
                                    box
                                    (
                                        id="t_3_3_1_test",
                                        imageOutput("t_3_3_1_testPlot"),
                                        div
                                        (
                                            selectInput("t_3_3_1_testPopSel", "Select Annotated Groups", choices=NULL, multiple = T),
                                            style="margin-top:10%"
                                        ),
                                        width=5
                                    ),
                                    box
                                    (
                                        id="t_3_3_1_markers",
                                        selectInput("t_3_3_1_m1", "Marker 1", choices = NULL),
                                        selectInput("t_3_3_1_m2", "Marker 2", choices = NULL),
                                        width=2
                                    )
                                )
                            ),
                            tabPanel
                            (
                                "F score details",id="t_3_3_2",
                                box
                                (
                                    width=12,
                                    selectInput("t_3_3_2_fileSel", "Select file", choices = NULL, selected = NULL),
                                    selectInput("t_3_3_2_methodSel", "Select Algorithm", choices = NULL, selected = NULL),
                                    actionButton("t_3_3_2_plotButton", "Plot")
                                ),
                                fluidRow
                                (
                                    box
                                    (
                                        id="t_3_3_2_sumScoreByParam", width=6,style="overflow:auto;height:65vh;max-height:75vh",
                                        selectInput("t_3_3_2_paramSel", "Select Variable Parameter", choices = NULL, selected = NULL),
                                        selectInput("t_3_3_2_FixedparamSel", "Select Fixed Values", choices = NULL, selected = NULL),
                                        imageOutput("t_3_3_2_paramSelPlot")
                                    ),
                                    box
                                    (
                                        id="t_3_3_2_scoreByPOP", width=6,style="overflow:auto;height:65vh;max-height:75vh",
                                        selectInput("t_3_3_2_runSel", "Select run", choices = NULL, selected = NULL),
                                        imageOutput("t_3_3_2_runSelPlot")
                                    )
                                )
                            ),
                            tabPanel
                            (
                                "Summarize table",id="t_3_3_3",
                                actionButton("t_3_3_3_plotButton", "Plot"),
                                fluidRow
                                (
                                    id="t_3_3_3_fr",
                                    box
                                    (
                                        width=12, style="overflow:auto",
                                        tableOutput("t_3_3_3_table")
                                    )
                                )
                            ),
                            tabPanel
                            (
                                "Advanced Options",id="t_3_3_4",
                                fluidRow
                                (
                                    box
                                    (
                                        selectInput("t_3_3_4_fileSel", "Select file", choices = NULL, selected = NULL),
                                        selectInput("t_3_3_4_methodSel", "Select Algorithm", choices = NULL, selected = NULL),
                                        actionButton("t_3_3_4_plotButton", "Plot")
                                    ),
                                    box
                                    (
                                        selectInput("t_3_3_4_runSel", "Select Run", choices = NULL, selected = NULL),
                                        downloadButton("t_3_3_4_exportButton", "Export Data")
                                    )
                                ),
                                
                                fluidRow
                                (
                                    id="t_3_3_4_purityByAnnot_fr",
                                    box
                                    (
                                        imageOutput("t_3_3_4_purityByAnnot"),
                                        div
                                        (
                                            sliderInput("t_3_3_4_purityByAnnot_slider", "Min Purity Threshold", min=0, max=1, val=0.5),
                                            style="margin-top:10%;display:inline-block;width:100%"
                                        ),
                                        style="height:65vh;overflow:auto", collapsible=T
                                    ),
                                    box
                                    (
                                        div
                                        (
                                            tableOutput("t_3_3_4_clustersDetailsBelow"),
                                            style="height:48%;overflow:auto"
                                        ),
                                        div
                                        (
                                            tableOutput("t_3_3_4_clustersDetailsAbove"),
                                            style="height:48%;overflow:auto"
                                        ),
                                        style="height:65vh;overflow:auto", collapsible=T
                                    ),
                                    box
                                    (
                                        DT::dataTableOutput("t_3_3_4_MapAnnotPop"),
                                        style="height:40vh;overflow:auto", collapsible=T
                                    ),
                                    box
                                    (
                                        DT::dataTableOutput("t_3_3_4_populationsDetails"),
                                        style="height:40vh;overflow:auto", collapsible=T
                                    )
                                )
                            )
                        )
                    ))
                )
            )
            # tabItem
            # (
            #     tabName="t_4",
            #     h2("Visualization")
            # )
            
        )
        
    )
    
)