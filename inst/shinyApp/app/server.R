library(shiny)
library(shinydashboard)
library(flowCore)
library(shinyjs)
library(DT)
library(parallel)
library(doSNOW)
library(ggplot2)

server <- function(input, output, session)
{
    "%not.in%" = Negate("%in%")
    useShinyjs()
    # options(shiny.reactlog=TRUE) 

    #======================================================================================================================
    #======================REACTIVE VALUES=================================================================================
    #======================================================================================================================
    
    current.project <- reactiveValues(
        name = NULL,
        fcs.files = NULL,
        fcs.files.ui.colnames = NULL,
        mapping.files = NULL,
        time.results = NULL,
        scores.results = NULL,
        modified.files = NULL,
        ref.files.populations.col  = NULL,
        test.files.clusters.col = NULL,
        verbose = T,
        nmb.cores = detectCores()
    )
    
    clustering.algorithms <- reactiveValues(
        algorithms = NULL,
        parameters = NULL,
        run.analysis = FALSE
    )
    
    computed.values <- reactiveValues(
        purity.matrix.annot = NULL,
        purity.matrix.clust = NULL,
        FG.matrices.annot = NULL,
        FG.matrices.clust = NULL,
        prec.rec.matrices.annot = NULL,
        prec.rec.matrices.clust = NULL,
        annot.sizes = NULL,
        clust.sizes = NULL,
        pop.sizes = NULL,
        summary.table = NULL,
        list.pop.points = NULL,
        list.pop.points.xval = NULL,
        ordered.table = NULL,
        fixed.parameters.values = NULL,
        fixed.parameters.ids = NULL
    )
    
    env.var <- reactiveValues(
        tool.wd = getwd()
    )
    
    write.enriched.FCS <- function(fcs, fcs.path)
    {
        keywords.to.save <- names(get.keywords.with.keypart.FCS(fcs, "MAPOP_pop_label"))
        keywords.to.save <- c(unlist(keywords.to.save), names(get.keywords.with.keypart.FCS(fcs, "EXPPUR__")))
        keywords.to.save <- c(unlist(keywords.to.save), names(get.keywords.with.keypart.FCS(fcs, "RF_pop_label")))
        keywords.to.save <- c(unlist(keywords.to.save), names(get.keywords.with.keypart.FCS(fcs, "CLMETH__")))
        
        write.FCS.CIPHE(fcs,fcs.path, keywords.to.save = keywords.to.save)
    }
    
    
    observe(
    {
        if(current.project$verbose)
        {
            options(warn = -1)
        }
    })
    
    
    
    
    #======================================================================================================================
    #======================================================================================================================
    #==========================================LOAD FILES==================================================================
    #======================================================================================================================
    
    
    
    load.project <- function(loaded.project = NULL)
    {
        shinyjs::disable("t_1_1_load_project")
        progress <- shiny::Progress$new()
        progress$set(message="LOADING PROJECT",value=0)
        on.exit(progress$close())
        
        clustering.algorithms$run.analysis <- FALSE
        
        current.project$fcs.files <- list()
        current.project$fcs.files.ui.colnames <- list()
        current.project$mapping.files <- list()
        current.project$modified.files <- list()
        removeUI("#t_1_3")
        insertUI("#t_1_fr",
                 "beforeEnd",
                 fluidRow
                 (
                     id="t_1_3"
                 ))
        
        if(is.null(loaded.project))
        {
            if( is.defined(input$t_1_1_select_project) && (input$t_1_1_select_project != "") && (input$t_1_1_select_project != " "))
            {
                current.project$name <- input$t_1_1_select_project
                if(dir.exists(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/")))
                {
                    files <- list.files(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/"), full.names = F)
                    lapply(files, function(f.name)
                    {
                        x <- NULL
                        nx <- list()
                        ext <- substr(f.name, nchar(f.name)-2, nchar(f.name))
                        if( ext == "fcs" )
                        {
                            x <- read.FCS(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/",f.name),emptyValue = FALSE)
                            for(i in 1:ncol(x@exprs))
                            {
                                d <- x@description[[paste0("$P",i,"S")]]
                                if(is.null(d) || is.na(d) || d == "" || d == " " || d == "NA" || d == "<NA>" || d == "'<NA>'")
                                {
                                    d <- colnames(x)[i]
                                }
                                nx[[i]] <- d
                            }
                        }
                        
                        if(is.defined(x))
                        {
                            register.name <- substr(f.name,1,nchar(f.name)-4)
                            
                            current.project$mapping.files[[register.name]] <<- NA
                            current.project$fcs.files[[register.name]] <<- x
                            current.project$fcs.files.ui.colnames[[register.name]] <<- nx
                            current.project$modified.files[[register.name]] <<- TRUE
                            
                            t <- substr(f.name,1,nchar(f.name)-4)
                            if(file.exists(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/",t,".csv")))
                            {
                                map.name <- paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/",t,".csv")
                                current.project$mapping.files[[register.name]] <<- as.matrix(read.csv(map.name))
                            }
                            
                        }
                    })
                }
            }
        }
        else
        {
            current.project$name <- loaded.project
            if(dir.exists(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/")))
            {
                files <- list.files(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/"), full.names = F)
                lapply(files, function(f.name)
                {
                    x <- NULL
                    nx <- list()
                    ext <- substr(f.name, nchar(f.name)-2, nchar(f.name))
                    if( ext == "fcs" )
                    {
                        x <- read.FCS(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/",f.name),emptyValue = FALSE)
                        lapply(1:ncol(x@exprs), function(i)
                        {
                            d <- x@description[[paste0("$P",i,"S")]]
                            if(is.null(d) || is.na(d) || d == "" || d == " " || d == "NA" || d == "<NA>" || d == "'<NA>'")
                            {
                                d <- colnames(x)[i]
                            }
                            nx[[i]] <<- d
                        })
                    }
                    
                    if(is.defined(x))
                    {
                        register.name <- substr(f.name,1,nchar(f.name)-4)
                        
                        current.project$mapping.files[[register.name]] <<- NA
                        current.project$fcs.files[[register.name]] <<- x
                        current.project$fcs.files.ui.colnames[[register.name]] <<- nx
                        current.project$modified.files[[register.name]] <<- TRUE
                        
                        t <- substr(f.name,1,nchar(f.name)-4)
                        if(file.exists(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/",t,".csv")))
                        {
                            map.name <- paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/",t,".csv")
                            current.project$mapping.files[[register.name]] <<- as.matrix(read.csv(map.name))
                        }
                        
                    }
                })
            }
        }
        
        if(length(current.project$fcs.files)==0)
        {
            current.project$fcs.files <- NULL
            current.project$fcs.files.ui.colnames <- NULL
        }
        progress$set(message="DONE",value=1)
        
        shinyjs::delay(500,shinyjs::enable("t_1_1_load_project"))
    }
    
    update.project.list <- function()
    {
        projects <- list.dirs(paste0(env.var$tool.wd,"/Projects/"), full.names = F, recursive = F)
        if(length(projects) > 0)
        {
            shinyjs::enable("t_1_1_load_project")
            names(projects) <- projects
            updateSelectInput(session,"t_1_1_select_project", "Select Project", choices=projects, selected = projects[[1]])
        }
        else
        {
            shinyjs::disable("t_1_1_load_project")
        }
    }
    
    remove.files <- function(files.id = NULL)
    {
        if(!is.null(files.id))
        {
            for(i in files.id)
            {
                file.remove(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/",
                                   names(current.project$fcs.files)[i],".fcs"))
                file.remove(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/",
                                   names(current.project$fcs.files)[i],".csv"))
                current.project$fcs.files[[i]] <- NA
                current.project$fcs.files.ui.colnames[[i]] <- NA
                current.project$mapping.files[[i]] <- NA
                current.project$modified.files[[i]] <- T
                removeUI(paste0("#t_1_3_fr_",i))
            }
        }
        else
        {
            if( !is.null(current.project$name) && length(current.project$fcs.files) >0 )
            {
                rm.files.ids <- c()
                for(i in 1:length(current.project$fcs.files))
                {
                    file.remove(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/",
                                       names(current.project$fcs.files)[i],".fcs"))
                    file.remove(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/",
                                       names(current.project$fcs.files)[i],".csv"))
                    current.project$fcs.files[[i]] <- NA
                    current.project$fcs.files.ui.colnames[[i]] <- NA
                    current.project$mapping.files[[i]] <- NA
                    current.project$modified.files[[i]] <- T
                    removeUI(paste0("#t_1_3_fr_",i))
                }
            }
        }
    }
    
    observe(#SEARCH FOR PROJECTS
    {
        update.project.list()
    })
    
    observe(#ACTIVATE SAVE UI
    {
        if(!is.null(current.project$name) && length(current.project$fcs.files)>0)
        {
            nmb.na.files <- 0
            lapply(1:length(current.project$fcs.files), function(f)
            {
                fcs <- current.project$fcs.files[[f]]
                
                if(!is.defined(fcs))
                {
                    nmb.na.files <<- nmb.na.files+1
                }
            })
            
            if(nmb.na.files >= length(current.project$fcs.files))
            {
                shinyjs::disable("t_1_1_save_project")
            }
            else
            {
                shinyjs::enable("t_1_1_save_project")
            }
        }
        else
        {
            shinyjs::disable("t_1_1_save_project")
        }
    })
    
    observe(#ACTIVATE LOAD/REMOVE UI
    {
        if(is.defined(current.project$name))
        {
            shinyjs::enable("t_1_1_load_project")
            shinyjs::enable("t_1_1_remove_project")
        }
        else
        {
            shinyjs::disable("t_1_1_laod_project")
            shinyjs::disable("t_1_1_remove_project")
        }
    })
    
    observeEvent(input$t_1_1_save_project,#SAVE PROJECT
    {
        shinyjs::disable("t_1_1_saveproject")
        progress <- shiny::Progress$new()
        progress$set(message="SAVING PROJECT",value=0)
        on.exit(progress$close())
        if(!is.null(current.project$name) && length(current.project$fcs.files)>0)
        {
            if(!dir.exists(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/")))
            {
                dir.create(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/"))
            }
            lapply(1:length(current.project$fcs.files), function(f)
            {
                if(is.defined(current.project$fcs.files[[f]]))
                {
                    t <- names(current.project$fcs.files)[f]
                    if(!file.exists(paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/",t,".fcs")))
                    {
                        t <- paste0(names(current.project$fcs.files)[f], "__", f,"_",length(current.project$fcs.files)+1)
                    }
                    print(t)
                    
                    if( is.defined(input[[paste0("t_1_3_",current.project$name,"_",f,"_pop_col")]]))
                    {
                        current.project$fcs.files[[f]] <<- add.keyword.to.fcs(current.project$fcs.files[[f]], input[[paste0("t_1_3_",current.project$name,"_",f,"_pop_col")]], "RF_pop_label")
                    }

                    if( length(current.project$mapping.files)>0 && is.defined(current.project$mapping.files[[f]]) )
                    {
                        if( is.defined(input[[paste0("t_1_3_",current.project$name,"_",f,"_lab_col")]]) &&
                            input[[paste0("t_1_3_",current.project$name,"_",f,"_lab_col")]] != "" &&
                            input[[paste0("t_1_3_",current.project$name,"_",f,"_lab_col")]] != " ")
                        {
                            new.nmb.row <- length(unique(current.project$fcs.files[[f]]@exprs[,as.numeric(input[[paste0("t_1_3_",current.project$name,"_",f,"_pop_col")]])]))
                            pop.written.names <- sapply(1:new.nmb.row, function(cl)
                            {
                                return(input[[paste0("t_1_3_",current.project$name,"_",f,"_1_pop_",cl)]])
                            })

                            current.project$fcs.files[[f]] <<- add.keyword.to.fcs(current.project$fcs.files[[f]],
                                                                                  input[[paste0("t_1_3_",current.project$name,"_",f,"_lab_col")]],
                                                                                  "MAPOP_pop_label")

                            if(new.nmb.row != nrow(current.project$mapping.files[[f]]))
                            {
                                current.project$mapping.files[[f]] <<- matrix(nrow=new.nmb.row,ncol=1)
                                mat <- matrix(1:new.nmb.row, ncol=1)
                                mat <- cbind(mat, mat)
                                colnames(mat) <- c("popID","popName")
                                current.project$mapping.files[[f]] <<- mat
                            }
                            current.project$mapping.files[[f]][,as.integer(input[[paste0("t_1_3_",current.project$name,"_",f,"_lab_col")]])] <<- pop.written.names
                        }
                        else
                        {
                            pop.written.names <- FPH.get.labels.from.mapping.file(current.project$mapping.files[[f]],1)
                            current.project$fcs.files[[f]] <<- add.keyword.to.fcs(current.project$fcs.files[[f]], 1, "MAPOP_pop_label")
                            current.project$mapping.files[[f]] <<- cbind(current.project$mapping.files[[f]],pop.written.names)
                            colnames(current.project$mapping.files[[f]])[ncol(current.project$mapping.files[[f]])] <-
                                paste0("added_",ncol(current.project$mapping.files[[f]]))
                        }
                        write.csv(current.project$mapping.files[[f]],paste0(env.var$tool.wd,"/Projects/",
                                                                            current.project$name,
                                                                            "/data/",
                                                                            t,
                                                                            ".csv"),
                                  row.names = F)
                    }
                    else
                    {
                        if( is.defined(input[[paste0("t_1_3_",current.project$name,"_",f,"_pop_col")]]) )
                        {
                            mat <- matrix(1:length(FPH.get.file.clusters(current.project$fcs.files[[f]],
                                                                         as.numeric(input[[paste0("t_1_3_",current.project$name,"_",f,"_pop_col")]]))),
                                          ncol=1)
                            new.nmb.row <- length(unique(current.project$fcs.files[[f]]@exprs[,as.numeric(input[[paste0("t_1_3_",current.project$name,"_",f,"_pop_col")]])]))
                            pop.written.names <- sapply(1:new.nmb.row, function(cl)
                            {
                                return(input[[paste0("t_1_3_",current.project$name,"_",f,"_1_pop_",cl)]])
                            })

                            mat <- cbind(mat, pop.written.names)
                            colnames(mat) <- c("popID","popName")
                            current.project$mapping.files[[names(current.project$fcs.files)[f]]] <<- mat
                            write.csv(current.project$mapping.files[[names(current.project$fcs.files)[f]]],paste0(env.var$tool.wd,"/Projects/",
                                                                                                                  current.project$name,
                                                                                                                  "/data/",
                                                                                                                  t,
                                                                                                                  ".csv"),
                                      row.names = F)

                            current.project$fcs.files[[f]] <<- add.keyword.to.fcs(current.project$fcs.files[[f]], 2, "MAPOP_pop_label")
                        }
                    }
                    
                    
                    write.enriched.FCS(current.project$fcs.files[[f]], paste0(env.var$tool.wd,"/Projects/",
                                                                           current.project$name,"/data/",
                                                                           t,
                                                                           ".fcs"))
                }
                progress$inc(1/length(current.project$fcs.files), detail = paste0("File ", f, " saved"))
            })
        }
        progress$set(message="Done",value=1)
        delay(500, progress$close())
        load.project()
    })
    
    observeEvent(input$t_1_1_load_project,#LOAD PROJECT
    {
        load.project()
    })
    
    observeEvent(input$t_1_1_remove_project,#REMOVE PROJECT
    {
        unlink(paste0(env.var$tool.wd,"/Projects/",current.project$name), force = T, recursive = T)
        remove.files()
        update.project.list()
        
    })
    
    observe(#LOAD FILES INFORMATION
    {
        input$t_1_1_load_project
        if(!is.null(current.project$name) && length(current.project$fcs.files)>0)
        {
            lapply(1:length(current.project$fcs.files), function(f)
            {
                fcs <- current.project$fcs.files[[f]]
                if(is.defined(fcs))
                {
                    idf <- names(current.project$fcs.files)[f]
                    
                    #POP COL LOADING---------------------
                    pop.col.sel <- 1:ncol(fcs@exprs)
                    names(pop.col.sel) <- lapply(1:ncol(fcs@exprs), function(j)
                    {
                        d <- fcs@description[[paste0("$P",j,"S")]]
                        if(is.null(d) || is.na(d) || d == "" || d == " " || d == "NA" || d == "<NA>")
                        {
                            d <- current.project$fcs.files.ui.colnames[[f]][[j]]
                        }
                        names(d) <- NULL
                        
                        return(unlist(d))
                    })
                    map.col.sel <- NULL
                    curr.map.label <- list()
                    if(is.defined(current.project$mapping.files[[f]]))
                    {
                        map.col.sel <- 1:ncol(current.project$mapping.files[[f]])
                        names(map.col.sel) <- colnames(current.project$mapping.files[[f]])
                        if(keyword.exists.FCS(fcs,"MAPOP_pop_label"))
                        {
                            curr.map.label <- as.numeric(get.keywords.with.keypart.FCS(fcs,"MAPOP_pop_label")[[1]][[1]])
                        }
                    }
                    curr.file.label <- list()
                    if(keyword.exists.FCS(fcs,"RF_pop_label"))
                    {
                        curr.file.label <- as.numeric(get.keywords.with.keypart.FCS(fcs,"RF_pop_label")[[1]][[1]])
                    }
                    
                    #UI CREATION------------------------
                    if(current.project$modified.files[[f]])
                    {
                        removeUI(paste0("#t_1_3_fr_",f))
                        insertUI("#t_1_3",
                                 "beforeEnd",
                                 fluidRow
                                 (
                                     style="margin-left:1.7vw",id=paste0("t_1_3_fr_",f),
                                     box
                                     (
                                         title=names(current.project$fcs.files)[f],collapsible=TRUE,width=10,collapsed=F,
                                         id=paste0("t_1_3_",f), 
                                         tabBox
                                         (
                                             side="right",width=12,
                                             tabPanel
                                             (
                                                 "Populations",id=paste0("t_1_3_",current.project$name,"_",f,"_1"),
                                                 fluidRow
                                                 (
                                                     id=paste0("t_1_3_",current.project$name,"_",f,"_1_fr"),
                                                     box
                                                     (
                                                         id=paste0("t_1_3_",current.project$name,"_",f,"_1_1"), width=3,
                                                         selectInput(paste0("t_1_3_",current.project$name,"_",f,"_pop_col"), "Population Column", choices=pop.col.sel, 
                                                                     selected=curr.file.label),
                                                         selectInput(paste0("t_1_3_",current.project$name,"_",f,"_lab_col"), "Labels Column (Mapping File)", choices=map.col.sel, 
                                                                     selected=curr.map.label)
                                                     ),
                                                     box
                                                     (
                                                         id=paste0("t_1_3_",current.project$name,"_",f,"_1_2"), width=9
                                                     )
                                                 )
                                             ),
                                             tabPanel
                                             (
                                                 "Previous Analyses", 
                                                 fluidRow
                                                 (
                                                     id=paste0("t_1_3_",current.project$name,"_",f,"_2")
                                                 )
                                             )
                                         )
                                     ),
                                     box
                                     (
                                         width=2,height="12vh",
                                         checkboxInput(paste0("t_1_3_",current.project$name,"_",f,"_cbox"), "Select", value = F),
                                         actionButton(paste0("t_1_3_",current.project$name,"_",f,"_mfile"), "Mapping File")
                                     )
                                 )
                        )
                    }
                    
                    #PREVIOUS ANALYSES LOADING---------------------------
                    if(current.project$modified.files[[f]])
                    {
                        previous.analyses <- FPH.retrieve.clusters.data.from.file(fcs)
                        prev.an.algo <- previous.analyses[[1]]
                        prev.an.markers <- previous.analyses[[2]]
                        prev.an.param <- previous.analyses[[3]]
                        if(!is.null(prev.an.algo))
                        {
                            lapply(1:length(prev.an.algo), function(k)
                            {
                                curr.algorithms <- prev.an.algo[[k]]
                                curr.parameters <- prev.an.param[[k]]
                                curr.markers <- prev.an.markers[[k]]

                                if(!is.null(curr.algorithms))
                                {
                                    run.choices <- 1:length(curr.algorithms)
                                    names(run.choices) <- curr.algorithms

                                    insertUI(paste0("#t_1_3_",current.project$name,"_",f,"_2"),
                                             "beforeEnd",
                                             box
                                             (
                                                 id=paste0("t_1_3_",current.project$name,"_",f,"_2_b_",k), width=12,
                                                 collapsible=T, title=names(prev.an.algo)[k],
                                                 selectInput(paste0("t_1_3_",current.project$name,"_",f,"_2_b_",k,"_run"),"Select analysis",choices = run.choices),
                                                 box
                                                 (
                                                     title = "Markers",id=paste0("t_1_3_",current.project$name,"_",f,"_2_b_",k,"_mark"),
                                                     div
                                                     (
                                                         id=paste0("t_1_3_",current.project$name,"_",f,"_2_b_",k,"_mark_content")
                                                     )
                                                 ),
                                                 box
                                                 (
                                                     title = "Parameters",id=paste0("t_1_3_",current.project$name,"_",f,"_2_b_",k,"_param"),
                                                     div
                                                     (
                                                         id=paste0("t_1_3_",current.project$name,"_",f,"_2_b_",k,"_param_content")
                                                     )
                                                 )
                                             )
                                    )
                                }
                            })
                        }
                    }


                    #LOAD MAPPING FILE-------------------------------------
                    if(current.project$modified.files[[f]])
                    {
                        if(is.defined(input[[paste0("t_1_3_",current.project$name,"_",f,"_mfile")]]))
                        {
                            observeEvent(input[[paste0("t_1_3_",isolate(current.project$name),"_",f,"_mfile")]],
                            {
                                 m <- matrix(nrow=1,ncol=2)
                                 m[1,1] = "Mapping File"
                                 m[1,2] = "*.csv"
                                 temp.files <- choose.files(filters = m,multi = F)

                                 if(length(temp.files) != 0)
                                 {
                                     lapply(temp.files, function(fname)
                                     {
                                         l <- length(fname)
                                         x <- NULL
                                         if(grepl("csv",fname))
                                         {
                                             x <- as.matrix(read.csv(fname))

                                         }
                                         current.project$mapping.files[[names(current.project$fcs.files)[f]]] <<- x
                                     })
                                 }
                                 shinyjs::delay(500,shinyjs::enable(paste0("t_1_3_",current.project$name,"_",f,"_mfile")))
                            })
                            current.project$modified.files[[f]] <<- FALSE
                        }
                    }
                }
            })
        }
    })
    
    observeEvent(input$t_1_1_create_project,#CREATE PROJECT
    {
        if(!is.null(input$t_1_1_name_project) && !is.na(input$t_1_1_name_project) && (input$t_1_1_name_project!="")
           && (input$t_1_1_name_project!=" "))
        {
            if(!dir.exists(paste0(env.var$tool.wd,"/Projects/",input$t_1_1_name_project)))
            {
                dir.create(paste0(env.var$tool.wd,"/Projects/",input$t_1_1_name_project))
            }
        }
        projects <- list.dirs(paste0(env.var$tool.wd,"/Projects/"), full.names = F, recursive = F)
        if(length(projects) > 0)
        {
            names(projects) <- projects
            updateSelectInput(session,"t_1_1_select_project", "Select Project", choices=projects, selected = input$t_1_1_name_project)
        }
        current.project$name <- input$t_1_1_name_project
        
        load.project(current.project$name)
    })

    observe(#LOAD PREVIOUS ANALYSES - CONTENT
    {
        if(!is.null(current.project$name) && length(current.project$fcs.files)>0)
        {
            lapply(1:length(current.project$fcs.files), function(f)
            {
                fcs <- current.project$fcs.files[[f]]
                if(is.defined(fcs))
                {
                    idf <- names(current.project$fcs.files)[f]

                    #POP COL LOADING---------------------
                    pop.col.sel <- 1:ncol(fcs@exprs)
                    names(pop.col.sel) <- lapply(1:ncol(fcs@exprs), function(j)
                    {
                        d <- fcs@description[[paste0("$P",j,"S")]]
                        if(is.null(d) || is.na(d) || d == "" || d == " " || d == "<NA>" || d == "NA")
                        {
                            d <- current.project$fcs.files.ui.colnames[[f]][[j]]
                        }
                        names(d) <- NULL

                        return(unlist(d))
                    })


                    #PREVIOUS ANALYSES LOADING---------------------------
                    previous.analyses <- FPH.retrieve.clusters.data.from.file(fcs)
                    prev.an.algo <- previous.analyses[[1]]
                    prev.an.markers <- previous.analyses[[2]]
                    prev.an.param <- previous.analyses[[3]]
                    if(!is.null(prev.an.algo))
                    {
                        lapply(1:length(prev.an.algo), function(k)
                        {
                            curr.algorithms <- prev.an.algo[[k]]
                            curr.parameters <- prev.an.param[[k]]
                            curr.markers <- prev.an.markers[[k]]

                            if(!is.null(curr.algorithms))
                            {
                                run.choices <- 1:length(curr.algorithms)
                                names(run.choices) <- curr.algorithms

                                if(is.defined(input[[paste0("t_1_3_",current.project$name,"_",f,"_2_b_",k,"_run")]]) &&
                                   input[[paste0("t_1_3_",current.project$name,"_",f,"_2_b_",k,"_run")]] != "" &&
                                   input[[paste0("t_1_3_",current.project$name,"_",f,"_2_b_",k,"_run")]] != " ")
                                {
                                    l <- as.numeric(input[[paste0("t_1_3_",current.project$name,"_",f,"_2_b_",k,"_run")]])

                                    removeUI(paste0("#t_1_3_",current.project$name,"_",f,"_2_b_",k,"_mark_content"))
                                    insertUI(paste0("#t_1_3_",current.project$name,"_",f,"_2_b_",k,"_mark"),
                                             "beforeEnd",
                                             div(id=paste0("t_1_3_",current.project$name,"_",f,"_2_b_",k,"_mark_content")))

                                    removeUI(paste0("#t_1_3_",current.project$name,"_",f,"_2_b_",k,"_param_content"))
                                    insertUI(paste0("#t_1_3_",current.project$name,"_",f,"_2_b_",k,"_param"),
                                             "beforeEnd",
                                             div(id=paste0("t_1_3_",current.project$name,"_",f,"_2_b_",k,"_param_content")))


                                    lapply(1:length(curr.markers[[l]]), function(m)
                                    {
                                        insertUI(paste0("#t_1_3_",current.project$name,"_",f,"_2_b_",k,"_mark_content"),
                                                 "beforeEnd",
                                                 h5(names(pop.col.sel)[[as.integer(curr.markers[[l]][m])]])
                                        )
                                    })

                                    if(!is.null(curr.parameters[[l]]))
                                    {
                                        lapply(1:length(curr.parameters[[l]]), function(m)
                                        {
                                            par.name <- strsplit(curr.parameters[[l]][m],"-")[[1]][1]
                                            par.val <- strsplit(curr.parameters[[l]][m],"-")[[1]][2]
                                            insertUI(paste0("#t_1_3_",current.project$name,"_",f,"_2_b_",k,"_param_content"),
                                                     "beforeEnd",
                                                     h5(paste0(par.name,": ",par.val))
                                            )
                                        })
                                    }
                                }
                            }
                        })
                    }
                }
            })
        }
    })

    observe(#UPDATE LABEL POP LIST
    {
        if(!is.null(current.project$name) && length(current.project$fcs.files)>0)
        {
            lapply(1:length(current.project$fcs.files), function(f)
            {
                    fcs <- current.project$fcs.files[[f]]

                    if(is.defined(fcs))
                    {
                        f.name <- names(current.project$fcs.files)[f]
                        if(is.defined(current.project$mapping.files))
                        {
                            if(is.defined(current.project$mapping.files[[f.name]]))
                            {
                                #LABEL COL LOADING
                                mapping.file <- current.project$mapping.files[[f.name]]
                                lab.col.sel <- 1:ncol(mapping.file)
                                names(lab.col.sel) <- colnames(mapping.file)
                                curr.file.label <- NULL
                                if(keyword.exists.FCS(fcs,"MAPOP_pop_label"))
                                {
                                    curr.file.label <- get.keywords.with.keypart.FCS(fcs,"MAPOP_pop_label")[[1]]
                                }
                                updateSelectInput(session,paste0("t_1_3_",current.project$name,"_",f,"_lab_col"), "Labels Column (Mapping File)", choices=lab.col.sel,
                                            selected=curr.file.label)
                            }
                        }
                    }
            })
        }
    })

    observe(#LOAD FILES CLUSTER
    {
        progress <- shiny::Progress$new()
        progress$set(message="RETRIEVING POPULATIONS LIST",value=0)
        on.exit(progress$close())
        if(length(current.project$fcs.files)>0)
        {
            lapply(1:length(current.project$fcs.files), function(f)
            {
                fcs <- current.project$fcs.files[[f]]
                if(is.defined(fcs))
                {
                    if(is.defined(input[[paste0("t_1_3_",current.project$name,"_",f,"_pop_col")]]) &&
                       input[[paste0("t_1_3_",current.project$name,"_",f,"_pop_col")]]!="")
                    {
                        if(length(unique(fcs@exprs[,as.integer(input[[paste0("t_1_3_",current.project$name,"_",f,"_pop_col")]])])) < 100)
                        {
                            f.name <- names(current.project$fcs.files)[f]
                            fcs.clusters <- FPH.get.file.clusters(fcs,as.integer(input[[paste0("t_1_3_",current.project$name,"_",f,"_pop_col")]]))
                            nmb.ev <- sum(sapply(1:length(fcs.clusters),function(cl){return(fcs.clusters[[cl]][[1]])}))
                            map.labels <- 1:length(fcs.clusters)
                            if(is.defined(current.project$mapping.files))
                            {
                                if(is.defined(input[[paste0("t_1_3_",current.project$name,"_",f,"_lab_col")]]) &&
                                   is.defined(input[[paste0("t_1_3_",current.project$name,"_",f,"_lab_col")]]) &&
                                   input[[paste0("t_1_3_",current.project$name,"_",f,"_lab_col")]]!="")
                                {
                                    t <- FPH.get.labels.from.mapping.file(current.project$mapping.files[[f.name]],
                                                                          input[[paste0("t_1_3_",current.project$name,"_",f,"_lab_col")]])
                                    map.labels[1:length(t)] <- t
                                }
                            }
    
                            removeUI(paste0("#t_1_3_",current.project$name,"_",f,"_1_2_fr"))
                            insertUI(paste0("#t_1_3_",current.project$name,"_",f,"_1_2"),
                                     "beforeEnd",
                                     fluidRow
                                     (
                                         id=paste0("t_1_3_",current.project$name,"_",f,"_1_2_fr"),style="margin-left:2%;max-height:30vh;overflow:auto"
                                     )
                            )
                            if(length(fcs.clusters)<50)
                            {
                                lapply(1:length(fcs.clusters), function(cl)
                                {
                                    insertUI(paste0("#t_1_3_",current.project$name,"_",f,"_1_2_fr"),
                                             "beforeEnd",
                                             textInput(paste0("t_1_3_",current.project$name,"_",f,"_1_pop_",cl),
                                                       paste0("Population ",cl," - Label:"),
                                                       value = map.labels[[cl]],
                                                       width = "70%")
                                    )
                                    progress$inc(1/length(fcs.clusters), detail=paste0("Population ", cl, " retrieved"))
                                })
                            }
                        }
                    }
                }
            })
        }
        progress$set(message="Done", value=1)
    })


    

    observe(#ACTIVATE ADD/RM/DL UI
    {
        if(!is.null(current.project$name))
        {
            shinyjs::enable("t_1_2_add")
            if(length(current.project$fcs.files)>0)
            {
                nmb.na.files <- 0
                lapply(1:length(current.project$fcs.files), function(f)
                {
                    fcs <- current.project$fcs.files[[f]]
                    
                    if(!is.defined(fcs))
                    {
                        nmb.na.files <<- nmb.na.files+1
                    }
                })
                
                if(nmb.na.files >= length(current.project$fcs.files))
                {
                    shinyjs::disable("t_1_2_rm")
                    shinyjs::disable("t_1_2_dl")
                }
                else
                {
                    shinyjs::enable("t_1_2_rm")
                    shinyjs::enable("t_1_2_dl")
                }
            }
            else
            {
                shinyjs::disable("t_1_2_rm")
                shinyjs::disable("t_1_2_dl")
            }
        }
        else
        {
            shinyjs::disable("t_1_2_add")
            shinyjs::disable("t_1_2_rm")
            shinyjs::disable("t_1_2_dl")
        }
    })
    
    observeEvent(input$t_1_2_add,#ADD FILES TO PROJECT
    {
        progress <- shiny::Progress$new()
        progress$set(message="ADDING FILES TO PROJECT",value=0)
        on.exit(progress$close())
        shinyjs::disable("t_1_2_add")
        
        if( !is.null(current.project$name) )
        {
            m <- matrix(nrow=1,ncol=2)
            m[1,1] = "FlowFrames"
            m[1,2] = "*.csv;*.fcs"
            temp.files <- choose.files(filters = m,multi = T)

            if(length(temp.files) != 0)
            {
                lapply(temp.files, function(f)
                {
                    l <- length(f)
                    x <- NULL
                    nx <- list()
                    if(grepl("csv",f))
                    {
                        x <- as.matrix(read.csv(f))
                        x <- flowFrame(x)
                        lapply(1:ncol(x@exprs), function(i)
                        {
                            d <- x@description[[paste0("$P",i,"S")]]
                            if(is.null(d) || is.na(d) || d == "" || d == " " || d == "NA" || d == "<NA>" || d == "'<NA>'")
                            {
                                d <- colnames(x)[i]
                            }
                            nx[[i]] <<- d
                        })
                    }
                    else
                    {
                        x <- read.FCS(f,emptyValue = FALSE)
                        nx <- list()
                        lapply(1:ncol(x@exprs), function(i)
                        {
                            d <- x@description[[paste0("$P",i,"S")]]
                            if(is.null(d) || is.na(d) || d == "" || d == " " || d == "NA" || d == "<NA>" || d == "'<NA>'")
                            {
                                d <- colnames(x)[i]
                            }
                            nx[[i]] <<- d
                        })
                    }
                    if( is.null(current.project$fcs.files) )
                    {
                        current.project$fcs.files <- list()
                        current.project$fcs.files.ui.colnames <- list()
                    }
                    current.project$fcs.files[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] <<- x
                    current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] <<- nx
                    current.project$modified.files[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] <- TRUE
                    current.project$mapping.files[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] <- NA

                    computed.values$purity.matrix.annot <- list()
                    computed.values$purity.matrix.clust <- list()
                    computed.values$FG.matrices.annot <- list()
                    computed.values$FG.matrices.clust <- list()
                    computed.values$prec.rec.matrices.annot <- list()
                    computed.values$prec.rec.matrices.clust <- list()
                    computed.values$pop.sizes <- list()
                    computed.values$clust.sizes <- list()
                    computed.values$annot.sizes <- list()
                    
                    progress$inc(1/length(temp.files), detail = paste0("File ", f, " added"))
                })
            }
        }
        
        shinyjs::delay(500, progress$set(message="DONE", value=1))
        shinyjs::enable("t_1_2_add")
    })

    observeEvent(input$t_1_2_rm,#REMOVE SELECTED FILES FROM PROJECT
    {
        if( !is.null(current.project$name) && length(current.project$fcs.files) >0 )
        {
            rm.files.ids <- c()
            for(i in 1:length(current.project$fcs.files))
            {
                if(is.defined(current.project$fcs.files[[i]]))
                {
                    if(input[[paste0("t_1_3_",current.project$name,"_",i,"_cbox")]])
                    {
                        rm.files.ids <- c(rm.files.ids,i)
                    }
                }
            }
            
            remove.files(rm.files.ids)
        }
    })

    output$t_1_2_dl <- downloadHandler(#DOWNLOAD SELECTED FILES
        filename = function()
        {
            paste0(current.project$name,"_SEL.zip")
        },
        content = function(file)
        {
            f.names <- c()
            if( !is.null(current.project$name) && length(current.project$fcs.files) >0 )
            {
                dl.files.ids <- c()
                lapply(1:length(current.project$fcs.files), function(i)
                {
                    if(input[[paste0("t_1_3_",current.project$name,"_",i,"_cbox")]])
                    {
                        dl.files.ids <<- c(dl.files.ids,i)
                    }
                })
                lapply(dl.files.ids, function(i)
                {
                    f.names <<- c(f.names,paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/",
                                       names(current.project$fcs.files)[i],".fcs"))
                    f.names <<- c(f.names,paste0(env.var$tool.wd,"/Projects/",current.project$name,"/data/",
                                                names(current.project$fcs.files)[i],".csv"))
                })
            }
            zip(zipfile = file, files = f.names)
        }
    )








    #======================================================================================================================
    #======================================================================================================================
    #==========================================SCORING=====================================================================
    #======================================================================================================================

    
    update.files.list <- function(current.section = NULL)
    {
        if(!is.null(current.section))
        {
            #FILE SELECTION-------------------------------------------------------------------
            f.list <- list()
            lapply(1:length(current.project$fcs.files), function(f.id)
            {
                if(is.defined(current.project$fcs.files[[f.id]]))
                {
                    f.name <- names(current.project$fcs.files)[f.id]
                    f.list[[f.name]] <<- NULL
                    f.list[[f.name]] <<- f.id
                }
            })
            if(is.defined(f.list))
            {
                updateSelectInput(session, paste0(current.section,"_fileSel"), choices=f.list, selected=f.list[[1]])
            }
        }
    }
    
    update.algorithms.list <-function(current.section = NULL)
    {
        fcs <- current.project$fcs.files[[as.integer(input[[paste0(current.section,"_fileSel")]])]]
        if(is.defined(fcs))
        {
            f.name <- names(current.project$fcs.files)[as.integer(input[[paste0(current.section,"_fileSel")]])]
            t <- names(computed.values$FG.matrices.annot[[as.numeric(input[[paste0(current.section,"_fileSel")]])]])
            algo.list <- lapply(1:length(t),function(i){return(i)})
            names(algo.list) <- t
            updateSelectInput(session, paste0(current.section,"_methodSel"), choices=algo.list, selected=algo.list[[1]])
        }
    }
    
    update.runs.list <- function(current.section = NULL)
    {
        if(!is.null(current.section))
        {
            fcs <- current.project$fcs.files[[as.integer(input[[paste0(current.section,"_fileSel")]])]]
            if(is.defined(fcs))
            {
                analyses.list <- FPH.retrieve.clusters.data.from.file(fcs)
                analyses.algorithms <- analyses.list[[1]]
                analyses.parameters <- analyses.list[[3]]
                current.algo <- analyses.algorithms[[as.integer(input[[paste0(current.section,"_methodSel")]])]]
                current.params <- analyses.parameters[[as.integer(input[[paste0(current.section,"_methodSel")]])]]
                
                available.runs <- 1:length(current.algo)
                for(current.algo.run.id in 1:length(current.algo))
                {
                    current.algo.run <- current.algo[[current.algo.run.id]]
                    tmp.run.name <- paste0(strsplit(current.algo.run,"__", fixed = T)[[1]][2],": ")
                    tmp.run.parameters <- extract.run.parameters(current.algo.run)
                    if(length(tmp.run.parameters)>0)
                    {
                        for(par.id in 1:length(tmp.run.parameters))
                        {
                            tmp.run.name <- paste0(tmp.run.name, names(tmp.run.parameters)[par.id], "=", tmp.run.parameters[[par.id]], ", ")
                        }
                    }
                    names(available.runs)[current.algo.run.id] <- tmp.run.name
                }
                
                if(is.defined(available.runs))
                {
                    updateSelectInput(session, paste0(current.section,"_runSel"), choices=available.runs, selected=available.runs[[1]])
                }
            }
        }
    }
    
    
    observe(#ACTIVATE RUN BUTTON
    {
        tmp <- 0
        map.tmp <- 0
        if(!is.null(current.project$name) && length(current.project$fcs.files)>0)
        {
            lapply(1:length(current.project$fcs.files), function(f.id)
            {
                if(is.defined(current.project$fcs.files[[f.id]]))
                {
                    tmp <<- tmp+1
                    if(is.defined(current.project$mapping.files[[f.id]]))
                    {
                        map.tmp <<- map.tmp+1
                    }
                }
            })
        }
        if(tmp>0 && map.tmp==tmp)
        {
            shinyjs::enable("t_3_2_run")
        }
        else
        {
            shinyjs::disable("t_3_2_run")
        }
    })

    observeEvent(input$t_3_2_run,#RUN SCORE COMPUTING
    {
        if( !is.null(current.project$name) && length(current.project$fcs.files) > 0)
        {
            showNotification("RUNNING ANALYSES",duration=NULL,type = "message",id="score_compute_message")
            #RETRIEVE CLUSTERS COLUMNS--------------------------------------------------------------------
            
            #Store reactive and global values
            current.project$ref.files.populations.col <- list()
            current.project$test.files.clusters.col <- list()
            lapply(1:length(current.project$fcs.files), function(f.id)
            {
                fcs <- current.project$fcs.files[[f.id]]
                if(is.defined(fcs))
                {
                    fcs.name <- names(current.project$fcs.files)[f.id]
                    analyses.list <- FPH.retrieve.clusters.data.from.file(fcs)
                        analyses.algorithms <- analyses.list[[1]]
                        analyses.markers <- analyses.list[[2]]
                        analyses.parameters <- analyses.list[[3]]
                        analyses.columns <- analyses.list[[4]]

                    if(is.null(current.project$ref.files.populations.col[[fcs.name]]))
                    {
                        current.project$ref.files.populations.col[[fcs.name]] <<- list()
                        current.project$test.files.clusters.col[[fcs.name]] <<- list()
                    }
                    if( length(analyses.algorithms)>0 )
                    {
                        lapply(1:length(analyses.algorithms), function(a)
                        {
                            algo.name <- names(analyses.algorithms)[a]
                            if(is.null(current.project$ref.files.populations.col[[fcs.name]][[algo.name]]))
                            {
                                current.project$ref.files.populations.col[[fcs.name]][[algo.name]] <<- list()
                                current.project$test.files.clusters.col[[fcs.name]][[algo.name]] <<- list()
                            }

                            if( length(analyses.algorithms[[a]])>0 )
                            {
                                lapply(1:length(analyses.algorithms[[a]]), function(a.run)
                                {
                                    current.project$ref.files.populations.col[[fcs.name]][[algo.name]] <<- c(current.project$ref.files.populations.col[[fcs.name]][[algo.name]],
                                                                                                            get.keywords.with.keypart.FCS(fcs,"RF_pop_label")[[1]])
                                    current.project$test.files.clusters.col[[fcs.name]][[algo.name]] <<- c(current.project$test.files.clusters.col[[fcs.name]][[algo.name]],
                                                                                                          analyses.columns[[a]][[a.run]])
                                })
                            }
                        })
                    }
                }
            })

            #COMPUTE MATRICES----------------------------------------------------------------------------
            computed.values$purity.matrix.annot <- list()
            computed.values$purity.matrix.clust <- list()
            computed.values$FG.matrices.annot <- list()
            computed.values$FG.matrices.clust <- list()
            computed.values$prec.rec.matrices.annot <- list()
            computed.values$prec.rec.matrices.clust <- list()
            computed.values$pop.sizes <-   list()
            computed.values$clust.sizes <-   list()
            computed.values$annot.sizes <-   list()
            
            tmp.fct <- function()
            {
                lapply(1:length(current.project$fcs.files), function(f.id)
                {
                    fcs <- current.project$fcs.files[[f.id]]
                    f.name <<- names(current.project$fcs.files)[f.id]
                    if(is.defined(fcs))
                    {
                        fcs.populations <- FPH.get.file.clusters(fcs,
                                                              as.numeric(current.project$ref.files.populations.col[[f.name]][[1]][[1]]))
                        pop.sizes <- sapply(fcs.populations, function(pop)
                        {
                            return(pop[[1]])
                        })
                        pop.sizes <- unlist(pop.sizes)/ sum(unlist(pop.sizes))
                        for( i in 1:length(pop.sizes) )
                        {
                            pop.ev <- as.integer(unlist(fcs.populations[[i]][[2]]))
                            pop.col <- as.numeric(current.project$ref.files.populations.col[[f.name]][[1]][[1]])
                            pop.real.ID <- as.integer(unique(fcs@exprs[pop.ev,pop.col]))
                            ta <- table(fcs@exprs[,pop.col])
                            pop.real.ID <- as.integer(which(as.integer(names(ta))==pop.real.ID)[[1]])
                            
                            names(pop.sizes)[i] <- current.project$mapping.files[[f.id]][pop.real.ID,2]
                        }
                        
                        
                        showNotification(paste("    File",f.id),id="score_compute_message",duration=NULL,type="message")
                        if(current.project$verbose)
                        {
                            print(paste("    File",f.id))
                        }
                        
                        computed.values$purity.matrix.clust[[f.name]] <<-   list()
                        computed.values$purity.matrix.annot[[f.name]] <<-   list()
                        #computed.values$prec.rec.matrices.clust[[f.name]] <<-   list()
                        computed.values$prec.rec.matrices.annot[[f.name]] <<-   list()
                        #computed.values$FG.matrices.clust[[f.name]] <<-   list()
                        computed.values$FG.matrices.annot[[f.name]] <<-   list()
                        computed.values$pop.sizes[[f.name]] <<-   list()
                        computed.values$clust.sizes[[f.name]] <<-   list()
                        computed.values$annot.sizes[[f.name]] <<-   list()
    
                        if(length(current.project$ref.files.populations.col[[f.name]])>0)
                        {
                            lapply(1:length(current.project$ref.files.populations.col[[f.name]]), function(current.algo)
                            {
                                algo.name <- names(current.project$ref.files.populations.col[[f.name]])[current.algo]
                                showNotification(paste("        Algo:",algo.name),id="score_compute_message_2",duration=NULL,type="message")
                                if(current.project$verbose)
                                {
                                    print(paste("        Algo:",algo.name))
                                }
                                computed.values$purity.matrix.clust[[f.name]][[algo.name]] <<-   list()
                                computed.values$purity.matrix.annot[[f.name]][[algo.name]] <<-   list()
                                #computed.values$prec.rec.matrices.clust[[f.name]][[algo.name]] <<-   list()
                                computed.values$prec.rec.matrices.annot[[f.name]][[algo.name]] <<-   list()
                                #computed.values$FG.matrices.clust[[f.name]][[algo.name]] <<-   list()
                                computed.values$FG.matrices.annot[[f.name]][[algo.name]] <<-   list()
                                computed.values$pop.sizes[[f.name]][[algo.name]] <<-   list()
                                computed.values$clust.sizes[[f.name]][[algo.name]] <<-   list()
                                computed.values$annot.sizes[[f.name]][[algo.name]] <<-   list()
    
                                if(length(current.project$ref.files.populations.col[[f.name]][[algo.name]])>0)
                                {
                                    tmp.curr.proj <- isolate(reactiveValuesToList(current.project))
                                    runs.list <- 1:length(tmp.curr.proj$ref.files.populations.col[[f.name]][[algo.name]])
                                    clust.col.list <- tmp.curr.proj$test.files.clusters.col[[f.name]][[algo.name]]
                                    
                                    file.size <- object.size(fcs)
                                    nmb.cl <- get.nmb.cores.max(file.size, 
                                                                available.cores = current.project$nmb.cores, 
                                                                x.cores = 0.5, x.ram = 0.4, correction.coef = 1.05)
                                    cl <- makeCluster(nmb.cl)
                                    registerDoSNOW(cl)
                                    
                                    notification.3.fct <- function(i)
                                    {
                                        print(paste0("            Run: ",i))
                                        showNotification(paste0("Run: ",i),id="score_compute_message_3",duration=NULL,type="message")
                                    }
                                    showNotification("Preparing computation",id="score_compute_message_3",duration=NULL,type="message")
                                    temp.out <- foreach(run.id=runs.list,clust.col=clust.col.list,
                                                        .options.snow = list(progress=notification.3.fct),
                                                        .packages=c("flowCore"),
                                                        .export = c("is.defined","runs.list","pop.sizes","fcs","fcs.populations",
                                                                    "FPH.get.file.clusters", "FPH.get.purity.matrix",
                                                                    "FPH.get.prec.rec.matrices", "FPH.compute.F.G.matrix",
                                                                    "FPH.annotate.clusters.to.fcs")) %dopar%
                                    {
                                        #run.name <- current.project$ref.files.populations.col[[f.name]][[algo.name]]
                                        
                                        fcs.clusters <-  FPH.get.file.clusters(fcs, as.numeric(clust.col))
                                        clust.sizes <- sapply(fcs.clusters, function(cl)
                                        {
                                            return(cl[[1]])
                                        })
                                        clust.sizes <- clust.sizes / sum(unlist(clust.sizes))
                                        mat.clust <- FPH.get.purity.matrix(fcs.populations, fcs.clusters)
                                        
                                        
                                        
                                        fcs.annotations.file <- FPH.annotate.clusters.to.fcs(fcs, mat.clust, as.numeric(clust.col))
                                        fcs.annotations <- FPH.get.file.clusters(fcs.annotations.file,as.numeric(ncol(fcs.annotations.file@exprs)))
                                        annot.sizes <- sapply(fcs.annotations, function(an)
                                        {
                                            return(an[[1]])
                                        })
                                        annot.sizes <<- unlist(annot.sizes) / sum(unlist(annot.sizes))
                                        mat.annot <- FPH.get.purity.matrix(fcs.populations, fcs.annotations)
                                        
                                        
                                        
                                        #prec.rec.clust <- FPH.get.prec.rec.matrices(fcs.populations, fcs.clusters)
                                        prec.rec.annot <- FPH.get.prec.rec.matrices(fcs.populations, fcs.annotations)
                                        
                                        
                                        
                                        #FG.clust <- FPH.compute.F.G.matrix(prec.rec.clust)
                                        FG.annot <- FPH.compute.F.G.matrix(prec.rec.annot)
                                        
                                        return(list(pop.sizes,
                                                    annot.sizes, mat.annot, prec.rec.annot, FG.annot,
                                                    clust.sizes, mat.clust#, prec.rec.clust, FG.clust
                                                    ))
                                    }
                                    
                                    for (current.run in unlist(runs.list) )
                                    {
                                        computed.values$pop.sizes[[f.name]][[algo.name]][[current.run]] <<- temp.out[[current.run]][[1]]
                                        computed.values$annot.sizes[[f.name]][[algo.name]][[current.run]] <<- temp.out[[current.run]][[2]]
                                        computed.values$purity.matrix.annot[[f.name]][[algo.name]][[current.run]] <<- temp.out[[current.run]][[3]]
                                        computed.values$prec.rec.matrices.annot[[f.name]][[algo.name]][[current.run]] <<- temp.out[[current.run]][[4]]
                                        computed.values$FG.matrices.annot[[f.name]][[algo.name]][[current.run]] <<- temp.out[[current.run]][[5]] 
                                        computed.values$clust.sizes[[f.name]][[algo.name]][[current.run]] <<- temp.out[[current.run]][[6]]
                                        computed.values$purity.matrix.clust[[f.name]][[algo.name]][[current.run]] <<- temp.out[[current.run]][[7]]
                                        # computed.values$prec.rec.matrices.clust[[f.name]][[algo.name]][[current.run]] <<- temp.out[[8]]
                                        # computed.values$FG.matrices.clust[[f.name]][[algo.name]][[current.run]] <<- temp.out[[9]] 
                                    }
                                    stopCluster(cl)
                                }
                            })
                        }
                    }
                })
            }
            showNotification("RUNNING ANALYSES",duration=NULL,type = "message",id="score_compute_message")
            if(current.project$verbose)
            {
                print("COMPUTING MATRICES:")
                tmp.time <- (as.numeric(microbenchmark(tmp.fct(), times=1, unit="ns")$time) - 1500000)*(10^(-9))
                print(tmp.time)
                print(paste("********MATRICES COMPUTING TIME********",tmp.time,"s"))
                print("==========================================================")
            }
            else
            {
                tmp.fct()
            }
            #GENERATE UI--------------------------------------------------------------------
            shinyjs::show("t_3_3_fr")
            clustering.algorithms$run.analysis <- TRUE
            removeNotification("score_compute_message")
            removeNotification("score_compute_message_2")
            removeNotification("score_compute_message_3")
        }
    })






    observe(#UPDATE ANNOTATIONS VIZ UI - RUN SELECTION
    {
        if(!is.null(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(length(computed.values$FG.matrices.annot)>0)
            {
                #RUN SELECTION-------------------------------------------------------------------
                if(is.defined(input[["t_3_3_1_fileSel"]]) && input[["t_3_3_1_fileSel"]]!="" && input[["t_3_3_1_fileSel"]]!=" " &&
                   is.defined(input[["t_3_3_1_methodSel"]]) && input[["t_3_3_1_methodSel"]]!="" && input[["t_3_3_1_methodSel"]]!=" ")
                {
                    update.runs.list("t_3_3_1")
                }
            } 
        }
    }) 

    observe(#UPDATE ANNOTATIONS VIZ UI - ALGORITHM AND MARKERS SELECTION
    {
        if(!is.null(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(length(computed.values$FG.matrices.annot)>0)
            {
                #MARKERS SELECTION-------------------------------------------------------------------
                if(is.defined(input[["t_3_3_1_fileSel"]]) && input[["t_3_3_1_fileSel"]]!="" && input[["t_3_3_1_fileSel"]]!=" ")
                {
                    update.algorithms.list("t_3_3_1")
                    fcs <- current.project$fcs.files[[as.integer(input[["t_3_3_1_fileSel"]])]]
                    if(is.defined(fcs))
                    {
                        markers.list <- 1:length(colnames(fcs))
                        names(markers.list) <- colnames(fcs)
                        updateSelectInput(session, paste0("t_3_3_1_m1"), choices=markers.list, selected=markers.list[[1]])
                        updateSelectInput(session, paste0("t_3_3_1_m2"), choices=markers.list, selected=markers.list[[2]])
                    }
                }
            }
        }
    })

    observe(#UPDATE ANNOTATIONS VIZ UI - FILE SELECTION
    {
        if(!is.null(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(length(computed.values$FG.matrices.annot)>0)
            {
                update.files.list("t_3_3_1")

            }
        }
    })

    observeEvent(input$t_3_3_1_plotButton,#UPDATE ANNOTATIONS VIZ UI - PLOT
    {
        if(!is.null(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(length(computed.values$FG.matrices.annot)>0)
            {
                if(is.defined(input[["t_3_3_1_fileSel"]]) && input[["t_3_3_1_fileSel"]]!="" && input[["t_3_3_1_fileSel"]]!=" ")
                {
                    fcs <- current.project$fcs.files[[as.integer(input[["t_3_3_1_fileSel"]])]]
                    if(is.defined(fcs))
                    {
                        f.name <- names(current.project$fcs.files)[as.integer(input[["t_3_3_1_fileSel"]])]

                        if(is.defined(input[["t_3_3_1_methodSel"]]) && input[["t_3_3_1_methodSel"]]!="" && input[["t_3_3_1_methodSel"]]!=" ")
                        {
                            fcs <- current.project$fcs.files[[as.integer(input[["t_3_3_1_fileSel"]])]]
                            f.name <- names(current.project$fcs.files)[as.integer(input[["t_3_3_1_fileSel"]])]
                            if(is.defined(input[["t_3_3_1_runSel"]]) && input[["t_3_3_1_runSel"]]!="" && input[["t_3_3_1_runSel"]]!=" ")
                            {
                                #POP AND ANNOT LIST--------------------------------------------------------
                                tmp <- computed.values$pop.sizes[[f.name]][[as.integer(input[["t_3_3_1_methodSel"]])]][[as.integer(input[["t_3_3_1_runSel"]])]]
                                pop.names <- 1:length(tmp)
                                names(pop.names) <- names(tmp)
                                updateSelectInput(session, "t_3_3_1_refPopSel", "Select Highlighted Populations", choices=pop.names, selected = pop.names)


                                mat <- computed.values$FG.matrices.annot[[f.name]][[as.integer(input[["t_3_3_1_methodSel"]])]][[as.integer(input[["t_3_3_1_runSel"]])]][[1]]
                                tmp <- FPH.map.test.to.ref(mat)
                                annot.names <- 1:length(tmp)
                                names(annot.names) <- paste0(names(pop.names)[tmp],"__ANNOTATED")
                                updateSelectInput(session, "t_3_3_1_testPopSel", "Select Highlighted Annotations", choices=annot.names, selected = annot.names)
                            }
                        }
                    }
                }
                #OUTPUT HANDLING-------------------------------------------------------------------
                output$t_3_3_1_refPlot <- renderImage(#UPDATE F SCORE DETAILS UI - LEFT PLOT
                {
                    outfile <- tempfile(fileext = ".jpg")
                    #PLOTS-------------------------------------------------------------
                    if(is.defined(input[["t_3_3_1_runSel"]]) && input[["t_3_3_1_runSel"]]!="" && input[["t_3_3_1_runSel"]]!=" " &&
                       is.defined(input[["t_3_3_1_refPopSel"]]) && input[["t_3_3_1_refPopSel"]]!="" && input[["t_3_3_1_refPopSel"]]!=" " )
                    {
                        fcs <- current.project$fcs.files[[as.integer(input[["t_3_3_1_fileSel"]])]]
                        if(is.defined(fcs))
                        {
                            f.name <- names(current.project$fcs.files)[as.integer(input[["t_3_3_1_fileSel"]])]
                            mat <- fcs@exprs
                            #RETRIEVE POP NAMES AND SIZES-------------------------------------------------------------
                            pop.names <- computed.values$pop.sizes[[f.name]][[as.integer(input[["t_3_3_1_methodSel"]])]][[as.integer(input[["t_3_3_1_runSel"]])]]
                            pop.names <- names(pop.names)
                            pop.col <- current.project$ref.files.populations.col[[f.name]][[as.integer(input[["t_3_3_1_methodSel"]])]][[as.integer(input[["t_3_3_1_runSel"]])]]
                            pop.col <- as.integer(pop.col)

                            highlighted.pop <- as.integer(as.list(input[["t_3_3_1_refPopSel"]]))
                            names(highlighted.pop) <- pop.names[unlist(highlighted.pop)]
                            selected.pop <- lapply(highlighted.pop, function(p)
                            {
                                return(as.numeric(unlist(which(mat[,pop.col]==p))))
                            })


                            #PLOT---------------------------------------------------------------------------------------------------
                            jpeg(outfile)
                            plot.selected.clusters( mat, selected.pop, c(as.integer(input[["t_3_3_1_m1"]]),as.integer(input[["t_3_3_1_m2"]])) )
                            dev.off()
                        }
                    }
                    list(src=outfile)
                })
                output$t_3_3_1_testPlot <- renderImage(#UPDATE F SCORE DETAILS UI - RIGHT PLOT
                {
                    outfile <- tempfile(fileext = ".jpg")
                    #PLOTS-------------------------------------------------------------
                    if(is.defined(input[["t_3_3_1_runSel"]]) && input[["t_3_3_1_runSel"]]!="" && input[["t_3_3_1_runSel"]]!=" " &&
                       is.defined(input[["t_3_3_1_refPopSel"]]) && input[["t_3_3_1_refPopSel"]]!="" && input[["t_3_3_1_refPopSel"]]!=" " )
                    {
                        fcs <- current.project$fcs.files[[as.integer(input[["t_3_3_1_fileSel"]])]]
                        if(is.defined(fcs))
                        {
                            f.name <- names(current.project$fcs.files)[as.integer(input[["t_3_3_1_fileSel"]])]
                            mat <- fcs@exprs
                            #RETRIEVE POP NAMES AND SIZES-------------------------------------------------------------
                            pop.names <- computed.values$pop.sizes[[f.name]][[as.integer(input[["t_3_3_1_methodSel"]])]][[as.integer(input[["t_3_3_1_runSel"]])]]
                            pop.names <- names(pop.names)

                            clust.col <- current.project$test.files.clusters.col[[f.name]][[as.integer(input[["t_3_3_1_methodSel"]])]][[as.integer(input[["t_3_3_1_runSel"]])]]
                            clust.col <- as.integer(clust.col)
                            pur.mat <- computed.values$purity.matrix.clust[[f.name]][[as.integer(input[["t_3_3_1_methodSel"]])]][[as.integer(input[["t_3_3_1_runSel"]])]]
                            fcs.annot <- FPH.annotate.clusters.to.fcs(fcs, pur.mat, clust.col)

                            highlighted.annot <- as.integer(as.list(input[["t_3_3_1_testPopSel"]]))
                            names(highlighted.annot) <- paste0(pop.names[unlist(highlighted.annot)],"__AN")

                            mat <- computed.values$FG.matrices.annot[[f.name]][[as.integer(input[["t_3_3_1_methodSel"]])]][[as.integer(input[["t_3_3_1_runSel"]])]][[1]]
                            tmp <- FPH.map.test.to.ref(mat)
                            fcs.annot.cl <- FPH.get.file.clusters(fcs.annot,as.numeric(ncol(fcs.annot@exprs)))

                            selected.annot <- lapply(highlighted.annot, function(p)
                            {
                                return(as.numeric(unlist(fcs.annot.cl[[as.numeric(p)]][[2]])))
                            })
                            names(selected.annot) <- paste0(pop.names[as.numeric(tmp[highlighted.annot])],"__AN")

                            #PLOT---------------------------------------------------------------------------------------------------
                            jpeg(outfile)
                            plot.selected.clusters( fcs.annot@exprs, selected.annot, c(as.integer(input[["t_3_3_1_m1"]]),as.integer(input[["t_3_3_1_m2"]])) )
                            dev.off()
                        }
                    }
                    list(src=outfile)
                })
            }
        }
    })

    observe(#CLEAR UI
    {
        if(!clustering.algorithms$run.analysis)
        {
            updateSelectInput(session, "t_3_3_1_fileSel", "Select file", choices = list(), selected = NULL)
            updateSelectInput(session, "t_3_3_1_methodSel", "Select Algorithm", choices = list(), selected = NULL)
            updateSelectInput(session, "t_3_3_1_runSel", "Select Run", choices = list(), selected = NULL)
            updateSelectInput(session, "t_3_3_1_m1", "Marker 1", choices = list(), selected = NULL)
            updateSelectInput(session, "t_3_3_1_m2", "Marker 2", choices = list(), selected = NULL)
            updateSelectInput(session, "t_3_3_1_testPopSel", "Select Annotated Groups", choices=list())
            updateSelectInput(session, "t_3_3_1_refPopSel", "Select Populations", choices=list())
        }
    })






    generate.ordered.table <- function()
    {
        fcs <- current.project$fcs.files[[as.integer(input[["t_3_3_2_fileSel"]])]]
        if(is.defined(fcs))
        {
            f.name <- names(current.project$fcs.files)[as.integer(input[["t_3_3_2_fileSel"]])]
            analyses.list <- FPH.retrieve.clusters.data.from.file(fcs)
            analyses.algorithms <- analyses.list[[1]]
            analyses.markers <- analyses.list[[2]]
            analyses.parameters <- analyses.list[[3]]
            analyses.columns <- analyses.list[[4]]
            
            current.algo <- analyses.algorithms[[as.integer(input[["t_3_3_2_methodSel"]])]]
            current.params <- analyses.parameters[[as.integer(input[["t_3_3_2_methodSel"]])]]
            if(length(current.params)>0)
            {
                lapply(1:length(current.params), function(run)
                {
                    if(length(current.params[[run]])>0)
                    {
                        lapply(1:length(current.params[[run]]), function(p)
                        {
                            val <- strsplit(current.params[[run]][[p]],"-")[[1]][[2]]
                            val.name <- strsplit(current.params[[run]][[p]],"-")[[1]][[1]]
                            current.params[[run]][[p]] <<- val
                            names(current.params[[run]])[p] <<- val.name
                        })
                    }
                })
            }
            
            current.markers <- analyses.markers[[as.integer(input[["t_3_3_2_methodSel"]])]]
            mark.names <- current.project$fcs.files.ui.colnames[[as.integer(input[["t_3_3_2_fileSel"]])]]
            mark.param.concatenated.list <- lapply(1:max(length(current.params),length(current.markers)), function(run)
            {
                t <- list()
                if(is.defined(current.params[[run]]))
                {
                    t <- c(t,unlist(current.params[[run]]))
                    names(t) <- names(current.params[[run]])
                }
                if(is.defined(current.markers[[run]]))
                {
                    mark.run <- as.integer(unlist(current.markers[[run]]))
                    mark.run <- mark.run[order(mark.run)]
                    temp.val <- ""
                    lapply(1:length(mark.run), function(m)
                    {
                        temp.val <<- paste0(temp.val,mark.run[[m]])
                        return(temp.val)
                    })
                    t <- c(t,temp.val)
                    names(t) <- c(names(current.params[[run]]),"MARKERS")
                }
                
                return(t)
            })
            
            all.params.names <- get.params.list(mark.param.concatenated.list)
            available.params <- sapply(1:length(all.params.names), function(k)
            {
                ordered.table <- get.ordered.table.from.runs(mark.param.concatenated.list, all.params.names[[k]])
                tmp <- get.IDs.identical.sub.params(ordered.table)
                t.sizes <- sapply(tmp, function(run){return(length(run))})
                
                pot.run <- which(t.sizes==max(t.sizes))[[1]]
                if(nrow(as.matrix(ordered.table[unlist(tmp[[pot.run]]),]))>0)
                {
                    return(k)
                }
                else
                {
                    return(NULL)
                }
                
            })
            
            ordered.table <- get.ordered.table.from.runs(mark.param.concatenated.list, all.params.names[[as.integer(input[["t_3_3_2_paramSel"]])]])
            if(ncol(ordered.table)<2)
            {
                ordered.table <- t(ordered.table)
            }
            
            tmp.runs.set <- get.IDs.identical.sub.params(ordered.table)
            runs.set.values <- list()
            for(k in 1:length(tmp.runs.set))
            {
                runs.set.values[[k]] <- ""
                for(l in 3:ncol(ordered.table))
                {
                    runs.set.values[[k]] <- paste0(runs.set.values[[k]], 
                                                   colnames(ordered.table)[l], "=", ordered.table[tmp.runs.set[[k]][[1]],l],", ")
                }
            }
            null.set <- length(tmp.runs.set)+1
            for(tmp.set.id in 1:length(tmp.runs.set))
            {
                if(ordered.table[as.numeric(tmp.runs.set[[tmp.set.id]][[1]]),1] == "NULL")
                {
                    null.set <- tmp.set.id
                }
            }
            computed.values$fixed.parameters.ids <<- tmp.runs.set[-null.set]
            computed.values$fixed.parameters.values <<- runs.set.values[-null.set]
            
            pot.run <- 0
            if(is.defined(input$t_3_3_2_FixedparamSel))
            {
                pot.run <- as.integer(input$t_3_3_2_FixedparamSel)
            }
            else
            {
                t.sizes <- sapply(tmp.runs.set, function(run){return(length(run))})
                pot.run <- which(t.sizes==max(t.sizes))[[1]]
            }
            
            ordered.table <- as.matrix(ordered.table[unlist(computed.values$fixed.parameters.ids[[pot.run]]),])
            if(ncol(ordered.table)<2)
            {
                ordered.table <- t(ordered.table)
            }
            
            
            #REMOVE IDENTICAL RUNS
            reduced.runs.list <- get.IDs.identical.main.param(ordered.table)
            for(l in 1:length(reduced.runs.list))
            {
                reduced.runs.list[[l]] <- reduced.runs.list[[l]][[1]]
            }
            ordered.table <- as.matrix(ordered.table[unlist(as.integer(reduced.runs.list)),])
            if(ncol(ordered.table)<2)
            {
                ordered.table <- t(ordered.table)
            }
            
            computed.values$ordered.table <<- ordered.table
        }
    }
    
    observeEvent(input$t_3_3_2_plotButton,#UPDATE F SCORE DETAILS UI - PLOT
    {
        if(!is.null(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(length(computed.values$FG.matrices.annot)>0)
            {
                output[["t_3_3_2_paramSelPlot"]] <- renderImage(#UPDATE F SCORE DETAILS UI - LEFT PLOT
                {
                    outfile <- tempfile(fileext = ".jpg")
                    if(is.defined(input[["t_3_3_2_paramSel"]]) && input[["t_3_3_2_paramSel"]]!="" && input[["t_3_3_2_paramSel"]]!=" ")
                    {
                        #COMPUTE PLOTS POINTS-------------------------------------------------------------
                        generate.ordered.table()
                        
                        runs.columns <- as.integer(computed.values$ordered.table[,1])
                        list.pop.points <- list()
                        list.pop.points.xval <- list()
                        pop.ids <- list()
                        for(i in 1:length(runs.columns))
                        {
                            mat <- computed.values$FG.matrices.annot[[f.name]][[as.integer(input[["t_3_3_2_methodSel"]])]][[runs.columns[i]]][[1]]
                            pop.names <- names(computed.values$pop.sizes[[f.name]][[as.integer(input[["t_3_3_2_methodSel"]])]][[1]])
                            
                            tmp <- FPH.map.test.to.ref(mat)
                            pop.ids[[i]] <- paste0(pop.names[tmp],"__AN")
                            for(j in 1:nrow(mat))
                            {
                                F.max <- max(mat[j,])
                                if(!(pop.ids[[i]][[j]]%in%names(list.pop.points)))
                                {
                                    list.pop.points[[pop.ids[[i]][[j]]]] <- NULL
                                    list.pop.points.xval[[pop.ids[[i]][[j]]]] <- NULL
                                }
                                
                                list.pop.points[[pop.ids[[i]][[j]]]] <- c(unlist(list.pop.points[[pop.ids[[i]][[j]]]]), F.max)
                                list.pop.points.xval[[pop.ids[[i]][[j]]]] <- c(unlist(list.pop.points.xval[[pop.ids[[i]][[j]]]]), 
                                                                               as.numeric(computed.values$ordered.table[i,2]))
                            }
                            
                        }
                        
                        #PLOT F SCORES---------------------------------------------------------------------------------------------------
                        if(length(list.pop.points)>1)
                        {
                            max.height <- length(list.pop.points)+0.5
                            ordered.table <- as.matrix(computed.values$ordered.table)
                            values.range <- c(min(as.numeric(ordered.table[,2]))-1, max(as.numeric(ordered.table[,2]))+1)
                            jpeg(outfile, width=640, height=480)
                            draw.cumulated.filled.plots(list.pop.points, list.pop.points.xval, max.height, x.values.range = values.range,
                                                        y.lab="Cumulated F-Score", x.lab = colnames(ordered.table)[2])
                            dev.off()
                        }
                    }
                    list(src=outfile)
                })

                output[["t_3_3_2_runSelPlot"]] <- renderImage(#UPDATE F SCORE DETAILS UI - RIGHT PLOT
                {
                    outfile <- tempfile(fileext = ".jpg")
                    if(is.defined(input[["t_3_3_2_runSel"]]) && input[["t_3_3_2_runSel"]]!="" && input[["t_3_3_2_runSel"]]!=" ")
                    {
                        fcs <- current.project$fcs.files[[as.integer(input[["t_3_3_2_fileSel"]])]]
                        if(is.defined(fcs))
                        {
                            f.name <- names(current.project$fcs.files)[as.integer(input[["t_3_3_2_fileSel"]])]
                            #RETRIEVE POP NAMES AND SIZES-------------------------------------------------------------
                            F.mat <- computed.values$FG.matrices.annot[[f.name]][[as.integer(input[["t_3_3_2_methodSel"]])]][[as.integer(input[["t_3_3_2_runSel"]])]][[1]]
                            pop.sizes <- computed.values$pop.sizes[[f.name]][[as.integer(input[["t_3_3_2_methodSel"]])]][[as.integer(input[["t_3_3_2_runSel"]])]]
                            pop.names <- names(pop.sizes)
                            
                            analyses.list <- FPH.retrieve.clusters.data.from.file(fcs)
                            analyses.run <- analyses.list[[1]][[as.integer(input[["t_3_3_2_methodSel"]])]][[as.integer(input[["t_3_3_2_runSel"]])]]
                            
                            tmp.run.name <- paste0(strsplit(analyses.run,"__", fixed = T)[[1]][2],": ")
                            tmp.run.parameters <- extract.run.parameters(analyses.run)
                            if(length(tmp.run.parameters)>0)
                            {
                                for(par.id in 1:(length(tmp.run.parameters)-1))
                                {
                                    tmp.run.name <- paste0(tmp.run.name, names(tmp.run.parameters)[par.id], "=", tmp.run.parameters[[par.id]], ", ")
                                }
                            }
                            tmp.run.name <- paste0(tmp.run.name, "\n",names(tmp.run.parameters)[length(tmp.run.parameters)], "=", 
                                                   tmp.run.parameters[[length(tmp.run.parameters)]])

                            #PLOT F SCORES---------------------------------------------------------------------------------------------------
                            jpeg(outfile, width=720, height=480)
                            draw.F.score.barplot(F.mat, pop.names, pop.sizes, 
                                                 plot.title = tmp.run.name)
                            dev.off()
                        }
                    }
                    list(src=outfile)
                })

            }
        }
    })

    observe(#UPDATE F SCORE DETAILS UI - FIXED PARAMETERS
    {
        if(!is.null(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(!is.null(computed.values$fixed.parameters.ids))
            {
                choices.list <- 1:length(computed.values$fixed.parameters.ids)
                names(choices.list) <- computed.values$fixed.parameters.values
                
                t.sizes <- sapply(computed.values$fixed.parameters.ids, function(run){return(length(run))})
                pot.run <- which(t.sizes==max(t.sizes))[[1]]
                
                updateSelectInput(session, "t_3_3_2_FixedparamSel", choices=choices.list, selected=pot.run)
            }
        }
    })
    
    observe(#UPDATE F SCORE DETAILS UI - RUN SELECTION
    {
        if(!is.null(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(length(computed.values$FG.matrices.annot)>0)
            {
                if(is.defined(input[["t_3_3_2_fileSel"]]) && input[["t_3_3_2_fileSel"]]!="" && input[["t_3_3_2_fileSel"]]!=" ")
                {
                    if(is.defined(input[["t_3_3_2_methodSel"]]) && input[["t_3_3_2_methodSel"]]!="" && input[["t_3_3_2_methodSel"]]!=" ")
                    {
                        update.runs.list("t_3_3_2")
                    }
                }
            }
        }
    })

    observe(#UPDATE F SCORE DETAILS UI - PARAMETER SELECTION
    {
        if(!is.null(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(length(computed.values$FG.matrices.annot)>0)
            {
                if(is.defined(input[["t_3_3_2_fileSel"]]) && input[["t_3_3_2_fileSel"]]!="" && input[["t_3_3_2_fileSel"]]!=" ")
                {
                    fcs <- current.project$fcs.files[[as.integer(input[["t_3_3_2_fileSel"]])]]
                    if(is.defined(fcs))
                    {
                        if(is.defined(input[["t_3_3_2_methodSel"]]) && input[["t_3_3_2_methodSel"]]!="" && input[["t_3_3_2_methodSel"]]!=" ")
                        {
                            f.name <- names(current.project$fcs.files)[as.integer(input[["t_3_3_2_fileSel"]])]
                            analyses.list <- FPH.retrieve.clusters.data.from.file(fcs)
                            analyses.algorithms <- analyses.list[[1]]
                            analyses.markers <- analyses.list[[2]]
                            analyses.parameters <- analyses.list[[3]]
                            analyses.columns <- analyses.list[[4]]

                            current.algo <- analyses.algorithms[[as.integer(input[["t_3_3_2_methodSel"]])]]
                            current.params <- analyses.parameters[[as.integer(input[["t_3_3_2_methodSel"]])]]
                            if(length(current.params)>0)
                            {
                                lapply(1:length(current.params), function(run)
                                {
                                    if(length(current.params[[run]])>0)
                                    {
                                        lapply(1:length(current.params[[run]]), function(p)
                                        {
                                            val <- strsplit(current.params[[run]][[p]],"-")[[1]][[2]]
                                            val.name <- strsplit(current.params[[run]][[p]],"-")[[1]][[1]]
                                            current.params[[run]][[p]] <<- val
                                            names(current.params[[run]])[p] <<- val.name
                                        })
                                    }
                                })
                            }

                            current.markers <- analyses.markers[[as.integer(input[["t_3_3_2_methodSel"]])]]
                            mark.names <- current.project$fcs.files.ui.colnames[[as.integer(input[["t_3_3_2_fileSel"]])]]
                            mark.param.concatenated.list <- lapply(1:max(length(current.params),length(current.markers)), function(run)
                            {
                                t <- list()
                                if(is.defined(current.params[[run]]))
                                {
                                    t <- c(t,unlist(current.params[[run]]))
                                    names(t) <- names(current.params[[run]])
                                }
                                if(is.defined(current.markers[[run]]))
                                {
                                    mark.run <- as.integer(unlist(current.markers[[run]]))
                                    mark.run <- mark.run[order(mark.run)]
                                    temp.val <- ""
                                    lapply(1:length(mark.run), function(m)
                                    {
                                        temp.val <<- paste0(temp.val,mark.run[[m]])
                                        return(temp.val)
                                    })
                                    t <- c(t,temp.val)
                                    names(t) <- c(names(current.params[[run]]),"MARKERS")
                                }

                                return(t)
                            })

                            all.params.names <- get.params.list(mark.param.concatenated.list)
                            if(length(all.params.names)>0)
                            {
                                available.params <- sapply(1:length(all.params.names), function(k)
                                {
                                    ordered.table <- get.ordered.table.from.runs(mark.param.concatenated.list, all.params.names[[k]])
                                    tmp <- get.IDs.identical.sub.params(ordered.table)
                                    t.sizes <- sapply(tmp, function(run){return(length(run))})
    
                                    pot.run <- which(t.sizes==max(t.sizes))[[1]]
                                    if(nrow(as.matrix(ordered.table[unlist(tmp[[pot.run]]),]))>0)
                                    {
                                        return(k)
                                    }
                                    else
                                    {
                                        return(NULL)
                                    }
    
                                })
                                if(is.defined(available.params))
                                {
                                    names(available.params) <- all.params.names[available.params]
                                    updateSelectInput(session, "t_3_3_2_paramSel", choices=available.params, selected=available.params[[1]])
                                }
                            }
                        }
                    }
                }
            }
        }
    })

    observe(#UPDATE F SCORE DETAILS UI - ALGORITHM SELECTION
    {
        if(!is.null(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(length(computed.values$FG.matrices.annot)>0)
            {
                if(is.defined(input[["t_3_3_2_fileSel"]]) && input[["t_3_3_2_fileSel"]]!="" && input[["t_3_3_2_fileSel"]]!=" ")
                {
                    update.algorithms.list("t_3_3_2")
                }
            }
        }
    })

    observe(#UPDATE F SCORE DETAILS UI - FILE SELECTION
    {
        if(!is.null(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(length(computed.values$FG.matrices.annot)>0)
            {
                update.files.list("t_3_3_2")
            }
        }
    })

    observe(#CLEAR UI
    {
        if(!clustering.algorithms$run.analysis)
        {
            updateSelectInput(session, "t_3_3_2_fileSel", "Select file", choices = list(), selected = NULL)
            updateSelectInput(session, "t_3_3_2_methodSel", "Select Algorithm", choices = list(), selected = NULL)
            updateSelectInput(session, "t_3_3_2_runSel", "Select Run", choices = list(), selected = NULL)
            updateSelectInput(session, "t_3_3_2_paramSel", "Select Parameter", choices=list())
        }
    })





    observeEvent(input$t_3_3_3_plotButton,#UPDATE SUMMARIZING TABLE UI
    {
        if(!is.null(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(length(computed.values$FG.matrices.annot)>0)
            {
                algo.list <- NULL
                runs.list <- list()
                lapply(1:length(current.project$fcs.files), function(f.id)
                {
                    if(is.defined(current.project$fcs.files[[f.id]]))
                    {
                        fcs <- current.project$fcs.files[[f.id]]
                        f.name <- names(current.project$fcs.files)[f.id]
                        clusters.data <- FPH.retrieve.clusters.data.from.file(fcs)
                        used.algo <- clusters.data[[1]]
                        used.markers<- clusters.data[[2]]
                        used.param <- clusters.data[[3]]
                        lapply(1:length(used.algo), function(alg.id)
                        {
                            alg <- used.algo[[alg.id]]
                            alg.name <- names(used.algo)[alg.id]
                            if(!(alg.name%in%algo.list))
                            {
                                algo.list <<- c(algo.list, alg.name)
                                runs.list[[alg.name]] <<- list()
                            }
                            for(run.id in 1:length(alg))
                            {
                                run <- ""
                                if(is.defined(used.markers[[alg.id]][[run.id]]))
                                {
                                    run <- "Markers: "
                                    lapply(1:length(used.markers[[alg.id]][[run.id]]), function(mark.id)
                                    {
                                        run <<- paste0(run,used.markers[[alg.id]][[run.id]][[mark.id]],",")
                                    })
                                    run <- paste0(run,"    ")
                                }
                                if(is.defined(used.param[[alg.id]][[run.id]]))
                                {
                                    run <- paste0(run, "Parameters: ")
                                    lapply(1:length(used.param[[alg.id]][[run.id]]), function(mark.id)
                                    {
                                        run <<- paste0(run,used.param[[alg.id]][[run.id]][[mark.id]],",")
                                    })
                                }
                                if(!(run%in%runs.list[[alg.name]]))
                                {
                                    runs.list[[alg.name]] <<- c(runs.list[[alg.name]], run)
                                }
                            }
                        })
                    }
                })
                mat.ncol <- 1+sum(unlist(sapply(runs.list, function(alg){return(length(alg))})))
                summ.table <- matrix(NA,ncol=mat.ncol,nrow=length(Filter(Negate(is.na), current.project$fcs.files)))
                if(is.defined(runs.list))
                {
                    colnames(summ.table) <- c("FILE",unlist(sapply(1:length(runs.list), function(alg.id)
                    {
                        return(unlist(runs.list[[alg.id]]))
                    })))
                }


                lapply(1:length(current.project$fcs.files), function(f.id)
                {
                    fcs <- current.project$fcs.files[[f.id]]
                    if(is.defined(fcs))
                    {
                        f.name <- names(current.project$fcs.files)[f.id]
                        summ.table[f.id,1] <<- f.name
                        clusters.data <- FPH.retrieve.clusters.data.from.file(fcs)
                        used.algo <- clusters.data[[1]]
                        used.markers<- clusters.data[[2]]
                        used.param <- clusters.data[[3]]

                        lapply(1:length(algo.list), function(alg.id)
                        {
                            alg <- algo.list[[alg.id]]
                            lapply(1:length(runs.list[[alg]]), function(run.id)
                            {
                                run <- ""
                                if(is.defined(used.markers[[alg.id]][[run.id]]))
                                {
                                    run <- "Markers: "
                                    lapply(1:length(used.markers[[alg.id]][[run.id]]), function(mark.id)
                                    {
                                        run <<- paste0(run,used.markers[[alg.id]][[run.id]][[mark.id]],",")
                                    })
                                    run <- paste0(run,"    ")
                                }
                                if(is.defined(used.param[[alg.id]][[run.id]]))
                                {
                                    run <- paste0(run,"Parameters: ")
                                    lapply(1:length(used.param[[alg.id]][[run.id]]), function(mark.id)
                                    {
                                        run <<- paste0(run,used.param[[alg.id]][[run.id]][[mark.id]],",")
                                    })
                                }
                                tmp.mat <- computed.values$FG.matrices.annot[[f.name]][[alg.id]][[run.id]][[1]]
                                summ.table[f.id,run] <<- trunc(sum(as.numeric(unlist(sapply(1:nrow(tmp.mat), function(r){return(max(tmp.mat[r,]))}))))*1000)/1000
                            })
                        })

                        computed.values$summary.table <<- summ.table
                        output[["t_3_3_3_table"]] <- renderTable(computed.values$summary.table)
                    }
                })
            }
        }
    })





    observe(#UPDATE ADVANCED OPTIONS - RUN SELECTION
    {
        if(is.defined(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(is.defined(current.project$fcs.files) && length(computed.values$FG.matrices.annot)>0)
            {
                if(is.defined(input[["t_3_3_4_fileSel"]]) && input[["t_3_3_4_fileSel"]]!="" && input[["t_3_3_4_fileSel"]]!=" " &&
                   is.defined(input[["t_3_3_4_methodSel"]]) && input[["t_3_3_4_methodSel"]]!="" && input[["t_3_3_4_methodSel"]]!=" ")
                {
                    update.runs.list("t_3_3_4")
                }

            }
        }
    })

    observe(#UPDATE ADVANCED OPTIONS - ALGORITHM SELECTION
    {
        if(is.defined(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(is.defined(current.project$fcs.files) && length(computed.values$FG.matrices.annot)>0)
            {
                if(is.defined(input[["t_3_3_4_fileSel"]]) && input[["t_3_3_4_fileSel"]]!="" && input[["t_3_3_4_fileSel"]]!=" ")
                {
                    update.algorithms.list("t_3_3_4")
                }

            }
        }
    })

    observe(#UPDATE ADVANCED OPTIONS - FILE SELECTION
    {
        if(is.defined(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(is.defined(current.project$fcs.files) && length(computed.values$FG.matrices.annot)>0)
            {
                update.files.list("t_3_3_4")
            }
        }
    })

    observeEvent(input$t_3_3_4_plotButton,#UPDATE ADVANCED OPTIONS  - PLOT
    {
        if(is.defined(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(is.defined(current.project$fcs.files) && length(computed.values$FG.matrices.annot)>0)
            {
                #OUTPUT PLOT----------------------------------------------------------------------

                if(is.defined(isolate(input[["t_3_3_4_fileSel"]])) && isolate(input[["t_3_3_4_fileSel"]])!="" && isolate(input[["t_3_3_4_fileSel"]])!=" " &&
                   is.defined(isolate(input[["t_3_3_4_methodSel"]])) && isolate(input[["t_3_3_4_methodSel"]])!="" && isolate(input[["t_3_3_4_methodSel"]])!=" " &&
                   is.defined(isolate(input[["t_3_3_4_runSel"]])) && isolate(input[["t_3_3_4_runSel"]])!="" && isolate(input[["t_3_3_4_runSel"]])!=" ")
                {
                    f.id <- as.integer(isolate(input[["t_3_3_4_fileSel"]]))
                    fcs <- current.project$fcs.files[[f.id]]
                    if(is.defined(fcs))
                    {
                        f.name <- names(current.project$fcs.files)[f.id]
                        alg.id <- as.integer(isolate(input[["t_3_3_4_methodSel"]]))
                        run.id <- as.integer(isolate(input[["t_3_3_4_runSel"]]))

                        clust.col <- as.integer(current.project$test.files.clusters.col[[f.name]][[alg.id]][[run.id]])
                        fcs.clusters <- FPH.get.file.clusters(fcs, clust.col)
                        clust.sizes <- as.integer(unlist(sapply(fcs.clusters, function(cl)
                        {
                            return(cl[[1]])
                        })))

                        pur.mat <- computed.values$purity.matrix.clust[[f.name]][[alg.id]][[run.id]]

                        fcs.annot.file <- FPH.annotate.clusters.to.fcs(fcs,pur.mat,clust.col)
                        fcs.annot <- FPH.get.file.clusters(fcs.annot.file, ncol(fcs.annot.file@exprs))
                        annot.sizes <- as.integer(unlist(sapply(fcs.annot, function(annot)
                        {
                            return(annot[[1]])
                        })))

                        f.mat.annot <- computed.values$FG.matrices.annot[[f.name]][[alg.id]][[run.id]][[1]]
                        pop.names <- computed.values$pop.sizes[[f.name]][[as.integer(isolate(input[["t_3_3_4_methodSel"]]))]][[as.integer(isolate(input[["t_3_3_4_runSel"]]))]]
                        pop.names <- names(pop.names)
                        tmp <- FPH.map.test.to.ref(f.mat.annot)
                        annot.names <- 1:length(tmp)
                        annot.names <- paste0(pop.names[tmp],"__AN")

                        names(annot.sizes) <- annot.names


                        #CLUSTERS DETAILS-------------------------------------------------------------------------------------------------
                        if(is.defined(input[["t_3_3_4_purityByAnnot_slider"]]))
                        {
                            val <- compute.purity.points(fcs.clusters, fcs.annot, pur.mat)
                            points.id <- val[[1]]
                            purity.points <- val[[2]]
                            associated.annot <- val[[3]]

                            selected.clusters <- c()

                            output$t_3_3_4_purityByAnnot <- renderImage(
                            {
                                out.file <- tempfile(fileext = ".jpg")
                                jpeg(out.file, width = 640, height = 480)
                                pur.thresh <- as.numeric(input[["t_3_3_4_purityByAnnot_slider"]])
                                selected.clusters <<- plot.purity.by.annot(points.id, purity.points, annot.sizes, pur.thresh)
                                dev.off()

                                list(src=out.file)
                            })

                            output$t_3_3_4_clustersDetailsBelow <- renderTable(
                            {
                                "%not.in%" = Negate("%in%")
                                outer.clusters <- 1:length(associated.annot)
                                outer.clusters <- outer.clusters[unlist(which(outer.clusters%notin%selected.clusters))]

                                cl.table.below <- matrix(NA, nrow=length(outer.clusters), ncol=4)
                                colnames(cl.table.below) <- c("BELOW THRESHOLD","precision","Relative Size (annotation)", "Relative Size (file)")
                                if(length(outer.clusters)>0)
                                {
                                    sapply(1:length(outer.clusters), function(i)
                                    {
                                        cl.id <- outer.clusters[[i]]
                                        annot.id <- associated.annot[[cl.id]]
                                        cl.size.file <- trunc(clust.sizes[[cl.id]]/sum(as.integer(unlist(clust.sizes)))*100000)/1000
                                        cl.size.annot <- trunc(clust.sizes[[cl.id]]/as.integer(annot.sizes[[annot.id]])*100000)/1000
                                        cl.pur <- max(pur.mat[cl.id,])

                                        cl.table.below[i,1] <<- paste0("Cluster ", cl.id)
                                        cl.table.below[i,2] <<- trunc(cl.pur*10000)/10000
                                        cl.table.below[i,3] <<- paste0(cl.size.annot, " %")
                                        cl.table.below[i,4] <<- paste0(cl.size.file, " %")

                                    })
                                }
                                pur.thresh <- as.numeric(input[["t_3_3_4_purityByAnnot_slider"]])
                                return(cl.table.below)
                            })
                            output$t_3_3_4_clustersDetailsAbove <- renderTable(
                            {
                                "%not.in%" = Negate("%in%")
                                outer.clusters <- 1:length(associated.annot)
                                outer.clusters <- outer.clusters[unlist(which(outer.clusters%notin%selected.clusters))]

                                cl.table.above <- matrix(NA, nrow=length(selected.clusters), ncol=4)
                                colnames(cl.table.above) <- c("ABOVE THRESHOLD","precision","Relative Size (annotation)", "Relative Size (file)")
                                if(length(selected.clusters)>0)
                                {
                                    sapply(1:length(selected.clusters), function(i)
                                    {
                                        cl.id <- selected.clusters[[i]]
                                        annot.id <- associated.annot[[cl.id]]
                                        cl.size.file <- trunc(clust.sizes[[cl.id]]/sum(as.integer(unlist(clust.sizes)))*100000)/1000
                                        cl.size.annot <- trunc(clust.sizes[[cl.id]]/as.integer(annot.sizes[[annot.id]])*100000)/1000
                                        cl.pur <- max(pur.mat[cl.id,])

                                        cl.table.above[i,1] <<- paste0("Cluster ", cl.id)
                                        cl.table.above[i,2] <<- trunc(cl.pur*10000)/10000
                                        cl.table.above[i,3] <<- paste0(cl.size.annot, " %")
                                        cl.table.above[i,4] <<- paste0(cl.size.file, " %")
                                    })
                                }
                                pur.thresh <- as.numeric(input[["t_3_3_4_purityByAnnot_slider"]])
                                return(cl.table.above)
                            })
                        }

                        #POP DETAILS TABLES-------------------------------------------------------------------------------------------------
                        pop.col <- as.integer(current.project$ref.files.populations.col[[f.name]][[alg.id]][[run.id]])
                        fcs.pop <- FPH.get.file.clusters(fcs, pop.col)
                        pop.sizes <- as.integer(unlist(sapply(fcs.pop, function(curr.pop)
                        {
                            return(curr.pop[[1]])
                        })))

                        annot.pop.mapping.table <- matrix(NA, ncol=length(fcs.pop), nrow=length(fcs.annot))

                        pur.mat.annot <- computed.values$purity.matrix.annot[[f.name]][[alg.id]][[run.id]]
                        fp.mat <- f.mat.annot
                        lapply(1:nrow(fp.mat), function(ro)
                        {
                            lapply(1:ncol(fp.mat), function(co)
                            {
                                fs <- f.mat.annot[ro,co]
                                p.v <- pur.mat.annot[ro,co]
                                fp.mat[ro,co] <<- paste0("(",trunc(fs*1000)/1000,"-",trunc(p.v*1000)/1000,")")
                            })
                        })
                        annot.pop.mapping.table[,] <- fp.mat
                        colnames(annot.pop.mapping.table) <- pop.names
                        rownames(annot.pop.mapping.table) <- annot.names
                        
                        tmp.table <- as.data.frame(t(annot.pop.mapping.table))
                        tmp.table.colors <- matrix("white",ncol=ncol(fp.mat),nrow=nrow(fp.mat))
                        tmp.table.cuts <- as.vector(fp.mat)
                        
                        for(k in 1:nrow(f.mat.annot))
                        {
                            tmp.max.id <- which(f.mat.annot[k,]==max(f.mat.annot[k,]))[[1]]
                            tmp.table.colors[k,tmp.max.id] <- "gray"
                        }

                        output$t_3_3_4_MapAnnotPop <- DT::renderDataTable({
                            datatable(tmp.table)%>%formatStyle(names(tmp.table),backgroundColor=styleEqual(tmp.table.cuts, as.vector(tmp.table.colors)))
                        })



                        pop.details.table <- matrix(NA, nrow=length(fcs.pop), ncol=2)
                        mapping.pop.to.annot <- FPH.map.test.to.ref(t(f.mat.annot))

                        annot.names <- paste0(pop.names,"__AN")

                        lapply(1:length(fcs.pop), function(p.id)
                        {
                            pop.size.file <- trunc(pop.sizes[[p.id]]/sum(as.integer(unlist(pop.sizes)))*100000)/1000
                            pop.details.table[p.id, 1] <<- paste0(pop.size.file," %")
                            pop.details.table[p.id, 2] <<- paste0(mapping.pop.to.annot[[p.id]], " - ", annot.names[[p.id]])
                        })
                        rownames(pop.details.table) <- pop.names
                        colnames(pop.details.table) <- c("Relative size (file)", "Associated annotated group")
                        
                        output$t_3_3_4_populationsDetails <- DT::renderDataTable({
                            datatable(as.data.frame(pop.details.table))
                        })
                    }
                }
            }
        }
    })

    observe(#CLEAR UI
    {
        if(!clustering.algorithms$run.analysis)
        {
            updateSelectInput(session, "t_3_3_4_fileSel", "Select file", choices = list(), selected = NULL)
            updateSelectInput(session, "t_3_3_4_methodSel", "Select Algorithm", choices = list(), selected = NULL)
            updateSelectInput(session, "t_3_3_4_runSel", "Select Run", choices = list(), selected = NULL)
        }
    })

    observe(#EXPORT DATA
    {
        if(is.defined(current.project$name) && clustering.algorithms$run.analysis)
        {
            if(is.defined(current.project$fcs.files) && length(computed.values$FG.matrices.annot)>0)
            {
                if(is.defined(input[["t_3_3_4_fileSel"]]) && input[["t_3_3_4_fileSel"]]!="" && input[["t_3_3_4_fileSel"]]!=" " &&
                   is.defined(input[["t_3_3_4_methodSel"]]) && input[["t_3_3_4_methodSel"]]!="" && input[["t_3_3_4_methodSel"]]!=" " &&
                   is.defined(input[["t_3_3_4_runSel"]]) && input[["t_3_3_4_runSel"]]!="" && input[["t_3_3_4_runSel"]]!=" " &&
                   is.defined(input[["t_3_3_4_purityByAnnot_slider"]]))
                {
                    f.id <- as.integer(input[["t_3_3_4_fileSel"]])
                    fcs <- current.project$fcs.files[[f.id]]
                    if(is.defined(fcs))
                    {
                        f.name <- names(current.project$fcs.files)[f.id]
                        alg.id <- as.integer(input[["t_3_3_4_methodSel"]])
                        run.id <- as.integer(input[["t_3_3_4_runSel"]])

                        pur.mat <- computed.values$purity.matrix.clust[[f.name]][[alg.id]][[run.id]]
                        clust.col <- as.numeric(current.project$test.files.clusters.col[[f.name]][[alg.id]][[run.id]])
                        fcs.clusters <- FPH.get.file.clusters(fcs, clust.col)

                        added.col <- matrix(rep(F,nrow(fcs@exprs)),ncol=1)
                        colnames(added.col) <- paste0("purity.",ncol(fcs@exprs))

                        lapply(1:length(fcs.clusters), function(i)
                        {
                            cl.p <- max(pur.mat[i,])
                            if(cl.p >= as.numeric(input[["t_3_3_4_purityByAnnot_slider"]]))
                            {
                                added.col[as.integer(unlist(fcs.clusters[[i]][[2]]))] <<- T
                            }
                        })
                        fcs <- enrich.FCS(fcs, added.col)

                        purity.keyword <- as.numeric(input[["t_3_3_4_purityByAnnot_slider"]])
                        purity.keyword.name <- paste0("EXPPUR__",ncol(fcs@exprs),"__",clust.col)
                        fcs <- add.keyword.to.fcs(fcs, purity.keyword, purity.keyword.name)

                        fcs.out.name <- paste0(f.name, "_EXPPUR")

                        output$t_3_3_4_exportButton <- downloadHandler(
                            filename = function()
                            {
                                paste0(fcs.out.name,".fcs")
                            },
                            content = function(file)
                            {
                                write.enriched.FCS(fcs, file)
                            }
                        )
                    }
                }
            }
        }
    })
}

