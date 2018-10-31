#EVERY FCS FILE MUST HAVE A LAST COLUMN WITH THE CLUSTER ID OF EACH EVENT
FPH.load.folder <- function(path.to.fcs.folder)
{
    folder.content <- list.files(path.to.fcs.folder, pattern=".fcs")
    fcs.files <- lapply(folder.content, function(f.path)
    {
        return(read.FCS(paste0(path.to.fcs.folder,f.path), emptyValue = FALSE))
    })

    return(fcs.files)
}

FPH.get.file.name <- function(fcs.file)
{
    f.name <- strsplit(fcs.file@description$FILENAME, "/", fixed = T)
    f.name <- unlist(f.name[[1]][[length(f.name[[1]])]])
    f.name <- substr(f.name,1,nchar(f.name)-4)

    return(f.name)
}

FPH.get.file.clusters <- function(fcs.file, clust.col, maxclust=5000)
{
    f <- fcs.file
    f.mat <- f@exprs
    clust.dim <- min(ncol(f.mat),clust.col)
    clust.names <- unique(f@exprs[,clust.dim])
    clust.names <- clust.names[clust.names!=-1]
    nmb.clust <- length(clust.names)
    
    if(nmb.clust>maxclust)
    {
        nmb.clust <- maxclust
    }
    
    clusters.data <- lapply(1:nmb.clust, function(c)
    {
        events.id <- which(f.mat[,clust.dim] %in% clust.names[c])
        nmb.events <- length(events.id)
        
        out.list <- NULL
        if(nmb.events > 0)
        {
            out.list <- list(nmb.events, events.id)
        }
        return(out.list)
    })
    clusters.data[ sapply(clusters.data, is.null) ] <- NULL
    return(clusters.data)
}

FPH.map.test.to.ref <- function(scoring.matrix)
{
    output <- sapply(1:nrow(scoring.matrix), function(r)
    {
        return(which(scoring.matrix[r,]==max(scoring.matrix[r,]))[[1]])
    })
    return(output)
}

FPH.get.purity.matrix <- function(clusters.ref, clusters.test)
{
    purity.matrix <- matrix(0,nrow = length(clusters.test), ncol = length(clusters.ref))
    
    lapply(c(1:length(clusters.test)), function(l)
    {
        lapply(c(1:length(clusters.ref)), function(c)
        {
            purity.matrix[l,c] <<- length(which(clusters.test[[l]][[2]] %in% clusters.ref[[c]][[2]])) / clusters.test[[l]][[1]]
        })
    })
    
    return(purity.matrix)
}

FPH.get.prec.rec.matrices <- function(clusters.ref, clusters.test)
{
    prec.matrix <- matrix(0,nrow = length(clusters.test), ncol = length(clusters.ref))
    rec.matrix <- matrix(0,nrow = length(clusters.test), ncol = length(clusters.ref))
    
    lapply(c(1:length(clusters.test)), function(l)
    {
        lapply(c(1:length(clusters.ref)), function(c)
        {
            prec.matrix[l,c] <<- length(which(clusters.test[[l]][[2]] %in% clusters.ref[[c]][[2]])) / clusters.test[[l]][[1]]
            rec.matrix[l,c] <<- length(which(clusters.test[[l]][[2]] %in% clusters.ref[[c]][[2]])) / clusters.ref[[c]][[1]]
        })
    })
    
    return(list(prec.matrix,rec.matrix))
}

FPH.compute.F.G.matrix <- function(prec.rec.matrices)
{
    p <- prec.rec.matrices[[1]]
    r <- prec.rec.matrices[[2]]
    F.mat <- 2*p*r/(p+r)
    F.mat[unlist(which(is.na(F.mat)))] <- 0
    G.mat <- sqrt(p*r)
    
    return(list(F.mat,G.mat))
}

FPH.annotate.clusters.to.fcs <- function(fcs.ref, purity.matrix, cl.col, threshold=0)
{
    fcs.out <- fcs.ref
    clusters.test <- FPH.get.file.clusters(fcs.ref, cl.col)
    if(nrow(purity.matrix)>ncol(purity.matrix))
    {
        fcs.out.mat <- cbind(fcs.ref@exprs,0)
        nc <- ncol(fcs.out.mat)
        lapply(1:length(clusters.test), function(i)
        {
            clust.id <- which(purity.matrix[i,] == max(purity.matrix[i,]))[[1]]
            lapply(clusters.test[[i]][[2]], function(ev)
            {
                fcs.out.mat[ev,nc] <<- 0
                if(purity.matrix[i,clust.id] >= threshold)
                {
                    fcs.out.mat[ev,nc] <<- clust.id
                }
            })
        })
        colnames(fcs.out.mat)[nc] <- "Annotations"
        fcs.out <- flowFrame(fcs.out.mat)
    }
    else
    {
        fcs.out.mat <- cbind(fcs.ref@exprs,fcs.ref@exprs[,cl.col])
        nc <- ncol(fcs.out.mat)
        colnames(fcs.out.mat)[nc] <- "Annotations"
        fcs.out <- flowFrame(fcs.out.mat)
    }
    
    return(fcs.out)
}

FPH.retrieve.clusters.data.from.file <- function(fcs.file)
{
    analyses.list <- get.keywords.with.keypart.FCS(fcs.file,"CLMETH__")
    analyses.algo <- NULL
    analyses.markers <- NULL
    analyses.parameters <- NULL
    analyses.column <- NULL
    
    if(length(analyses.list)>0)
    {
        analyses.algo <- list()
        analyses.markers <- list()
        analyses.parameters <- list()
        analyses.column <- list()
        for(k in 1:length(analyses.list))
        {
            if(!is.null(analyses.list[[k]]) && analyses.list[[k]] != "NULL")
            {
                curr.analysis <- analyses.list[[k]]
                run.name <- strsplit(curr.analysis,"__", fixed = T)[[1]][2]
                
                if(is.null(analyses.algo[[run.name]]))
                {
                    analyses.algo[[run.name]] <- NULL
                    analyses.markers[[run.name]] <- NULL
                    analyses.parameters[[run.name]] <- NULL
                    analyses.column[[run.name]] <- NULL
                }
                analyses.algo[[run.name]] <- c(analyses.algo[[run.name]], curr.analysis)
                analyses.column[[run.name]] <- c(analyses.column[[run.name]], strsplit(curr.analysis,"__", fixed = T)[[1]][3])
                
                used.markers <- strsplit(curr.analysis,"__", fixed = T)[[1]][4]
                if(used.markers != "NULL")
                {
                    used.markers <- list(unlist(strsplit(used.markers, ".-.", fixed = T)[[1]]))
                }
                analyses.markers[[run.name]] <- c(analyses.markers[[run.name]], used.markers)
                
                used.parameters <- strsplit(curr.analysis,"__", fixed = T)[[1]][5]
                if(used.parameters != "NULL")
                {
                    used.parameters <- list(unlist(strsplit(used.parameters, ".-.", fixed = T)[[1]]))
                }
                analyses.parameters[[run.name]] <- c(analyses.parameters[[run.name]], used.parameters)
            }
        }
    }
    return(list(analyses.algo,analyses.markers,analyses.parameters, analyses.column))
}

FPH.retrieve.scores.data.from.file <- function(fcs.file)
{
    
    f.scores <- NULL
    g.scores <- NULL
    t.scores <- NULL
    col.list <- NULL
    
    analyses.list <- get.keywords.with.keypart.FCS(fcs.file,"CLFMETH__")
    if(length(analyses.list)>0)
    {
        f.scores <- list()
        col.list <- list()
        for(k in 1:length(analyses.list))
        {
            curr.analysis <- analyses.list[[k]]
            run.name <- strsplit(curr.analysis,"__", fixed = T)[[1]][2]
            
            if(is.null(f.scores[[run.name]]))
            {
                f.scores[[run.name]] <- NULL
                col.list[[run.name]] <- NULL
            }
            f.scores[[run.name]] <- c(f.scores[[run.name]], strsplit(curr.analysis,"__", fixed = T)[[1]][6])
            col.list[[run.name]] <- c(col.list[[run.name]], strsplit(curr.analysis,"__", fixed = T)[[1]][3])
        }
    }
    
    analyses.list <- get.keywords.with.keypart.FCS(fcs.file,"CLGMETH__")
    if(length(analyses.list)>0)
    {
        g.scores <- list()
        for(k in 1:length(analyses.list))
        {
            curr.analysis <- analyses.list[[k]]
            run.name <- strsplit(curr.analysis,"__", fixed = T)[[1]][2]
            
            if(is.null(g.scores[[run.name]]))
            {
                g.scores[[run.name]] <- NULL
            }
            g.scores[[run.name]] <- c(g.scores[[run.name]], strsplit(curr.analysis,"__", fixed = T)[[1]][6])
        }
    }
    
    analyses.list <- get.keywords.with.keypart.FCS(fcs.file,"CLTMETH__")
    if(length(analyses.list)>0)
    {
        g.scores <- list()
        for(k in 1:length(analyses.list))
        {
            curr.analysis <- analyses.list[[k]]
            run.name <- strsplit(curr.analysis,"__", fixed = T)[[1]][2]
            
            if(is.null(t.scores[[run.name]]))
            {
                t.scores[[run.name]] <- NULL
            }
            t.scores[[run.name]] <- c(t.scores[[run.name]], strsplit(curr.analysis,"__", fixed = T)[[1]][6])
        }
    }
    
    return(list(f.scores,g.scores,t.scores,col.list))
}

FPH.get.dim.names <- function(fcs.file)
{
    out.list <- lapply(1:ncol(fcs.file@exprs), function(j)
    {
        d <- fcs@description[[paste0("$P",j,"S")]]
        if(is.null(d) || !is.na(d) || d != "" || d != " ")
        {
            d <- colnames(fcs@exprs)[j]
        }
        names(d) <- NULL
        return(unlist(d))
    })
    names(out.list) <- unlist(out.list)
    
    return(out.list)
}

FPH.get.labels.from.mapping.file <- function(map.file, label.col)
{
    x <- map.file[,as.integer(label.col)]
    return(x)
}

FPH.link.labels.to.populations <- function(fcs.file, map.file, fcs.col, map.col)
{
    fcs.clusters <- FPH.get.file.clusters(fcs.file, fcs.col)
    map.labels <- FPH.get.labels.from.mapping.file(map.file, map.col)
    ordered.clusters.list <- rep(NA,length(fcs.clusters))
    
    for(i in 1:length(fcs.clusters))
    {
        id <-  as.integer( unique(fcs.file[fcs.clusters[[i]],fcs.col])[[1]] )
        ordered.clusters.list[[i]] <- map.labels[id]
    }
}