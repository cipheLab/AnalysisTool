add.keyword.to.fcs <- function(fcs, added.keyword, added.keyword.name)
{
    fcs.out <- fcs
    fcs.out@description[[added.keyword.name]] <- added.keyword
    
    return(fcs.out)
} 

is.defined <- function(obj)
{
    defined <- T
    a <- !is.na(obj)
    b <- !is.null(obj)
    
    if(length(obj)>1)
    {
        a <- T
        b <- T
    }
    
    if((length(a)+length(b))<2)
    {
        defined <- F
    }
    else
    {
        if(a==T && b==T)
        {
            defined <- T
        }
        else
        {
            defined <- F
        }
    }
    return(defined)
}

keyword.exists.FCS <- function(fcs, key.part)
{
    key.exists <- F
    i <- 1
    while(!key.exists && i <= length(fcs@description))
    {
        if(grepl(key.part, names(fcs@description)[i], fixed = T))
        {
           key.exists <- T
        }
        i <- i+1
    }
    
    return(key.exists)
}

get.keywords.with.keypart.FCS <- function(fcs, key.part)
{
    keywords.list <- c()
    i <- 1
    while(i <= length(fcs@description))
    {
        if(grepl(key.part, names(fcs@description)[i], fixed = T))
        {
            keywords.list <- c(keywords.list,fcs@description[[i]])
            names(keywords.list)[length(keywords.list)] <- names(fcs@description)[i]
        }
        i <- i+1
    }
    
    return(keywords.list)
}

modify.keyword.value <- function(old.keyword, separator, old.position, new.value)
{
    tmp <- strsplit(old.keyword, separator, T)[[1]]
    new.keyword <- ""
    if(length(tmp)>1)
    {
        lapply(1:length(tmp), function(curr.position)
        {
            if(curr.position != old.position)
            {
                new.keyword <<- paste0(new.keyword,tmp[[curr.position]])
            }
            else
            {
                new.keyword <<- paste0(new.keyword, new.value)
            }
            if(curr.position < length(tmp))
            {
                new.keyword <<- paste0(new.keyword, separator)
            }
        })
    }
    
    return(new.keyword)
}

enrich.FCS <- function(target.fcs, new.column, values.range = 262144, flow.data = T)
{
    if(!is.matrix(new.column))
    {
        new.column <- matrix(new.column,ncol=1)
    }
    if(is.null(colnames(new.column)))
    {
        colnames(new.column) <- paste0("new_",ncol(target.fcs@exprs))
    }
    new.mat <- cbind(target.fcs@exprs, new.column) 
    new.fcs <- flowFrame(new.mat, description = target.fcs@description)
    
    levels(new.fcs@parameters@data[[2]]) <- c(target.fcs@parameters@data[[2]], colnames(new.column))
    new.fcs@parameters@data[[2]] <- c(target.fcs@parameters@data[[2]], colnames(new.column))
    if(flow.data)
    {
        new.fcs@parameters@data[[3]] <- values.range
        new.fcs@parameters@data[[4]] <- (0*as.numeric(!flow.data) - 111*as.numeric(flow.data))
        new.fcs@parameters@data[[5]] <- values.range-1
    }
    else
    {
        new.fcs@parameters@data[[3]] <- c(target.fcs@parameters@data[[3]], max(as.numeric(new.column)))
        new.fcs@parameters@data[[4]] <- c(target.fcs@parameters@data[[4]], 0)
        new.fcs@parameters@data[[5]] <- c(target.fcs@parameters@data[[5]], max(as.numeric(new.column))-1)
    }
    
    
    if(flow.data)
    {
        if(!is.null(target.fcs@description[["APPLY COMPENSATION"]]))
        {
            new.fcs@description[["APPLY COMPENSATION"]] <- target.fcs@description[["APPLY COMPENSATION"]]
        }
        if(!is.null(target.fcs@description[["SPILL"]]))
        {
            new.fcs@description[["SPILL"]] <- target.fcs@description[["SPILL"]]
        }
    }
    if(!is.null(target.fcs@description[["$TIMESTEP"]]))
    {
        new.fcs@description[["$TIMESTEP"]] <- target.fcs@description[["$TIMESTEP"]]
    }
    
    return(new.fcs)
}

write.FCS.CIPHE <- function(fcs, fcs.path, values.range = 262144, flow.data = T, keywords.to.save = NULL)
{
    fcs.out <- flowFrame(fcs@exprs)
    
    descR <- description(fcs.out)
    if(!is.null(keywords.to.save))
    {
        lapply(keywords.to.save, function(key)
        {
            descR[[key]] <<- fcs@description[[key]]
        })
    }
    
    lapply(1:ncol(fcs@exprs),function(x)
    {
        if(flow.data)
        {
            descR[[paste0("$P",x,"R")]] <<- values.range
        }
        else
        {
            descR[[paste0("$P",x,"R")]] <<- description(fcs)[[paste0("$P",x,"R")]]
        }
        descR[[paste0("$P",x,"S")]] <<- description(fcs)[[paste0("$P",x,"S")]]
    })
    pd <- pData(parameters(fcs))
    for(p in seq_along(pd[["name"]]))
    {
        descR[[sprintf("flowCore_$P%sRmax", p)]] <- pd[p,"minRange"]
        descR[[sprintf("flowCore_$P%sRmin", p)]] <- pd[p,"maxRange"]
    }
    levels(fcs.out@parameters@data[[2]]) <- fcs@parameters@data[[2]]
    fcs.out@parameters@data[[2]] <- fcs@parameters@data[[2]]
    if(flow.data)
    {
        fcs.out@parameters@data[[3]] <- values.range
        fcs.out@parameters@data[[4]] <- 0
        fcs.out@parameters@data[[5]] <- values.range-1
    }
    else
    {
        fcs.out@parameters@data[[3]] <- c(fcs@parameters@data[[3]], max(as.numeric(new.column)))
        fcs.out@parameters@data[[4]] <- c(fcs@parameters@data[[4]], 0)
        fcs.outs@parameters@data[[5]] <- c(fcs@parameters@data[[5]], max(as.numeric(new.column))-1)
    }
    
    fcs.out <- flowFrame(fcs@exprs, description = descR, parameters = fcs.out@parameters)
    if(flow.data)
    {
        if(!is.null(fcs@description[["SPILL"]]))
        {
            fcs.out@description[["SPILL"]] <- fcs@description[["SPILL"]]
        }
        if(!is.null(fcs@description[["APPLY COMPENSATION"]]))
        {
            fcs.out@description[["APPLY COMPENSATION"]] <- fcs@description[["APPLY COMPENSATION"]]
        }
    }
    if(!is.null(fcs@description[["$TIMESTEP"]]))
    {
        fcs.out@description[["$TIMESTEP"]] <- fcs@description[["$TIMESTEP"]]
    }
    
    write.FCS(fcs.out, fcs.path, delimiter = '#')
}