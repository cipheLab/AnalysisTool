m.no.transf.param = c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","I515-A","Time")
"%not.in%" = Negate("%in%")

m.transform.asinh <- function(fcs.file, channels, cofactor)
{
    fcs <- fcs.file
    fcs@exprs[,channels] <- asinh(fcs@exprs[,channels]/cofactor)
    
    return(fcs)
}

m.inv.transform.asinh <- function(fcs.file, channels, cofactor)
{
    fcs <- fcs.file
    fcs@exprs[,channels] <- sinh(fcs@exprs[,channels])*cofactor
    
    return(fcs)
}

m.transform.logicle <- function(flow.frame, markers = NULL, value = NULL) #logiclTransformCIPHE
{
    
    if(is.null(markers)){
        if(is.null(flow.frame@description[["SPILL"]])){
            markers.transform <- colnames(flow.frame)
        } else {
            markers.transform <- colnames(flow.frame@description[["SPILL"]])
        }
    } else {
        markers.transform <- markers
    }
    
    list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
    list.index <- gsub("N","", list.index)
    list.index <- gsub("\\$P","", list.index)
    
    if(is.null(value)){
        if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]])){
            r.values <- unlist(lapply(list.index, function(x)
                as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
            )
        } else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]])) {
            r.values <- unlist(lapply(list.index, function(x)
                as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
            )
        } else {
            r.values <- rep(90, length(list.index))
        }
    }
    else {
        r.values <- rep(value, length(list.index))
    }
    
    w.values <- (4.5-log10(262143/abs(r.values)))/2
    w.values[which(w.values<0)] <- 0.5
    w.values[which(is.infinite(w.values))] <- 0.5
    
    for(t in 1:length(markers.transform)){
        lgcl <- logicleTransform(w=w.values[t])
        flow.frame <- transform(flow.frame, transformList(markers.transform[t],lgcl))
    }
    
    return(flow.frame)
}

m.inv.transform.logicle <- function(flow.frame, markers = NULL, value = NULL) #inversLogiclTransformCIPHE
{
    if(is.null(markers)){
        if(is.null(flow.frame@description[["SPILL"]])){
            markers.transform <- colnames(flow.frame)
        } else {
            markers.transform <- colnames(flow.frame@description[["SPILL"]])
        }
    } else {
        markers.transform <- markers
    }
    
    list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
    list.index <- gsub("N","", list.index)
    list.index <- gsub("\\$P","", list.index)
    
    if(is.null(value)){
        if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]])) {
            r.values <- unlist(lapply(list.index, function(x) 
                as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
            ) 
        } else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]])) {
            r.values <- unlist(lapply(list.index, function(x) 
                as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
            )   
        } else {
            r.values <- rep(90, length(list.index))
        }
    }
    else {
        r.values <- rep(value, length(list.index))
    }
    
    w.values <- (4.5-log10(262144/abs(r.values)))/2
    w.values[which(w.values<0)] <- 0.5
    w.values[which(is.infinite(w.values))] <- 0.5
    
    flow.frame.inv <- flow.frame
    
    for(t in 1:length(markers.transform)){
        invLgcl <- inverseLogicleTransform(trans = logicleTransform(w=w.values[t]))
        flow.frame.inv <- transform(flow.frame.inv, transformList(markers.transform[t],invLgcl))
    }
    
    return(flow.frame.inv)
}

m.compensate <- function(fcs.file)
{
    fcs <- fcs.file
    if(!is.null(fcs@description[["SPILL"]]))
    {
        fcs <- compensate(fcs, fcs@description[["SPILL"]])
    }
    
    return(fcs)
}

m.inv.compensate <- function(x, spillover)  #deCompensateFlowFrame
{
    if(!is.null(spillover)){
        cols <- colnames(spillover)
        sel <- cols %in% colnames(x)
        if(!all(sel)) {
            stop(keyword(x)[["FILENAME"]], "\\nThe following parameters in the spillover matrix are not present in the flowFrame:\\n",
                 paste(cols[!sel], collapse=", "), call.=FALSE)
        }
        e <- exprs(x)
        e[, cols] <- e[, cols] %*% spillover
        exprs(x) = e
        return(x)
    } else {
        return(x)
    }
}