library(RColorBrewer)

draw.cumulated.filled.plots <- function(points.list, points.list.xvalues, max.height=2, x.values.range=c(0.15,1), x.lab="par", y.lab="val")
{
    x.values <- sort(unique(as.numeric(unlist(unlist(points.list.xvalues)))))
    
    mat <- matrix(0, ncol=3, nrow=length(points.list)*length(x.values))
    colnames(mat) <- c("y.val","Population","x.val")
    for(i in 1:length(points.list))
    {
        mat.row.id <- (i-1)*length(x.values)
        mat[(mat.row.id+1):(mat.row.id+length(x.values)), 2] <- names(points.list)[i]
        mat[(mat.row.id+1):(mat.row.id+length(x.values)), 3] <- x.values
        for(j in 1:length(points.list.xvalues[[i]]))
        {
            x <- points.list.xvalues[[i]][[j]]
            x.id <- which(unlist(as.numeric(x.values))==as.numeric(x))[[1]]
            mat[mat.row.id+x.id, 1] <- as.numeric(points.list[[i]][[j]])
        }
    }
    df <- data.frame(mat)
    
    plt.colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    plt.colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    plot(
        ggplot(df, aes(x=as.numeric(as.character(x.val)))) + 
        geom_area(aes(y=as.numeric(as.character(y.val)), fill=Population), position="stack")+
        scale_fill_manual("Annotations", values=plt.colors[1:nrow(mat)]) +
        geom_area(aes(y=as.numeric(as.character(y.val)), fill=Population), position="stack", colour="black", show_guide=F)+
        ylim(0,length(points.list)+1) +
        xlab(x.lab) +
        ylab(y.lab) +
        geom_hline(yintercept=length(points.list)) + 
        annotate("text", min(as.numeric(x.values)), length(points.list)+1, label="Max value") + 
        labs(title = paste0("Cumulated F-score for parameter ",x.lab)) +
        theme_bw() 
    )
    
}


draw.F.score.barplot <- function(F.score.matrix, populations.names, populations.sizes, plot.title = "")
{
    scores.list.col <- sapply(1:nrow(F.score.matrix), function(r)
    {
        return(which(F.score.matrix[r,]==max(F.score.matrix[r,]))[[1]])
    })
    pop.order <- order(populations.sizes)
    
    scores.list <- rep(0,ncol(F.score.matrix))
    for(i in 1:length(scores.list))
    {
        tmp.val <- which(scores.list.col==i)
        if( length(tmp.val)>0 )
        {
            scores.list[i] <- F.score.matrix[tmp.val[[1]],i]
        }
    }
    names(scores.list) <- as.character(unlist(populations.names))
    scores.list <- scores.list[pop.order]
    scores.colors <- grDevices::terrain.colors(length(scores.list))
    
    
    pop.sizes.list <- sort(populations.sizes)
    names(pop.sizes.list) <- paste0(names(scores.list),": size")
    
    bar.names <- sapply(1:(2*length(scores.list)), function(i)
    {
        if(i%%2==1)
        {
            return(names(pop.sizes.list)[[(i+1)/2]])
        }
        else
        {
            return(names(scores.list)[[i/2]])
        }
    })
    
    df.score <- data.frame(name=names(scores.list), value=unlist(scores.list))
    df.size <- data.frame(name=names(pop.sizes.list), value=unlist(pop.sizes.list))
    plot(
        
        ggplot() + 
        geom_bar(data=df.score, aes(x=factor(name,level=names(scores.list)),y=value, fill=1:length(scores.colors)), stat="identity") +
            scale_fill_gradientn("F-score", colors = scores.colors, labels=NULL) +
        geom_bar(data=df.size, aes(x=factor(name,level=names(pop.sizes.list)),y=value, color="red"), stat="identity", fill="red") +
            scale_color_discrete("Relative Size", labels=NULL) +
        scale_x_discrete(limits=factor(unlist(bar.names), levels=unlist(bar.names))) +
        ylim(0,1.05) +
        xlab("Populations") +
        ylab("") +
        labs(title = paste0("RUN - ",plot.title)) +
        coord_flip() +
        theme_bw() + 
        theme(legend.direction="horizontal")
    )
}


plot.selected.clusters <- function(val.mat, clusters, markers)
{
    highlighted = rep("gray7",nrow(val.mat))
    lapply(clusters, function(cl)
    {
        highlighted[unlist(as.integer(cl))] <<- "firebrick"
    })
    highlited.explicit.ids <- which(highlighted=="firebrick")
    plot(val.mat[,markers], col=highlighted, xlim=c(-0.5,4.5), ylim=c(-0.5,4.5), pch=".")
    points(val.mat[as.numeric(highlited.explicit.ids),markers], pch=".", col=highlighted[highlited.explicit.ids])
    lapply(1:length(clusters), function(cl.id)
    {
        cl <- clusters[[cl.id]]
        xco <- mean(val.mat[unlist(as.integer(cl)), markers[1]])
        yco <- mean(val.mat[unlist(as.integer(cl)), markers[2]])
        text(xco,yco, names(clusters)[cl.id], col = "darkgreen", cex=1.7)
    })
}

plot.purity.by.annot <- function(annot.clusters, purity.val, annot.sizes, purity.threshold)
{
    annot.order <- order(annot.sizes)
    
    pts <- c()
    lapply(1:length(annot.clusters), function(an)
    {
        an.id <- annot.order[an]
        lapply(1:length(annot.clusters[[an.id]]), function(cl)
        {
            pts <<- c(pts, as.numeric(purity.val[as.integer(annot.clusters[[an.id]][[cl]])]))
            names(pts)[length(pts)] <<- as.integer(annot.clusters[[an.id]][[cl]])
        })
    })
    
    Y <- pts
    X <- as.integer(names(pts))
    above.points <- which(Y>=purity.threshold)
    Y <- Y[above.points]
    X <- X[above.points]
    
    plt.colors <- topo.colors(length(annot.clusters))
    par(xpd = T, mar = par()$mar + c(0,0,0,7))
    plot(X, Y, ylim=c(0,1.05), xlim=c(0,max(X)+1))
    x0 <- 0
    y0 <- 0
    lapply(1:length(annot.clusters), function(an)
    {
        an.id <- annot.order[an]
        if(an==1)
        {
            x1 <- x0+length(annot.clusters[[an.id]])+0.5
            y1 <- 1.05
            rect(x0,y0,x1,y1, col = plt.colors[an])
            x0 <<- x1
        }
        else
        {
            x1 <- x0+length(annot.clusters[[an.id]])
            y1 <- 1.05
            rect(x0,y0,x1,y1, col = plt.colors[an])
            x0 <<- x1
        }
    })
    pts.order <- order(as.numeric(names(pts)))
    lines(as.integer(names(pts))[pts.order], pts[pts.order], type="b")
    abline(h=purity.threshold)
    legend(max(as.integer(names(pts)))+11,0.8, names(annot.sizes)[annot.order], pch=15, col=plt.colors)
    par(mar=c(5, 4, 4, 2) + 0.1)
    
    return(X)
}