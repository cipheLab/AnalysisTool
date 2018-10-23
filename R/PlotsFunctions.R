draw.cumulated.filled.plots <- function(points.list, points.list.xvalues, max.height=2, x.values.range=c(0.15,1), x.lab="par", y.lab="val")
{
    par(xpd = T, mar = par()$mar + c(0,0,0,8))
    plot(x.values.range,c(0,max.height),xlab = x.lab, ylab = y.lab)
    
    x.vals <- sort(unique(as.numeric(unlist(unlist(points.list.xvalues)))))
    y.vals <- rep(0,length(x.vals))
    
    plotting.colors <- topo.colors(length(points.list))
    
    
    lapply(1:length(points.list), function(i)
    {
        y.pts <- y.vals
        
        lapply(1:length(points.list.xvalues[[i]]), function(j)
        {
            x <- points.list.xvalues[[i]][[j]]
            x.id <- which(unlist(as.numeric(x.vals))==as.numeric(x))[[1]]
            y.pts[[x.id]] <<- y.pts[[x.id]] + unlist(as.numeric(points.list[[i]][[j]]))
            
        })
        
        polygon(c(x.vals,rev(x.vals)), c(y.pts,rev(y.vals)),col=plotting.colors[i])
        points(c(x.vals,rev(x.vals)),c(y.pts,rev(y.vals)),pch="o")
        
        y.vals <<- y.pts
    })
    
    
    legend(x.values.range[[2]]+2,max.height-2, names(points.list), pch=15, col=plotting.colors)
    par(mar=c(5, 4, 4, 2) + 0.1)
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
    plotting.colors <- grDevices::terrain.colors(length(scores.list))
    
    barplot(scores.list, main=paste0("F scores by population (ordered by %events) \n ", plot.title), horiz = T, xlab = "F-score",
            names.arg = names(scores.list), cex.names = 0.8, xlim = c(0,1.05), col = plotting.colors)
    
    
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