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
        if(length(points.list.xvalues)>=i && length(points.list.xvalues[[i]])>0 && length(points.list[[i]])>0)
        {
            for(j in 1:length(points.list.xvalues[[i]]))
            {
                x <- points.list.xvalues[[i]][[j]]
                x.id <- which(unlist(as.numeric(x.values))==as.numeric(x))[[1]]
                mat[mat.row.id+x.id, 1] <- as.numeric(points.list[[i]][[j]])
            }
        }
    }
    df <- data.frame(mat)
    
    plt.colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    plt.colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    plot(
        ggplot(df, aes(x=as.numeric(as.character(x.val)))) + 
            geom_area(aes(y=as.numeric(as.character(y.val)), fill=Population), position="stack") +
            scale_fill_manual("Annotations", values=alpha(plt.colors[1:nrow(mat)], 0.7)) +
            geom_point(aes(y=as.numeric(as.character(y.val))), position="stack", colour="red") +
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

plot.purity.by.pop <- function(pop.clusters, purity.val, pop.sizes, purity.threshold)
{
    pop.order <- order(pop.sizes)
    plt.colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    plt.colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    pur.pts <- c()
    for(p in 1:length(pop.clusters))
    {
        pop.id <- pop.order[p]
        if(length(pop.clusters[[pop.id]])>0)
        {
            for(cl in 1:length(pop.clusters[[pop.id]]))
            {
                pur.pts <- c(pur.pts, purity.val[[pop.id]][[cl]])
                names(pur.pts)[length(pur.pts)] <- as.integer(pop.clusters[[pop.id]][[cl]])
            }
        }
    }
    
    Y <- pur.pts
    X <- as.integer(names(pur.pts))
    mat <- matrix(ncol = 2, nrow=length(Y))
    mat[,1] <- 1:length(Y)
    mat[,2] <- Y
    colnames(mat) <- c("x","y")
    df.clusters <- data.frame(mat)
    
    
    pop.rect.start <- c()
    pop.rect.end <- c()
    for(p in 1:length(pop.clusters))
    {
        pop.id <- pop.order[p]
        print(pop.id)
        if(length(pop.clusters[[pop.id]])>0)
        {
            if(p==1)
            {
                pop.rect.start <- c(-Inf)
                pop.rect.end <- c(pop.rect.end, length(pop.clusters[[pop.id]])+0.5)
            }
            else
            {
                x1 <- pop.rect.end[[length(pop.rect.end)]]
                
                pop.rect.start <- c(pop.rect.start, x1)
                pop.rect.end <- c(pop.rect.end, x1+length(pop.clusters[[pop.id]]))
            }
        }
        else
        {
            x1 <- 0
            if(length(pop.rect.end)>0)
            {
                x1 <- pop.rect.end[[length(pop.rect.end)]]
            }
            pop.rect.start <- c(pop.rect.start, x1)
            pop.rect.end <- c(pop.rect.end, x1)
        }
    }
    pop.rect.end[[length(pop.rect.end)]] <- Inf
    df.pop <- data.frame(xstart=pop.rect.start, xend=pop.rect.end, pop=factor(names(pop.sizes)[pop.order], levels=names(pop.sizes)[pop.order]))
    
    plot(
        ggplot() +
            geom_rect(data=df.pop, aes(xmin=xstart, xmax=xend, ymin=-Inf, ymax=Inf, fill=pop), alpha=0.4) +
            scale_fill_manual("Populations", values = plt.colors[1:length(pop.sizes)]) +
            geom_line(data=df.clusters, aes(x,y)) +
            geom_point(data=df.clusters, aes(x,y)) +
            geom_hline(yintercept=purity.threshold) +
            xlab("Clusters") +
            ylab("Purity") + 
            ylim(c(0,1.05)) +
            xlim(c(0,length(X))) +
            scale_x_discrete(labels=NULL) + 
            annotate("text", min(as.numeric(X))+4, purity.threshold+0.03, label="Purity Threshold")
    )
    
    
    above.points <- which(Y>=purity.threshold)
    Y <- Y[above.points]
    X <- X[above.points]
    
    return(X)
}