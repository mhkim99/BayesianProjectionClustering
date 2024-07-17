# # base R plot
# plot(Kcluster.obj$KL,
#      xlab='K: number of clusters',
#      ylab='KL',
#      main=sprintf('Case %d: fourier basis n from %d to %d\n select K = %d at threshold of %.1f', t, min(idxA),max(idxA), nCluster, thKL))
# abline(h=threshold, col='red')
#grid plot
kmeans_fitted.plot <- function(ls_sB,plotRow,subTitles,threshold,xText,yText){
  tb <- ls_sB$sBtb
  nClusters <- ls_sB$nClusters
  nPlot <- length(nClusters)
  plotCol <- nPlot/plotRow
  nCl <- nrow(tb)+1
  
  pushViewport(plotViewport(c(1,2,.5,1)))
  pushViewport(viewport(layout=grid.layout(nrow=plotRow,ncol=plotCol)))
  
  for(t in seq(nPlot)){
    pRow <- ceiling(t/plotCol)
    pCol <-  map_dbl(t%%plotCol, function(x) ifelse(x==0, plotCol, x))
    
    x_axis <- pretty(1:nCl)
    x_axis[1] <- 1
    x_rg <- range(min(x_axis)-1,max(x_axis)+1)
    y_axis <- pretty(c(0,1),n=4)
    y_rg <- range(-0.05,1.05)
    
    pushViewport(viewport(layout.pos.row=pRow,layout.pos.col=pCol))
    pushViewport(plotViewport(c(3,3,3,1),
                              xscale=range(x_rg), yscale=y_rg))
    
    grid.text(sprintf('(%s) %s',letters[t], subTitles[t]),
              x= unit(-1,'lines'), y=unit(1, 'npc')+unit(2,'lines'), just='left')
    
    grid.points(1:nCl, c(1,tb[,t]), default.units = 'native', size = unit(.5, "char"))
    grid.points(nClusters[t], tb[nClusters[t]-1,t], default.units = 'native', size = unit(.5, "char"), pch=19, gp=gpar(col='red'))
    grid.lines(x = x_rg,y = rep(threshold[t],2), default.units = 'native', gp= gpar(col='red',lty=2))
    grid.xaxis(x_axis,x_axis)
    grid.yaxis(y_axis,y_axis)
    grid.lines(x = c(min(x_rg), max(x_rg)), y = c(min(y_rg), min(y_rg)), default.units = "native")
    grid.lines(x = c(min(x_rg), min(x_rg)), y = c(min(y_rg), max(y_rg)), default.units = "native")
    if(pCol==1) grid.text(yText, rot=90, x= unit(-4, 'lines'))
    if(pRow==plotRow) grid.text(xText, y= unit(-3, 'lines'))
    
    popViewport()
    popViewport()
  }
  
  popViewport()
  popViewport()
}

