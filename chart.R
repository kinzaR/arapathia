
plot.heatmap <- function(path.vals,sample_type,colors=c("blue","gray","red"),sample.clust=T,variable.clust=F,labRow=NULL, labCol=NULL, sample_colors=NULL, scale=T){  
  if(sample.clust==F){
    colv <- NA
  } else {
    colv <- T
  }
  if(variable.clust==F){
    rowv <- NA
  } else {
    rowv <- T
  }
  if(is.null(labRow)){
    labRow <- rownames(path.vals)
  }
  if(is.null(labCol)){
    labCol <- colnames(path.vals)
  }
  if(is.null(sample_colors)){
    sample_colors <- rainbow(length(unique(sample_type)))
    names(sample_colors) <- unique(sample_type)
  }
  if(scale==T){
    path.vals <- t(apply(path.vals,1,function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))
  }
  heatmap(path.vals,margins = c(10,10),labRow = labRow, labCol = labCol, scale="none",Rowv=rowv,Colv=colv,ColSideColors = sample_colors[sample_type],col=colorRampPalette(colors)(256))  
}


plot.pca <- function(fit,cp1=1,cp2=2,colors=NULL,legend=NULL,legend_colors=NULL,cex=0.5,pch=20){
  if(is.null(colors)) colors <- "black"
  cpv1 <- fit$scores[,cp1]
  cpv2 <- fit$scores[,cp2]
  plot(cpv1,cpv2,xlab=paste("pc",cp1),ylab=paste("pc",cp2),col=colors,pch=pch,cex=cex)
  if(!is.null(legend)){
    legend("left",legend=legend,col=legend_colors,pch=pch,xpd=T,cex=cex,border=0)
  }
}


plot.multiple.pca <- function(fit,comps=1:3,colors,plot.variance=F,legend=NULL,legend_colors=NULL,cex=0.9,pch=20){
  combs <- combn(comps,2)
  ncombs <- ncol(combs)
  nn <- ncombs
  if(!is.null(legend)) nn <- nn + 1
  if(plot.variance==T) nn <- nn + 1  
  
  nr <- floor(sqrt(nn))
  nc <- ceiling((nn)/nr)
  par(mfrow=c(nr,nc))
  for(i in 1:(ncombs)){
    plot.pca(fit,cp1=combs[1,i],cp2=combs[2,i],colors=colors,cex=cex,pch=pch)    
  }
  if(!is.null(legend)){
    plot(1,type="n",axes=F,xlab="",ylab="")
    legend("center",legend=legend,col=legend_colors,pch=15,lwd=2,xpd=T,cex=0.8,border=NA,pt.cex=1.2)
  }
  if(plot.variance==T){
    plot.pca.variance(fit,acum=T,thresh=0.1)
  }
  par(mfrow=c(1,1))
}


