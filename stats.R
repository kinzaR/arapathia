
normalize.data <- function(exp.data,by.quantiles=T,by.gene=T,percentil=F, truncation_percentil=NULL){
  
  # Normalize data matrix
  norm.data <- data.matrix(exp.data)
  if(!is.null(truncation_percentil)){
    if( truncation_percentil >= 0 & truncation_percentil <= 1){
      # Guarantees that truncation_percentil is in [0.5,1]
      if(truncation_percentil < (1-truncation_percentil)) truncation_percentil <- 1-truncation_percentil
      # Truncates by the percentil
      norm.data <- t(apply(norm.data, 1, function(x){
        x[which(x < quantile(x, 1-truncation_percentil, na.rm=T))] <- quantile(x, 1-truncation_percentil, na.rm=T)
        x[which(x > quantile(x, truncation_percentil, na.rm=T))] <- quantile(x, truncation_percentil, na.rm=T)
        x
      }))
    }
    else{
      stop("Parameter truncation_percentil must be in [0,1]")
    }
  }
  if(by.quantiles){
    norm.data <- normalize.quantiles(norm.data)
  }
  if(by.gene){
    if(percentil){
      norm.data <- t(apply(norm.data,1, function(x){ecdf(x)(x)}))
    }else{
      norm.data <- t(apply(norm.data, 1, function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))
    }
  } else {
    if(percentil){
      emp <- ecdf(norm.data)
      norm.data <- t(apply(norm.data,1,emp))
    }else{
      norm.data <- (norm.data-min(norm.data,na.rm=T))/(max(norm.data,na.rm=T)-min(norm.data,na.rm=T))
    }
  }
  
  colnames(norm.data) <- colnames(exp.data)
  rownames(norm.data) <- rownames(exp.data)
  
  return(norm.data)
}


add.missing.genes <- function(exp.data,genes,default=0.5){
  missing_genes <- setdiff(genes,rownames(exp.data))
  if(length(missing_genes>0)){
    if(ncol(exp.data)==1){
      fakemat <- default + as.matrix(mat.or.vec(nr=length(missing_genes),nc=ncol(exp.data)))
    } else {
      fakemat <- default + mat.or.vec(nr=length(missing_genes),nc=ncol(exp.data))
    }
    rownames(fakemat) <- missing_genes
    colnames(fakemat) <- colnames(exp.data)
    exp.data <- rbind(exp.data,fakemat)
    cat(missing_genes)
  }
  return(exp.data)
}

do.anova.test <- function(data,group,adjust=T){
  tests <- apply(data,1,anova.test.fun,group=group)
  out <- do.call("rbind",lapply(tests,function(x) c(summary(x$fit)[[1]][[1,"Pr(>F)"]])))
  out[is.na(out)] <- 1
  if(adjust==T){
    out <- cbind(out,p.adjust(out[,1],method="fdr"))
  } else {
    out <- cbind(out,out[,1])
  }
  out <- as.data.frame(out,stringsAsFactors=F)
  colnames(out) <- c("p.value","adj.p.value")
  return(out)
}

anova.test.fun <- function(values,group,verbose=F){
  group <- factor(group)
  fit <- aov(values ~ group)
  tuckey <- TukeyHSD(fit)
  if(verbose==T){
    summary(fit)
  }
  return(list(fit=fit,tuckey=tuckey))
}


do.cor <- function(sel.vals, design,adjust=T){
  
  data <- sel.vals[,design$sample]
  
  testData <- do.call("rbind",apply(data, 1, cor.test.fun, design$value))
  if(adjust==T){
    fdrData <- p.adjust(testData[,1], method = "fdr")
  } else {
    fdrData <- testData[,1]
  }
  data2 <- data.frame(testData[,1:3], fdrData,stringsAsFactors=F)
  
  colnames(data2) <- c("p.value", "UP/DOWN", "correlation", "FDRp.value")
  data2[data2$statistic>0,"UP/DOWN"] <- "UP"
  data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
  return(data2)
}


cor.test.fun <- function(x, values){
  r <- try(cor.test(x=as.numeric(x), y=as.numeric(values)))
  result <- cor.data.frame(r)
  return(result)
}


cor.data.frame <- function(wilcox){
  if (class(wilcox) == "try-error" | is.na( wilcox$p.value )){
    pvalue <- 1
    class <- "0"
    esti <- 0
  } else{
    pvalue <- wilcox$p.value
    esti <- wilcox$estimate[1]
    if (esti < 0){
      class <- "DOWN" ## regarding DISEASE
    } else if (esti > 0){
      class <- "UP" ## regarding DISEASE
    } else if (esti == 0){
      if (wilcox$conf.int[1] == 0){
        class <- "UP"
      } else if (wilcox$conf.int[2] == 0){
        class <- "DOWN"
      } else{
        class <- 0
      }
    }
  }
  result <- data.frame(pvalue, class, esti, stringsAsFactors=F)
  return(result)
}
#########################################################
#################### T Test #############################
#########################################################
do.t.test <- function(sel.vals, group.value, g1, g2, paired=T, verbose=T, adjust=T){
  g1_indexes <- which(group.value==g1)
  g2_indexes <- which(group.value==g2)
  stat.vals <- calculate.ttest.test(sel.vals,g1_indexes,g2_indexes,paired=paired,adjust=adjust)
  return(stat.vals)

}
calculate.ttest.test <- function( data, control, treatment, paired, adjust=T){
	testData <- do.call("rbind",apply(data, 1, t.test.test.fun, control, treatment, paired))
	if(adjust==T){
      fdrData <- p.adjust(testData[,1], method = "fdr")
    } else {
      fdrData <- testData[,1]
    }
	if(paired == TRUE){
    #dat <- apply(data, 1, ttest.test.fun, control, treatment, paired)
    #testData <- do.call("rbind",dat)
		data2 <- data.frame(testData, fdrData, stringsAsFactors=F)
	}else{
    #testData <- do.call("rbind",apply(data, 1, t.test.test.fun, control, treatment, paired))

    # Standardize statistic
    m <- length(control)*length(treatment)/2
    sigma <- sqrt(length(control)*length(treatment)*(length(control)+length(treatment)+1)/12)
    z <- (testData[,3]-m)/sigma
    data2 <- data.frame(testData[,1:2], z, fdrData,stringsAsFactors=F)
    data2[which(data2$pvalue==1 & data2$class == "0" & data2$fdrData == 1), "z"] <- 0
  }
  colnames(data2) <- c("p.value", "UP/DOWN", "statistic", "FDRp.value")
  data2[data2$statistic>0,"UP/DOWN"] <- "UP"
  data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
  return(data2)
}
ttest.test.fun <- function(x, control, treatment, paired){
  r <- try(t.test( as.numeric(x[treatment]) ~ as.numeric(x[control]), paired=paired))
  if (class(r) == "try-error"){
    pvalue <- 1
    class <- "0"
    stat <- 0
  } else{
    pvalue <- pvalue(r)
    stat <- statistic(r, "standardized")[1]
    if(stat > 0){
      class <- "UP"
    }else if (stat < 0){
      class <- "DOWN"
    }else{
      class <- "0"
    }
  }
  result <- data.frame(pvalue=pvalue, class=class, stat=stat, stringsAsFactors=F)
  return(result)
}

t.test.test.fun <- function(x, control, treatment, paired){
  
  r <- try(t.test(x=as.numeric(x[treatment]), y=as.numeric(x[control]) , conf.int=TRUE, alternative="two.sided", paired=paired))
  
  result <- t.test.data.frame(r)
  return(result)
}

t.test.data.frame <- function(ttest){
  if (class(ttest) == "try-error"){
    pvalue <- 1
    class <- "0"
    stat <- 0
  } else{
    pvalue <- ttest$p.value
    esti <- ttest$estimate
    stat <- ttest$statistic[[1]]
    if(is.na(stat)){
        pvalue <- 1
        stat <- 0
        class <- "0"
      }else{
      if (stat < 0){
        class <- "DOWN" ## regarding DISEASE
      } else if (stat > 0){
        class <- "UP" ## regarding DISEASE
      } else if (stat == 0){
        if(ttest$conf.int[1] == "NaN"){
          class <- 0
        } else {
          if (ttest$conf.int[1] == 0 ){

            class <- "UP"
          } else if (ttest$conf.int[2] == 0){
            class <- "DOWN"
          } else{
            class <- 0
          }
      }
      }
    }
  }

  result <- data.frame(pvalue, class, stat,stringsAsFactors=F)
  return(result)
}
#########################################################
#########################################################

#########################################################
#################### limma  #############################
#########################################################
do.limma.test <- function(sel.vals, group.value, g1, g2, paired=F, verbose=T, adjust=T){
  design <- cbind(group.value==g1, group.value==g2)+0

  g1_indexes <- which(group.value==g1)
  g2_indexes <- which(group.value==g2)

  stat.vals <- calculate.limma.test(sel.vals,design[,2],g1_indexes,g2_indexes,paired=paired,adjust=adjust)
  return(stat.vals)

}
calculate.limma.test <- function( data, group.value, control, treatment, paired=F, adjust=T){
    #testData <- do.call("rbind",apply(data, 1, limma.test.fun, group.value, control, treatment))
    testData <- limma.test.fun(data, group.value, control, treatment)
    if(adjust==T){
      fdrData <- p.adjust(testData[,1], method = "fdr")
    } else {
      fdrData <- testData[,1]
    }
    if(paired == TRUE){
    	data2 <- data.frame(testData, fdrData, stringsAsFactors=F)
  	}else{  
    # Standardize statistic
	    m <- length(control)*length(treatment)/2
	    sigma <- sqrt(length(control)*length(treatment)*(length(control)+length(treatment)+1)/12)
	    z <- (testData[,3]-m)/sigma
	    data2 <- data.frame(testData[,1:2], z, fdrData,stringsAsFactors=F)
	    data2[which(data2$pvalue==1 & data2$class == "0" & data2$fdrData == 1), "z"] <- 0
	}
  colnames(data2) <- c("p.value", "UP/DOWN", "statistic", "FDRp.value")
  data2[data2$statistic>0,"UP/DOWN"] <- "UP"
  data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
  return(data2)
}

limma.test.fun <- function(x, group.value, control, treatment){

  fit <- try(lmFit(x, group.value))
  fit2 <- eBayes(fit)
  result <- limma.data.frame(fit2)
  return(result)
}
limma.data.frame <- function(fit2){
  if (class(fit2) == "try-error"){
    pvalue <- 1
    class <- "0"
    stat <- 0
  } else{
    pvalue <- fit2$p.value
    esti <- fit2$coefficients
    stat <- fit2$t
    if(is.na(stat)){
        pvalue <- 1
        stat <- 0
        class <- "0"
      }else{
      if (stat < 0){
        class <- "DOWN" ## regarding DISEASE
      } else if (stat > 0){
        class <- "UP" ## regarding DISEASE
      } else if (stat == 0){
          class <- 0
      }
    }
  }

  result <- data.frame(pvalue, class, stat,stringsAsFactors=F)
  return(result)
}
#########################################################
#########################################################
do.wilcox <- function(sel.vals, group.value, g1, g2, paired=F, verbose=T, adjust=T){
  
  g1_indexes <- which(group.value==g1)
  g2_indexes <- which(group.value==g2)
  stat.vals <- calculate.wilcox.test(sel.vals,g2_indexes,g1_indexes,paired=paired,adjust=adjust)
  
  return(stat.vals)
  
}


#calculate.wilcox.test <- function( data, control, disease, paired, adjust=T){
#  if(paired == TRUE){
#    dat <- apply(data, 1, wilcoxsign.test.fun, control, disease)
#    testData <- do.call("rbind",dat)
#    if(adjust==T){
#      fdrData <- p.adjust(testData[,1], method = "fdr")
#    } else {
#      fdrData <- testData[,1]
#    }
#   data2 <- data.frame(testData, fdrData, stringsAsFactors=F)
#  }else{
#    testData <- do.call("rbind",apply(data, 1, wilcox.test.fun, control, disease, paired))
#    if(adjust==T){
#      fdrData <- p.adjust(testData[,1], method = "fdr")
#    } else {
#      fdrData <- testData[,1]
#    }
#    # Standardize statistic
#    m <- length(control)*length(disease)/2
#    sigma <- sqrt(length(control)*length(disease)*(length(control)+length(disease)+1)/12)
#    z <- (testData[,3]-m)/sigma
#    data2 <- data.frame(testData[,1:2], z, fdrData,stringsAsFactors=F)
#    data2[which(data2$pvalue==1 & data2$class == "0" & data2$fdrData == 1), "z"] <- 0
#  }
#  colnames(data2) <- c("p.value", "UP/DOWN", "statistic", "FDRp.value")
#  data2[data2$statistic>0,"UP/DOWN"] <- "UP"
#  data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
#  return(data2)
#}
calculate.wilcox.test <- function( data, control, disease, paired, adjust=T){
  testData <- do.call("rbind",apply(data, 1, wilcox.test.fun, control, disease, paired))
  if(adjust==T){
    fdrData <- p.adjust(testData[,1], method = "fdr")
  } else {
    fdrData <- testData[,1]
  }
  if(paired == TRUE){
    data2 <- data.frame(testData, fdrData, stringsAsFactors=F)
  }else{    
    # Standardize statistic
    m <- length(control)*length(disease)/2
    sigma <- sqrt(length(control)*length(disease)*(length(control)+length(disease)+1)/12)
    z <- (testData[,3]-m)/sigma
    data2 <- data.frame(testData[,1:2], z, fdrData,stringsAsFactors=F)
    data2[which(data2$pvalue==1 & data2$class == "0" & data2$fdrData == 1), "z"] <- 0
  }
  colnames(data2) <- c("p.value", "UP/DOWN", "statistic", "FDRp.value")
  data2[data2$statistic>0,"UP/DOWN"] <- "UP"
  data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
  return(data2)
}



wilcoxsign.test.fun <- function(x, control, disease){
  r <- try(wilcoxsign_test( as.numeric(x[disease]) ~ as.numeric(x[control])))
  if (class(r) == "try-error"){
    pvalue <- 1
    class <- "0"
    stat <- 0
  } else{
    pvalue <- pvalue(r)
    stat <- statistic(r, "standardized")[1]
    if(stat > 0){
      class <- "UP"
    }else if (stat < 0){
      class <- "DOWN"
    }else{
      class <- "0"
    }
  }
  result <- data.frame(pvalue=pvalue, class=class, stat=stat, stringsAsFactors=F)
  return(result)
}


wilcox.test.fun <- function(x, control, disease, paired){
  r <- try(wilcox.test(x=as.numeric(x[disease]), y=as.numeric(x[control]) , conf.int=TRUE, alternative="two.sided", paired=paired))
  result <- wilcox.data.frame(r)
  return(result)
}


wilcox.data.frame <- function(wilcox){
  if (class(wilcox) == "try-error"){
    pvalue <- 1
    class <- "0"
    stat <- 0
  } else{
    pvalue <- wilcox$p.value
    esti <- wilcox$estimate
    stat <- wilcox$statistic[[1]]
    if (esti < 0){
      class <- "DOWN" ## regarding DISEASE
    } else if (esti > 0){
      class <- "UP" ## regarding DISEASE
    } else if (esti == 0){
      print(wilcox$conf.int[1])
      if (wilcox$conf.int[1] == 0){
        class <- "UP"
      } else if (wilcox$conf.int[2] == 0){
        class <- "DOWN"
      } else{
        class <- 0
      }
    }
  }
  result <- data.frame(pvalue, class, stat,stringsAsFactors=F)
  return(result)
}

#########################################################
#########################################################

do.pca <- function(data,cor=F){
  fit <- princomp(t(data), cor=cor)
  fit$var <- fit$sdev^2
  fit$explain_var <- fit$var/sum(fit$var)
  fit$acum_explain_var <- cumsum(fit$explain_var)
  return(fit)
}


compute.difexp <- function(vals, control.string, disease.string, experimental.design){
  
  control <- experimental.design[which(experimental.design[,2]==control.string),1]
  disease <- experimental.design[which(experimental.design[,2]==disease.string),1]
  sgroup <- experimental.design[experimental.design[,2] == control.string | experimental.design[,2]== disease.string,2]
  names(sgroup) <- experimental.design[experimental.design[,2] == control.string | experimental.design[,2]== disease.string,1]
  design <- cbind(sgroup==control.string, sgroup==disease.string)+0
  colnames(design) <- c("grupo1", "grupo2")
  
  fit <- lmFit(vals[,c(control, disease)], design[c(control, disease),])
  cont.matrix <- makeContrasts(grupo1-grupo2, levels=design[c(control, disease),])
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  result <- data.frame(statistic=as.numeric(fit2$t), p.value=as.numeric(fit2$p.value), laterality=as.factor(fit2$t>0))
  rownames(result) <- rownames(fit2)
  
  return(result)
}


compute.node.difexp <- function(results,group1.label,group2.label,groups){
  
  difexp <- list()
  for(pathway in names(results$by.path)){
    print(pathway)
    difexp[[pathway]] <- compute.difexp2(results$by.path[[pathway]]$nodes.vals,group1.label,group2.label,groups)
  }
  return(difexp)
}

###################################################
###################################################
###################################################
#added from stats.R : i add all the file without reading 
###################################################
###################################################
###################################################


do.anova.test <- function(data,group,adjust=T){
  tests <- apply(data,1,anova.test.fun,group=group)
  out <- do.call("rbind",lapply(tests,function(x) c(summary(x$fit)[[1]][[1,"Pr(>F)"]])))
  out[is.na(out)] <- 1
  if(adjust==T){
    out <- cbind(out,p.adjust(out[,1],method="fdr"))
  } else {
    out <- cbind(out,out[,1])
  }
  out <- as.data.frame(out,stringsAsFactors=F)
  colnames(out) <- c("p.value","adj.p.value")
  return(out)
}
  
anova.test.fun <- function(values,group,verbose=F){
  group <- factor(group)
  fit <- aov(values ~ group)
  tuckey <- TukeyHSD(fit)
  if(verbose==T){
    summary(fit)
  }
  return(list(fit=fit,tuckey=tuckey))
}
  

do.cor <- function(sel.vals, design,adjust=T){
  
  data <- sel.vals[,design$sample]

  testData <- do.call("rbind",apply(data, 1, cor.test.fun, design$value))
  if(adjust==T){
    fdrData <- p.adjust(testData[,1], method = "fdr")
  } else {
    fdrData <- testData[,1]
  }
  data2 <- data.frame(testData[,1:3], fdrData,stringsAsFactors=F)

  colnames(data2) <- c("p.value", "UP/DOWN", "correlation", "FDRp.value")
  data2[data2$statistic>0,"UP/DOWN"] <- "UP"
  data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
  return(data2)
}

cor.test.fun <- function(x, values){
  r <- try(cor.test(x=as.numeric(x), y=as.numeric(values)))
  result <- cor.data.frame(r)
  return(result)
}

cor.data.frame <- function(wilcox){
  if (class(wilcox) == "try-error" | is.na( wilcox$p.value )){
    pvalue <- 1
    class <- "0"
    esti <- 0
  } else{
    pvalue <- wilcox$p.value
    esti <- wilcox$estimate[1]
    if (esti < 0){
      class <- "DOWN" ## regarding DISEASE
    } else if (esti > 0){
      class <- "UP" ## regarding DISEASE
    } else if (esti == 0){
      if (wilcox$conf.int[1] == 0){
        class <- "UP"
      } else if (wilcox$conf.int[2] == 0){
        class <- "DOWN"
      } else{
        class <- 0
      }
    }
  }
  result <- data.frame(pvalue, class, esti, stringsAsFactors=F)
  return(result)
}


do.wilcox <- function(sel.vals, group.value, g1, g2, paired=F, verbose=T, adjust=T){
  
  g1_indexes <- which(group.value==g1)
  g2_indexes <- which(group.value==g2)
  
  stat.vals <- calculate.wilcox.test(sel.vals,g2_indexes,g1_indexes,paired=paired,adjust=adjust)
  
  return(stat.vals)
  
}



calculate.wilcox.test <- function( data, control, disease, paired, adjust=T){
  testData <- do.call("rbind",apply(data, 1, wilcox.test.fun, control, disease, paired))
  if(adjust==T){
    fdrData <- p.adjust(testData[,1], method = "fdr")
  } else {
    fdrData <- testData[,1]
  }
  if(paired == TRUE){
    data2 <- data.frame(testData, fdrData, stringsAsFactors=F)
  }else{    
    # Standardize statistic
    m <- length(control)*length(disease)/2
    sigma <- sqrt(length(control)*length(disease)*(length(control)+length(disease)+1)/12)
    z <- (testData[,3]-m)/sigma
    data2 <- data.frame(testData[,1:2], z, fdrData,stringsAsFactors=F)
    data2[which(data2$pvalue==1 & data2$class == "0" & data2$fdrData == 1), "z"] <- 0
  }
  colnames(data2) <- c("p.value", "UP/DOWN", "statistic", "FDRp.value")
  data2[data2$statistic>0,"UP/DOWN"] <- "UP"
  data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
  return(data2)
}


#
# COMENTAR CON MARTA!!!!!!!!!!!!!
#
# calculate.wilcox.test <- function( data, control, disease, paired, adjust=T){
#   if(paired == TRUE){
#     dat <- apply(data, 1, wilcoxsign.test.fun, control, disease)
#     testData <- do.call("rbind",dat)
#     if(adjust==T){
#       fdrData <- p.adjust(testData[,1], method = "fdr")
#     } else {
#       fdrData <- testData[,1]
#     }
#     data2 <- data.frame(testData, fdrData, stringsAsFactors=F)
#   }else{
#     testData <- do.call("rbind",apply(data, 1, wilcox.test.fun, control, disease, paired))
#     if(adjust==T){
#       fdrData <- p.adjust(testData[,1], method = "fdr")
#     } else {
#       fdrData <- testData[,1]
#     }
#     # Standardize statistic
#     m <- length(control)*length(disease)/2
#     sigma <- sqrt(length(control)*length(disease)*(length(control)+length(disease)+1)/12)
#     z <- (testData[,3]-m)/sigma
#     data2 <- data.frame(testData[,1:2], z, fdrData,stringsAsFactors=F)
#     data2[which(data2$pvalue==1 & data2$class == "0" & data2$fdrData == 1), "z"] <- 0
#   }
#   colnames(data2) <- c("p.value", "UP/DOWN", "statistic", "FDRp.value")
#   data2[data2$statistic>0,"UP/DOWN"] <- "UP"
#   data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
#   return(data2)
# }


wilcoxsign.test.fun <- function(x, control, disease){
  r <- try(wilcoxsign_test( as.numeric(x[disease]) ~ as.numeric(x[control])))
  if (class(r) == "try-error"){
    pvalue <- 1
    class <- "0"
    stat <- 0
  } else{
    pvalue <- pvalue(r)
    stat <- statistic(r, "standardized")[1]
    if(stat > 0){
      class <- "UP"
    }else if (stat < 0){
      class <- "DOWN"
    }else{
      class <- "0"
    }
  }
  result <- data.frame(pvalue=pvalue, class=class, stat=stat, stringsAsFactors=F)
  return(result)
}


wilcox.test.fun <- function(x, control, disease, paired){
  r <- try(wilcox.test(x=as.numeric(x[disease]), y=as.numeric(x[control]) , conf.int=TRUE, alternative="two.sided", paired=paired))
  result <- wilcox.data.frame(r)
  return(result)
}

wilcox.data.frame <- function(wilcox){
  if (class(wilcox) == "try-error"){
    pvalue <- 1
    class <- "0"
    stat <- 0
  } else{
    pvalue <- wilcox$p.value
    esti <- wilcox$estimate
    stat <- wilcox$statistic[[1]]
    if (esti < 0){
      class <- "DOWN" ## regarding DISEASE
    } else if (esti > 0){
      class <- "UP" ## regarding DISEASE
    } else if (esti == 0){
      if (wilcox$conf.int[1] == 0){
        class <- "UP"
      } else if (wilcox$conf.int[2] == 0){
        class <- "DOWN"
      } else{
        class <- 0
      }
    }
  }
  result <- data.frame(pvalue, class, stat,stringsAsFactors=F)
  return(result)
}



classify <- function(vals, group, valid=NULL, crossk=5, verbose=T, repe=100){
  
  # define default valid groups
  if(is.null(valid)){
    valid <- unique(group)
  }
  
  # select samples
  sel <- which(group %in% valid)
  vals2 <- vals[,sel]
  group2 <- group[sel]
  
  # classify
  classifs <- lapply(1:repe, function(x){ svm(x=t(vals2), y = as.factor(group2), cross=crossk, type="nu-classification",nu=0.1)})
  allaccus <- sapply(classifs, "[[", "tot.accuracy")
            
  accu.mean <- mean(allaccus)
  accu.sd <- sd(allaccus)
  
  # verbose
  if(verbose==T){ 
    cat("Accuracy: mean:",round(digits=2, accu.mean),"% sd:", round(digits=2, accu.sd),sep="")
  }  
  
  return(list( svm = classifs, accuracies = allaccus, mean = accu.mean, sd= accu.sd))
  
}


compare.classifs.wilcox <- function(cl1, cl2, verbose=T, plot=T){
  w <- wilcox.test(cl1$accuracies, cl2$accuracies, conf.int=T, alternative="two.sided")
  if(verbose)
    print(wilcox.data.frame(w))
  if(plot)
    boxplot(cl1$accuracies, cl2$accuracies)
  return(w)
}



check.accuracy <- function(results, experimental.design, control.string, disease.string, output.folder="", save=F, crossk=5, n=10){
  
  control <- experimental.design[which(experimental.design[,2]==control.string),1]
  disease <- experimental.design[which(experimental.design[,2]==disease.string),1]
  raw.exp <- as.factor(experimental.design[,2])
  names(raw.exp) <- experimental.design[,1]  
  exp <- raw.exp[c(control, disease)]
  
  path.accu <- unlist(lapply(results$by.path, function(x){
    m <- t(x$path.vals[,c(control,disease)])
    # Removing paths with the same value for each sample
    good <- sapply(colnames(m), function(y){sum(duplicated(m[,y]))!=(nrow(m)-1)})
    if(sum(good)>1){
      svm <- svm(x=m[,good], y = exp, cross=crossk)
      return(svm$tot.accuracy)
    }else{
      return(NA)
    }
  }))
  results$all$accuracy$by.path <- path.accu
  
  m <- t(results$all$path.vals[,c(control, disease)])
  good <- sapply(colnames(m), function(x){sum(duplicated(m[,x]))!=(nrow(m)-1)})
  if(sum(good)>1){
    all.svm <- svm(x=m[,good], y = exp, cross=crossk)
    mat <- results$all$path.vals[good,]
    p <- predict(all.svm, t(mat))
    percent.all <- sum(p==raw.exp[names(p)])/ncol(mat)*100
    results$all$accuracy$all$total <- all.svm$tot.accuracy
    results$all$accuracy$all$percent <- percent.all
  }else{
    results$all$accuracy$all$total <- NA
    results$all$accuracy$all$percent <- NA
  }

  if(nrow(results$all$wilcox)>=n){
    ord <- order(results$all$wilcox[,"p.value"] )
    sigs <- rownames(results$all$wilcox[ord,])[1:n]
    is.sig <- colnames(m) %in% sigs
    sig.svm <- svm(x=m[,is.sig], y = exp, cross=crossk)
    mat.sig <- results$all$path.vals[is.sig,]
    p.sig <- predict(sig.svm, t(mat.sig))
    percent.sig <- sum(p.sig==raw.exp[names(p.sig)])/ncol(mat.sig)*100
    results$all$accuracy$nsigs$total <- sig.svm$tot.accuracy
    results$all$accuracy$nsigs$percent <- percent.sig
  }else{
    results$all$accuracy$nsigs$total <- NA
    results$all$accuracy$nsigs$percent <- NA
  }

  sigs <- rownames(results$all$wilcox[results$all$wilcox$p.value<0.05,])
  if(length(sigs)>0){
    is.sig <- colnames(m) %in% sigs
    sig.svm <- svm(x=m[,is.sig], y = exp, cross=crossk)
    mat.sig <- results$all$path.vals[is.sig,]
    p.sig <- predict(sig.svm, t(mat.sig))
    percent.sig <- sum(p.sig==raw.exp[names(p.sig)])/ncol(mat.sig)*100
    results$all$accuracy$sigs$total <- sig.svm$tot.accuracy
    results$all$accuracy$sigs$percent <- percent.sig
  }else{
    results$all$accuracy$sigs$total <- NA
    results$all$accuracy$sigs$percent <- NA
  }
  
  return(results)
  
}


empirical.pval <- function( results, genes.vals, pathigraphs, experimental.design, control.string, disease.string, deep=100, paired=F, maxnum=1000, nodeMethod="pond", verbose=F ){
  
  results.list <- lapply(1:deep, function(x){
    print(x)
    new.experimental.design <- experimental.design
    sam <- sample(nrow(experimental.design), nrow(experimental.design)/2)
    new.experimental.design[sam,2] <- control.string
    new.experimental.design[setdiff(1:nrow(experimental.design), sam),2] <- disease.string
    prettyways( genes.vals, pathigraphs, new.experimental.design, control.string, disease.string, paired, maxnum, nodeMethod, verbose )
  })
  
  results$all$empirical.pval <- do.call("cbind", lapply(1:deep, function(x){results.list[[x]]$all$wilcox[rownames(results.list[[1]]$all$wilcox),"statistic"]}))
  rownames(results$all$empirical.pval) <- rownames(results.list[[1]]$all$wilcox)
  
  return(results)
}



plot.heatmap <- function(path.vals,sample_type,colors=c("blue","gray","red"),sample.clust=T,variable.clust=F,labRow=NULL, labCol=NULL, sample_colors=NULL){  
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
  path.vals <- t(apply(path.vals,1,function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))
  #heatmap(path.vals,margins = par("mar")[1:2],labRow = labRow, labCol = labCol, scale="none",Rowv=rowv,Colv=colv,ColSideColors = sample_colors[sample_type],col=colorRampPalette(colors)(256))  
  heatmap(path.vals,margins = c(8,8),labRow = labRow, labCol = labCol, scale="none",Rowv=rowv,Colv=colv,ColSideColors = sample_colors[sample_type],col=colorRampPalette(colors)(256))  
}


# Compute differential expression 
compute.difexp <- function(vals, control.string, disease.string, experimental.design){
  
  control <- experimental.design[which(experimental.design[,2]==control.string),1]
  disease <- experimental.design[which(experimental.design[,2]==disease.string),1]
  sgroup <- experimental.design[experimental.design[,2] == control.string | experimental.design[,2]== disease.string,2]
  names(sgroup) <- experimental.design[experimental.design[,2] == control.string | experimental.design[,2]== disease.string,1]
  design <- cbind(sgroup==control.string, sgroup==disease.string)+0
  colnames(design) <- c("grupo1", "grupo2")
  
  fit <- lmFit(vals[,c(control, disease)], design[c(control, disease),])
  cont.matrix <- makeContrasts(grupo1-grupo2, levels=design[c(control, disease),])
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  result <- data.frame(statistic=as.numeric(fit2$t), p.value=as.numeric(fit2$p.value), laterality=as.factor(fit2$t>0))
  rownames(result) <- rownames(fit2)

  return(result)
  
}

# Compute differential expression 
# compute.difexp2 <- function(vals, group1.label ,group2.label, groups){
#     
#   vals <- vals[,which(groups==group1.label | groups==group2.label)]
#   groups <- groups[which(groups==group1.label | groups==group2.label)]
# 
#   design <- cbind(groups==group1.label, groups==group2.label)+0
#   colnames(design) <- c("group1", "group2")
#     
#   fit <- lmFit(vals,design)
#   cont.matrix <- makeContrasts(group1-group2, levels=design)
#   fit2 <- contrasts.fit(fit, cont.matrix)
#   fit2 <- eBayes(fit2)
#   result <- data.frame(statistic=as.numeric(fit2$t), p.value=as.numeric(fit2$p.value), laterality=as.factor(fit2$t>0))
#   rownames(result) <- rownames(fit2)
# 
#   return(result)
# }

compute.node.difexp <- function(results,group1.label,group2.label,groups){
  difexp <- list()
  for(pathway in names(results$by.path)){
    print(pathway)
    difexp[[pathway]] <- compute.difexp2(results$by.path[[pathway]]$nodes.vals,group1.label,group2.label,groups)
  }
  return(difexp)
}


# Compute recurrence of a path at the gene level
# difexp: diferential expression of genes
subgraph.recurrence.by.genes <- function(subgraph, difexp){
  genes <- unique(unlist(V(subgraph)$genesList))
  genes <- genes[genes!="/" & !is.na(genes)]
  sigs <- genes[difexp[genes,2]<=0.05] 
  rec <- sum(abs(difexp[sigs,1]))/sum(abs(difexp[genes,1]))
  return(rec)
}

# Compute recurrence of a path at the node level
# difexp: diferential expression of nodes
subgraph.recurrence.by.nodes <- function(subgraph, difexp, pathway){
  nodes <- paste0(pathway, "__", V(subgraph)$label)
  sigs <- nodes[difexp[nodes,2]<=0.05] 
  rec <- sum(abs(difexp[sigs,1]))/sum(abs(difexp[nodes,1]))
  return(rec)
}

# Media de los estadisticos de los genes
gene.stat.mean <- function(subgraph, difexp){
  genes <- unique(unlist(V(subgraph)$genesList))
  genes <- genes[genes!="/" & !is.na(genes)]
  rec <- mean(abs(difexp[genes,1]))  
  return(rec)
}

# Media de los estadisticos de los nodos
node.stat.mean <- function(subgraph, difexp, pathway){
  nodes <- paste0(pathway, "__", V(subgraph)$label)
  rec <- sum(abs(difexp[nodes,1]))/length(nodes)
  return(rec)
}
# Media de los estadisticos de los genes
gene.stat.quantile <- function(subgraph, difexp,q=.95){
  genes <- unique(unlist(V(subgraph)$genesList))
  genes <- genes[genes!="/" & !is.na(genes)]  
  rec <- quantile(abs(difexp[genes,1]),q)
  return(rec)
}

# Media de los estadisticos de los nodos
node.stat.quantile <- function(subgraph, difexp, pathway,q=.95){
  nodes <- paste0(pathway, "__", V(subgraph)$label)  
  rec <- quantile(abs(difexp[nodes,1]),q)
  return(rec)
}


path.recurrence <- function(path.vals, subpath.name, control, disease, use.true.props=F,prop.thresh=.75){
  con <- (density(path.vals[subpath.name, control],from=0,to=1,n=1001))
  dis <- (density(path.vals[subpath.name, disease],from=0,to=1,n=1001))
  n <- length(control) + length(disease)
  if(use.true.props){
    cy <- con$y * length(control)/n
    dy <- dis$y * length(disease)/n    
  } else {
    cy <- con$y/sum(con$y)
    dy <- dis$y/sum(dis$y)    
  }
  den <- dy / (dy+cy)
  den[is.na(den)] <- 0
  recurrence <- sum(approx(con$x,den,path.vals[subpath.name, disease])$y>prop.thresh)/length(disease)    
  return(recurrence)
}

path.sample.recurrence <- function(path.vals, subpath.name, control, disease, use.true.props=F,prop.thresh=.75){
  con <- (density(path.vals[subpath.name, control],from=0,to=1,n=1001))
  dis <- (density(path.vals[subpath.name, disease],from=0,to=1,n=1001))
  n <- length(control) + length(disease)
  if(use.true.props){
    cy <- con$y * length(control)/n
    dy <- dis$y * length(disease)/n    
  } else {
    cy <- con$y/sum(con$y)
    dy <- dis$y/sum(dis$y)    
  }
  den <- dy / (dy+cy)
  den[is.na(den)] <- 0
  case_recurrence <- approx(con$x,den,path.vals[subpath.name, disease])$y
  contol_recurrence <- approx(con$x,den,path.vals[subpath.name, control])$y
  recurrence <- c(case_recurrence,contol_recurrence)
  names(recurrence) <- c(disease,control)
  return(recurrence)
}

print.recurrence.plots <- function(results, pathigraphs, control.string, disease.string, experimental.design){
  
  difexp.gen <- compute.difexp(genes.vals, control.string, disease.string, experimental.design)
  nodes.vals <- do.call("rbind", lapply(results$by.path, function(x){x$nodes.vals}))
  rownames(nodes.vals) <- unlist(lapply(names(results$by.path), function(x){paste0(x, "__", get.pretty.name.node(pathigraphs[[x]], rownames(results$by.path[[x]]$nodes.vals)))}))
  difexp.nod <- compute.difexp(nodes.vals, control.string, disease.string, experimental.design)

  control <- experimental.design[which(experimental.design[,2]==control.string),1]
  disease <- experimental.design[which(experimental.design[,2]==disease.string),1]
    
  recgen <- unlist(lapply(pathigraphs, function(x){lapply(x$subgraphs, subgraph.recurrence.by.genes, difexp.gen)}))
  recnod <- unlist(lapply(pathigraphs, function(x){lapply(x$subgraphs, subgraph.recurrence.by.nodes, difexp.nod, x$path.id)}))
  medgen <- unlist(lapply(pathigraphs, function(x){lapply(x$subgraphs, gene.stat.mean, difexp.gen)}))
  mednod <- unlist(lapply(pathigraphs, function(x){lapply(x$subgraphs, node.stat.mean, difexp.nod, x$path.id)}))
  quangen <- unlist(lapply(pathigraphs, function(x){lapply(x$subgraphs, gene.stat.quantile, difexp.gen)}))
  quannod <- unlist(lapply(pathigraphs, function(x){lapply(x$subgraphs, node.stat.quantile, difexp.nod, x$path.id)}))
  lonpath <- unlist(lapply(pathigraphs, function(x){lapply(x$subgraphs, function(subgraph) length(V(subgraph)))}))

  allnames <- gsub("\\.", "__", names(recgen))
  allstat <- results$all$wilcox[allnames, "statistic"]
  allpval <- results$all$wilcox[allnames, "FDRp.value"]
  recpath <- unlist(lapply(allnames, path.recurrence, path.vals = results$all$path.vals, control=control, disease=disease))
  recpathsample <- do.call("rbind",lapply(allnames, path.sample.recurrence, path.vals = results$all$path.vals, control=control, disease=disease))
  rownames(recpathsample) <- allnames
  
#   names(allstat) <- allnames
#   names(allpval) <- allnames
#   names(recpath) <- allnames
#   names(recgen) <- allnames
#   names(recnod) <- allnames
#   names(medgen) <- allnames
#   names(mednod) <- allnames
#   names(quangen) <- allnames
#   names(quannod) <- allnames
#   names(lonpath) <- allnames
  
  results$all$recurrence <- as.data.frame(cbind(recpath, recgen, recnod, medgen, mednod,quangen,quannod))
  rownames(results$all$recurrence) <- allnames
  results$all$sample_recurrence <- recpathsample
  
  plot(recgen[!allpval<0.05], abs(allpval[!allpval<0.05]), main="Genes", xlab="Recurrencia gen: Sum(estad. sigs)/Sum(estad.)", ylab="pvalor del path", cex=0.3, col="#00000045", xlim=c(0,1), ylim=c(0,1))
  points(recgen[allpval<0.05], abs(allpval[allpval<0.05]), col="#ff000060", cex=0.3, xlim=c(0,1), ylim=c(0,1))
 
  plot(recgen, allstat, main="Genes", xlab="Recurrencia gen: Sum(estad. sigs)/Sum(estad.)", ylab="estadistico del path", cex=0.2, col="#00000045")
  points(recgen[allpval<0.05], allstat[allpval<0.05], col="#ff000060", cex=0.3, xlim=c(0,1), ylim=c(0,1))
  
  plot(recnod, abs(allpval), main="Nodes", xlab="Recurrencia nodo: Sum(estad. sigs)/Sum(estad.)", ylab="pvalor del path", cex=0.2, col="#00000045")
  points(recnod[allpval<0.05], abs(allpval[allpval<0.05]), col="#ff000060", cex=0.3, xlim=c(0,1), ylim=c(0,1))
  
  plot(recnod, allstat, main="Nodes", xlab="Recurrencia nodo: Sum(estad. sigs)/Sum(estad.)", ylab="estadístico del path", cex=0.2, col="#00000045")
  points(recnod[allpval<0.05], allstat[allpval<0.05], col="#ff000060", cex=0.3, xlim=c(0,1), ylim=c(0,1))
  
  plot(medgen, abs(allpval), main="Probes. Media", ylab="pvalor del path", xlab="Media de los estadísticos de los probes del path", cex=0.2, col="#00000045")
  points(medgen[allpval<0.05], abs((allpval[allpval<0.05])), col="#ff000060", cex=0.3, xlim=c(0,1), ylim=c(0,1))
  
  plot(recgen, recpath, cex=0.2, main="Probes", xlab="Recurrencia probe: Sum(estad. sigs)/Sum(estad.)", ylab="Recurrencia del path: Proporción de casos seguros", col="#00000045")
  points(recgen[allpval<0.05], recpath[allpval<0.05], col="#ff000060", cex=0.3, xlim=c(0,1), ylim=c(0,1))
  
  plot(recnod, recpath, cex=0.2, main="Nodos", xlab="Recurrencia nodo: Sum(estad. sigs)/Sum(estad.)",ylab="Recurrencia del path: Proporción de casos seguros", col="#00000045")
  points(recnod[allpval<0.05], recpath[allpval<0.05], col="#ff000060", cex=0.3, xlim=c(0,1), ylim=c(0,1))
  
  plot(medgen, recpath, cex=0.2, main="Probes. Media", xlab="Media de los estadísticos de los probes del path", ylab="Recurrencia del path: Proporción de casos seguros", col="#00000045")
  points(medgen[allpval<0.05], recpath[allpval<0.05], col="#ff000060", cex=0.3, xlim=c(0,1), ylim=c(0,1))
  
  plot(quangen, recpath, cex=0.2, main="Genes. Media", xlab="Q95 de los estadísticos de los genes del path", ylab="Recurrencia del path: Proporción de casos seguros", col="#00000045")
  points(quangen[allpval<0.05], recpath[allpval<0.05], col="#ff000060", cex=0.3, xlim=c(0,1), ylim=c(0,1))
  
  plot(quannod, recpath, cex=0.2, main="Nodos. Media", xlab="Q95 de los estadísticos de los nodos del path", ylab="Recurrencia del path: Proporción de casos seguros", col="#00000045")
  points(quannod[allpval<0.05], recpath[allpval<0.05], col="#ff000060", cex=0.3, xlim=c(0,1), ylim=c(0,1))
  
  return(results)
}


# PCA 

do.pca <- function(data,cor=F){
  fit <- princomp(t(data), cor=cor)
  fit$var <- fit$sdev^2
  fit$explain_var <- fit$var/sum(fit$var)
  fit$acum_explain_var <- cumsum(fit$explain_var)
  return(fit)
}

plot.pca.variance <- function(fit,thresh=0,acum=F){
  if(acum==F){
    barplot(fit$explain_var[ fit$explain_var>thresh],ylab="explain variance",xlab="",las=2,cex.names=0.5,ylim=c(0,1))
  } else {
    barplot(fit$acum_explain_var[ fit$acum_explain_var<(1-thresh)],ylab="acum explain variance",xlab="",las=2,cex.names=0.5,ylim=c(0,1))
  }
}

plot.pca <- function(fit,cp1=1,cp2=2,colors=NULL,legend=NULL,legend_colors=NULL,cex=0.5,pch=20,text=F){
  if(is.null(colors)) colors <- "black"
  cpv1 <- fit$scores[,cp1]
  cpv2 <- fit$scores[,cp2]
  if(text==F){
    plot(cpv1,cpv2,xlab=paste("pc",cp1),ylab=paste("pc",cp2),col=colors,pch=pch,cex=cex)
  } else {
    plot(cpv1,cpv2,xlab=paste("pc",cp1),ylab=paste("pc",cp2),type="n")
    text(cpv1,cpv2,labels = rownames(fit$scores),col=colors,pch=pch,cex=cex)
  }
  if(!is.null(legend)){
    legend("left",legend=legend,col=legend_colors,pch=pch,xpd=T,cex=cex,border=0)
  }
}

plot.multiple.pca <- function(fit,comps=1:3,colors,plot.variance=F,legend=NULL,legend_colors=NULL,cex=0.9,pch=20,text=F){
  combs <- combn(comps,2)
  ncombs <- ncol(combs)
  nn <- ncombs
  if(!is.null(legend)) nn <- nn + 1
  if(plot.variance==T) nn <- nn + 1  
  
  nr <- floor(sqrt(nn))
  nc <- ceiling((nn)/nr)
  par(mfrow=c(nr,nc))
  for(i in 1:(ncombs)){
    plot.pca(fit,cp1=combs[1,i],cp2=combs[2,i],colors=colors,cex=cex,pch=pch,text=text)
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



compute.zdiff <- function(r1,r2,n1,n2){
  if(is.na(r1) | is.na (r2)) return(NA)
  if(abs(r1-r2)<0.0000001){
    return(0)
  }else{
    return(((0.5 * log((1+r1)/(1-r1+.Machine$double.eps))) - (0.5 * log((1+r2)/(1-r2+.Machine$double.eps)))) / sqrt( (1/(n1-3)) + (1/(n2-3))))
  }
}


compute.correlation <- function(pathigraphs, results, c1.names){ 
  cor1 <- unlist(sapply(names(pathigraphs), function(g){
    l <- lapply(names(pathigraphs[[g]]$subgraphs_funs), function(namesub){sapply(V(pathigraphs[[g]]$subgraphs_funs[[namesub]])$name, function(node){
      if(node %in% rownames(results$by.path[[g]]$nodes.vals)){
        return(cor(results$by.path[[g]]$nodes.vals[node,c1.names], results$by.path[[g]]$path.vals[paste(g,namesub, sep="__"),c1.names]))
      }
      else{
        return(1)
      }
    })})
    names(l) <- names(pathigraphs[[g]]$subgraphs_funs)
    return(l)
  }))
  return(cor1)
}


# Calcula el pvalor de la diferencia entre la correlación de los genes con el valor de su pathway en cada uno de los dos casos
# c1.names y c2.names es el vector de nombres de los dos grupos que queremos comparar
compute.gene.correlation.difference <- function(pathigraphs, results, c1.names, c2.names){
  cor1 <- compute.correlation(pathigraphs, results, c1.names)
  cor2 <- compute.correlation(pathigraphs, results, c2.names)
  difs <- sapply(1:length(cor1), function(n){compute.zdiff(cor1[n], cor2[n], length(c1.names), length(c2.names))} )
  pvals <- sapply(1:length(difs), function(d){2*pnorm(-abs(difs[d]))})
  adj.pvals <- p.adjust(pvals, method="fdr")
  res <- cbind(difs, pvalue=pvals, adj.pvalue=adj.pvals, cor1, cor2)
  rownames(res) <- names(cor1)
  return(res)
}


# c1.names y c2.names es el vector de nombres de los dos grupos que queremos comparar

node.color.per.correlation.difference <- function(fpathigraphs, cordiff){
  
  dff <- do.call("rbind", lapply(rownames(cordiff), function(x){
    names <- unlist(strsplit(x, split="\\."))
    data.frame(pathway=names[1], path=names[2], node=names[3], dif = cordiff[x,1], pval=cordiff[x,"adj.pvalue"], stringsAsFactors=F)
  }))
  rownames(dff) <- rownames(cordiff)
  
  df <- lapply(unique(dff$pathway), function(pathway){
    nodes <- unique(dff[dff$pathway==pathway, "node"])
    out <- lapply(nodes, function(node){
      whichpaths <- which(dff$pathway==pathway & dff$node==node & dff$pval< 0.05)
      if(length(whichpaths) > 0 ){
        sel <- whichpaths[which(dff$pval[whichpaths]==min(dff$pval[whichpaths]))]
        if(all(dff[whichpaths, "dif"]>0) || all(dff[whichpaths, "dif"]<0)){
          updown <- ifelse(dff[sel, "dif"]>0, "up", "down")
        }else{
          updown <- "both"
        }    
      }else{
        sel <- which(dff$pathway==pathway & dff$node==node)[1]
        updown <- ifelse(dff[sel, "dif"]>0, "up", "down")
      }    
      return(data.frame(updown=updown, pval = dff[sel, "pval"], stringsAsFactors=F))    
    })
    names(out) <- nodes
    out <- do.call("rbind", out)  
    return(out)
  })
  names(df) <- unique(dff$pathway)
  
  cols <- lapply(names(df), function(path){
    newdf <- df[[path]][V(fpathigraphs[[path]]$graph)$name,]
    col <- get.colors.from.pval(newdf$updown, newdf$pval)
    names(col) <- V(fpathigraphs[[path]]$graph)$name
    return(col)
  })
  names(cols) <- names(df)
  return(cols)
}


node.color.per.differential.expression <- function(results, fpgs, control.string, disease.string, experimental.design){ 
  cols <- lapply(names(results$by.path), function(path){
    difexp.nod <- compute.difexp(results$by.path[[path]]$nodes.vals, control.string, disease.string, experimental.design)
    updown <- sapply(difexp.nod$statistic, function(x){if(x <0){"down"}else if(x > 0){"up"}else{"both"}})
    col <- get.colors.from.pval(updown, difexp.nod$p.value)
    names(col) <- rownames(results$by.path[[path]]$nodes.vals)
    # Add function colors
    toadd <- V(fpgs[[path]]$graph)$name[!V(fpgs[[path]]$graph)$name %in% rownames(difexp.nod)]
    coltoadd <- rep("white", length(toadd))
    names(coltoadd) <- toadd
    col <- c(col, coltoadd)
    return(col)
  })
  names(cols) <- names(results$by.path)
  return(cols)
}

node.color.per.node.expression <- function(results, fpgs, control.string, experimental.design){
  control <- experimental.design[which(experimental.design[,2]==control.string),1]
  cols <- lapply(names(results$by.path), function(path){
    exp.nod <- apply(results$by.path[[path]]$nodes.vals[,control], 1, mean)
    col <- grey(1-ecdf(exp.nod)(exp.nod))
    names(col) <- rownames(results$by.path[[path]]$nodes.vals)
    # Add function colors
    toadd <- V(fpgs[[path]]$graph)$name[!V(fpgs[[path]]$graph)$name %in% names(col)]
    coltoadd <- rep("white", length(toadd))
    names(coltoadd) <- toadd
    col <- c(col, coltoadd)
    return(col)
  })
  names(cols) <- names(results$by.path)
  return(cols)
}


# att.list is listed by pathways, and represent attributes for the nodes, ordered as in V()$name
summarize.atts <- function(att.list, att.names){
  df.list <- c()
  for(pathway in names(att.list[[1]])){
    df <- sapply(att.list, function(l){l[[pathway]]})
    colnames(df) <- att.names
    df.list[[pathway]] <- df
  }
  return(df.list)
}


