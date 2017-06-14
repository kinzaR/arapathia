
read.expression.matrix <- function(file){
  data.matrix(read.table(file, header=T, sep="\t", stringsAsFactors=F, row.names=1, comment.char=""))  
}


read.experimental.design <- function(file){
  des <- read.table(file,header=F,stringsAsFactors=F,row.names=1)
  colnames(des)[1] <- c("group")
  return(des)
}


translate.ids <- function(ids,xref){
  
  # get translation
  ltids <- xref[ids]
  tids <- sapply(ltids,function(x){
    if(is.null(x)){
      return(NA)
    } else {
      return(x[1])
    }
  })
  names(tids) <- ids
  
  # check problems
  is_na <- is.na(tids)
  tids_count <- table(tids)
  #is_dup <- duplicated(tids)==T & is_na==F
  is_dup <- tids_count[as.character(tids)]>1 & is_na==F
  
  # separate ids
  translated_ids <- unique(ids[!is_na])
  untranslated_ids <- unique(ids[is_na])
  duplicated_ids <- unique(ids[is_dup])
  
  # compute counts and ratios
  translated_ids_count <- length(translated_ids)
  untranslated_ids_count <- length(untranslated_ids)
  duplicated_ids_count <- length(duplicated_ids)
  n <- length(ids)
  translated_ids_ratio <- translated_ids_count/n
  untranslated_ids_ratio <- untranslated_ids_count/n
  duplicated_ids_ratio <- duplicated_ids_count/n
  
  return(list(
    translation=tids,
    is_na=is_na,
    is_dup=is_dup,
    translated_ids=translated_ids,
    untranslated_ids=untranslated_ids,
    duplicated_ids=duplicated_ids,
    translated_ids_count=translated_ids_count,
    untranslated_ids_count=untranslated_ids_count,
    duplicated_ids_count=duplicated_ids_count,
    translated_ids_ratio=translated_ids_ratio,
    untranslated_ids_ratio=untranslated_ids_ratio,
    duplicated_ids_ratio=duplicated_ids_ratio
  )
  )
  
}


translate.matrix <- function(exp,xref,verbose=T){
  
  new_ids <- gsub("\\.[0123456789]$","",rownames(exp))
  
  # translate ids
  tt <- translate.ids(new_ids,xref)
  
  # filter untranslated ids
  exp2 <- exp[!tt$is_na,,drop=F]
  valid_translation <- tt$translation[!tt$is_na]
  
  # average duplicated ids
  raw_exp3 <- by(exp2,valid_translation,colMeans,na.rm=T)
  if(ncol(exp2)>1){
    exp3 <- do.call("rbind",raw_exp3)
  } else {
    exp3 <- matrix(raw_exp3,ncol=1)
    rownames(exp3) <- names(raw_exp3)
    colnames(exp3) <- colnames(exp2)
  }
  
  if(verbose==T){
    cat("translated ids = ",tt$translated_ids_count," (",format(digits=2,tt$translated_ids_ratio),") \n",sep="")
    cat("untranslated ids = ",tt$untranslated_ids_count," (",format(digits=2,tt$untranslated_ids_ratio),") \n",sep="")
    cat("multihit ids = ",sum(tt$duplicated_ids_count)," (",format(digits=2,tt$duplicated_ids_ratio),") \n",sep="")
  }
  
  attr(exp3,"translation") <- tt
  return(exp3)
  
}


get.pretty.path.name <- function(pathigraphs, pathnames,maxchar=NULL){
  
  prettynames <- unlist(lapply(pathnames, function(pathname){
    strname <- unlist(strsplit(pathname, "\\__"))
    pathway <- strname[1]
    labelid <- pathigraphs[[pathway]]$label.id
    ininode <- gsub("_", " ", strname[2])
    newini <- labelid[which(labelid[,1]==ininode),2]
    if(length(strname) == 3 ){
      endnode <- gsub("_", " ", (sub("_", "", strname[3])))
      newend <- labelid[which(labelid[,1]==endnode),2]
      name <- paste0(pathigraphs[[pathway]]$path.name, ": ", newini, " -> ", newend )
    }else{
      name <- paste0(pathigraphs[[pathway]]$path.name, ": ", newini )
    }
  }))
  
  if(!is.null(maxchar)) prettynames <- clip.names(prettynames,maxchar = maxchar)
  
  return(prettynames)
  
}

