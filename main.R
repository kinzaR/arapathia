

hipathia <- function(genes.vals, pathigraphs, maxnum = 100, nodeMethod = "pond", verbose = T, tol = 0.000001, divide=F, response.tol=0){
  
  results <- list()
  
  results$by.path <- lapply(pathigraphs,function(pathigraph){    
    
    if(verbose) print(paste(pathigraph$path.id, "-", pathigraph$path.name))
    
    res <- list()
    res$nodes.vals <- nodes.values.from.genes( genes.vals, pathigraph$graph )
    
    respaths <- all.path.values( res$nodes.vals, pathigraph$subgraphs, pathigraph$path.id, nodeMethod, maxnum, tol=tol, divide = divide, response.tol = response.tol )
    effector.respaths <- effector.path.values( res$nodes.vals, pathigraph$effector.subgraphs, pathigraph$path.id, nodeMethod, maxnum, tol=tol, divide = divide, response.tol = response.tol )
    
    res$path.vals <- respaths[[1]]
    res$convergence <- respaths[[2]]
    
    res$effector.path.vals <- effector.respaths[[1]]
    res$effector.convergence <- effector.respaths[[2]]
    
    return(res)
    
  })
  
  results$all$path.vals <- do.call("rbind", lapply(results$by.path, function(x){x$path.vals}))
  results$all$effector.path.vals <- do.call("rbind", lapply(results$by.path, function(x){x$effector.path.vals}))
  
  return(results)
}


nodes.values.from.genes <- function(genes.vals, ig, summ="per90"){
  genes.list <- V(ig)$genesList
  names(genes.list) <- V(ig)$name
  genes.list <- genes.list[!grepl("_func", names(genes.list))]
  nodes.vals <- matrix(NA, nrow=length(names(genes.list)), ncol=ncol(genes.vals), dimnames = list(names(genes.list), colnames(genes.vals)))
  for (node.name in names(genes.list)){
    genes <- genes.list[[node.name]]
    if( "/" %in% genes ){ #Then the node is a protein complex
      lists <- get.genes.lists( genes )
      probabilities.mat <- matrix( NA, nrow=0, ncol=ncol(genes.vals))
      for( list1 in lists ){
        if( length(list1) > 1 ){
          prob <- summarize.probabilities(genes.vals[list1,,drop=F], summ)
        }else{
          prob <- genes.vals[list1,,drop=F]
        }
        probabilities.mat <- rbind( probabilities.mat, prob )
      }
      nodes.vals[node.name,] <- apply(probabilities.mat, 2, min)
    }else{
      if (length(genes.list[[node.name]]) > 1){
        nodes.vals[node.name,] <- summarize.probabilities(genes.vals[genes.list[[node.name]],,drop=F], summ)
      }else if (length(genes.list[[node.name]]) == 1 && !is.na(genes.list[[node.name]])){
        nodes.vals[node.name,] <- data.matrix(genes.vals[genes.list[[node.name]],,drop=F])
      }else{
        nodes.vals[node.name,] <- rep(1, ncol(nodes.vals))
      }
    }
  }
  return(nodes.vals)
}


summarize.probabilities <- function(probabilities, summ="per90"){
  if (summ == "mean"){
    prob <- apply(probabilities, 2, mean, na.rm = TRUE)
  }else if(summ == "median"){
    prob <- apply(probabilities, 2, median, na.rm = TRUE)
  }else if (summ == "max"){
    prob <- apply(probabilities, 2, max, na.rm = TRUE)
  }else if (summ == "min"){
    prob <- apply(probabilities, 2, min, na.rm = TRUE)
  }else if (summ == "per90"){
    prob <- apply(probabilities, 2, quantile, 0.9, na.rm = TRUE)
  }else if (summ == "per95"){
    prob <- apply(probabilities, 2, quantile, 0.95, na.rm = TRUE)
  }else if (summ == "per99"){
    prob <- apply(probabilities, 2, quantile, 0.99, na.rm = TRUE)
  }else{
    print (paste ("The option for summarizing the probabilities", summ, "is not valid"))
  }
  return(prob)
}


get.genes.lists <- function( genes.list ){
  g.list <- NULL
  g.list.list <- list()
  while( length(genes.list) > 0 ){
    if( genes.list[[1]] != "/" ){
      g.list <- c( g.list, genes.list[[1]])
    }
    else{
      g.list.list[[length(g.list.list)+1]] <- g.list
      g.list <- NULL
    }
    genes.list <- genes.list[-1]
  }
  g.list.list[[length(g.list.list)+1]] <- g.list
  return(g.list.list)
}


all.path.values <- function( nodes.vals, subgraph.list, pathway.name, method="pond", maxnum=100, tol = 0.000001, divide = F, response.tol = 0 ){
  path.names <- paste0(pathway.name, "__", names(subgraph.list))
  path.vals <- matrix(0, ncol=ncol(nodes.vals), nrow = length(subgraph.list), dimnames = list(path.names, colnames(nodes.vals)) )
  signal.dif <- list()
  for( path in path.names){
    nodes <- unlist(strsplit(unlist(strsplit(path, "\\__"))[2], " - ", fixed=T))
    res <- path.value(nodes.vals, subgraph.list[[unlist(strsplit(path, "\\__"))[2]]], nodes[1], nodes[2], method, maxnum=maxnum, tol=tol, divide=divide, response.tol = response.tol)
    path.vals[path,] <- res[[1]]
    signal.dif[[path]] <- res[[2]]
  }
  return(list(path.vals, signal.dif))
}


effector.path.values <- function( nodes.vals, effector.subgraph.list, pathway.name, method="pond", maxnum=100, tol = 0.000001, divide = F, response.tol = 0 ){
  effector.path.names <- paste0(pathway.name, "__", names(effector.subgraph.list))
  effector.path.vals <- matrix(0, ncol=ncol(nodes.vals), nrow = length(effector.subgraph.list), dimnames = list(effector.path.names, colnames(nodes.vals)) )
  signal.dif <- list()
  for( path in effector.path.names){
    endnode <- unlist(strsplit(path, "\\__"))[2]
    ininodes <- V(effector.subgraph.list[[endnode]])$name[!V(effector.subgraph.list[[endnode]])$name%in%get.edgelist(effector.subgraph.list[[endnode]])[,2]]
    res <- path.value(nodes.vals, effector.subgraph.list[[endnode]], ininodes, endnode, method, maxnum=maxnum, tol=tol, divide=divide, response.tol = response.tol)
    effector.path.vals[path,] <- res[[1]]
    signal.dif[[path]] <- res[[2]]
  }
  return(list(effector.path.vals, signal.dif))
}


path.value <- function( nodes.vals, subgraph, ininodes, endnode, method="pond", maxnum = 100, tol = 0.000001, divide=F, response.tol = 0 ){
  
  # Initialize lists
  ready <- ininodes
  processed <- list()
  
  # Initialize node values
  node.signal <- matrix(NA, ncol=ncol(nodes.vals), nrow = length(V(subgraph)), dimnames = list(V(subgraph)$name, colnames(nodes.vals)))
  endnode.signal.dif <- 10
  
  num <- 0
  reached_last <- F
  while( length(ready) > 0 && num <= maxnum){
    num <- num + 1
    actnode <- ready[[1]]
    old.signal <- node.signal[actnode,]
    
    # Compute node signal
    if(divide && actnode != endnode){
      nfol <- length(incident(subgraph, actnode, mode="out"))
    }
    else{
      nfol <- 1
    }
    node.signal[actnode,] <- compute.node.signal2(actnode, nodes.vals[actnode,], node.signal, subgraph, method, response.tol) / nfol
    
    # Transmit signal
    nextnodes <- get.edgelist(subgraph)[incident(subgraph, actnode, mode="out"),2]
    dif <- old.signal - node.signal[actnode,]
    
    if(actnode==endnode){
      reached_last <- T
      if(!all(is.na(dif)))
        endnode.signal.dif <- c(endnode.signal.dif, sqrt(sum(dif^2)))
      #num <- num+1
    }
    if(all(is.na(old.signal)) || endnode.signal.dif[length(endnode.signal.dif)] > tol )
      ready <- unique(c(ready, nextnodes))
    ready <- ready[-1]
  }
  if(reached_last==F){
    endnode.signal.dif <- NA
  }
  return(list(node.signal[endnode,], endnode.signal.dif))
}


compute.node.signal2 <- function(actnode, node.val, node.signal, subgraph, method="pond", response.tol = 0){
  
  incis <- incident(subgraph, actnode, mode="in")
  
  if(length(incis)==0){
    signal <- rep(1, length(node.val))
    
  } else {
    
    # get activators and inhibitors signal
    prevs <- get.edgelist(subgraph)[incis,1]
    input_signals <- node.signal[prevs,,drop=F]
    nas <- is.na(input_signals[,1])
    prevs <- prevs[!nas]
    incis <- incis[!nas]
    input_signals <- input_signals[!nas,,drop=F]
    typeincis <- E(subgraph)$relation[incis]
    activators <- typeincis==1
    nactivators <- sum(activators)
    inhibitors <- typeincis==-1
    ninhibitors <- sum(inhibitors)
    activator_signals <- input_signals[activators,,drop=F]
    inhibitor_signals <- input_signals[inhibitors,,drop=F]
    
    if( method == "sum"){
      s1 <- prettyifelse(nactivators>0, colSums(activator_signals), rep(1,length(node.val)))
      s2 <- prettyifelse(ninhibitors>0, colSums(inhibitor_signals), rep(0,length(node.val)))
      signal <- s1-s2
    }
    else if( method == "pond"){
      s1 <- prettyifelse(nactivators>0, apply(1- activator_signals, 2, prod), rep(0,length(node.val)))
      s2 <- prettyifelse(ninhibitors>0, apply(1- inhibitor_signals, 2, prod), rep(1,length(node.val)))
      signal <- (1-s1)*s2
    }
    else if( method == "min"){
      s1 <- prettyifelse(nactivators>0, apply(activator_signals,2,min), rep(1,length(node.val)))
      s2 <- prettyifelse(ninhibitors>0, 1-apply(inhibitor_signals,2,max), rep(1,length(node.val)))
      signal <- s1*s2
    }
    else {
      stop("Unknown propagation rule")
    }
    
    # If signal too low, signal do not propagate
    if(sum(nas) == 0 && signal < response.tol)
      signal <- rep(0,length(node.val))
    
  }
  
  signal[signal>1] <- 1
  signal[signal<0] <- 0
  signal <- signal*node.val
  
  return(signal)
}


prettyifelse <- function(test,una,olaotra){
  if(test){
    return(una)
  } else {
    return(olaotra)
  }
}



