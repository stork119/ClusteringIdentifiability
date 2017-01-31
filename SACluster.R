#####################################
### Sensitivity Matrix Clustering
### GPL 2.0 license
### Authors:
### -- Michal Wlodarczyk
### -- Karol Nienaltowski - corresponding author
### -- Michal Komorowski
### email: karol.nienaltowski@gmail.com
#####################################

### libraries ###
try({package.list <- list("logging","spatstat", "fields")
package.load <- sapply(package.list, function(package.name){
  package.exist <- require(package.name, character.only = TRUE)
  if(!package.exist){
    install.packages(package.name)
    return(library(package.name, character.only = TRUE))
  }
  return(package.exist)
})
})


wd <- dirname(parent.frame(2)$ofile)
source(paste(wd,"/CanonicalCorrelation/CANCOR.R", sep = ""))
rm(wd)

#### SMC #### 
SMC <- function(S      = NULL,
                FIM    = NULL,
                labels = NULL,
                names  = NULL,
                zeta = 1,
                delta = 0.95,
                rules.delta = TRUE,
                K = 1,
                dir.output = "~/"
                ){
  # Description :
    # Class that implements functions for clustering and identifiability analysis of parameters
    # in the model for which we are able to compute Fisher Information Matrix and/or
    # Sensitivity Matrix.
  # Parameters  :
    # S          - optional - if S is null then FIM is required
    #            - sensitivity matrix; columns correspond to parameters; 
    #            - rows are realisations of sensitivities procedure
    #            - For ODE model with three parameters p1, p2 and p3 and
    #      two variables y1 and y2 defined as :
    #      dy1/dt = y1(t, y1, y2, p1, p2, p3)
    #      dy2/dt = y2(t, y1, y2, p1, p2, p3) 
    #      observed in times t1, t2, ..., tk
    #      S = [ dy1(t1)/dp1  dy1(t1)/dp2 dy1(t1)/dp3]
    #          [ dy2(t1)/dp1  dy2(t1)/dp2 dy2(t1)/dp3]
    #          [ dy1(t2)/dp1  dy1(t2)/dp2 dy1(t2)/dp3]
    #          [ dy2(t2)/dp1  dy2(t2)/dp2 dy2(t2)/dp3]
    #          ...
    #          [ dy1(tk)/dp1  dy1(tk)/dp2 dy1(tk)/dp3]
    #          [ dy2(tk)/dp1  dy2(tk)/dp2 dy2(tk)/dp3]
    # FIM         - optional - if FIM is null then S is required
    #             - Fisher Information Matrix = S^(T)*(S)
    # labels      - ids   of parameters 
    # names       - names of parameters 
    # zeta        - sensitivity threshold parameter according to article [1]
    # delta       - identifiability threshold parameter according to article [1]
    # K           - identifiability threshold parameter according to
    #               supplement to article [1] (not used in this version)
    # rules.delta - TRUE !
    # dir.output  - default folder of output
  # Return :
    # SMC class object with parameters
    # ind        - indices of sensitive parameters
    # S_all      - Sensitivity Matrix
    #            - null if in the input S was null
    # S          - Sensitivity Matrix reduced to sensitive parameters
    #            - null if in the input S was null
    #            - S = S_all[,ind]
    #            - null if in the input S was null
    # FIM_all    - Fisher Information Matrix
    # FIM        - Fisher Information Matrix reduced to sensitive parameters
    #            - FIM = FIM_all[ind, ind]
    # M          - S or FIM depending on class input
    #            - if as an input user set S then M = S,  M = FIM otherwise
    # labels     - parameters ids reduced to sensitive parameters
    # labels_all - parameters ids
    # names      - parameters names used in plots reduced to sensitive parameters
    # names_all  - parameters names used in plots 
    # zeta       - sensitivity threshold
    # delta      - identifiability threshold 
    # rules.delta - functions for identifiabilty analysis
    # K          - identifiability threshold 
    # rules.K    - functions for identifiabilty analysis
    # call       = sys.call(),
    # CANCOR     - object of class CANCOR  consisitng functions for calculating
    #              canonical correaltions 
  
  ### LOGGING
  removeHandler(names(getLogger()[['handlers']])) 
  dir.log.folder <- paste(dir.output, "logs/",sep = "")
  dir.create(dir.log.folder, recursive = TRUE, showWarnings = FALSE) 
  dir.log.file <- gsub(" ",
                       "-",
                       paste( dir.log.folder,
                              gsub( ":", 
                                    "-", 
                                    Sys.time(), 
                                    fixed = TRUE), 
                              ".log", 
                              sep = ""),
                       fixed = TRUE)

  basicConfig(level = "FINEST")
  addHandler(writeToFile, file = dir.log.file)
  removeHandler('basic.stdout') 
  
  
  ### MAIN CODE  
  if( is.null(S) && is.null(FIM) ){
    logerror("Enter value of FIM (Fisher Information Matrix) or S (Sensitivity Matrix)")
    stop("Enter value of FIM (Fisher Information Matrix) or S (Sensitivity Matrix)")
  }
  else if( is.null(FIM) ){
    names.default <- colnames(S)
    FIM <- t(S)%*%S
    CANCOR <- CANCOR.SM()
  }
  else {
    names.default <- colnames(FIM)
    CANCOR <- CANCOR.FIM()
  }
  
  if( is.null(labels) ){
    labels <- names.default
  }
  if( is.null(names) ){
    names <- names.default
  }
    
  S_all   <- S
  FIM_all <- FIM
  
  ind     <- (diag(FIM) > zeta)
  
  if( !is.null(S) ){
    S <- S[,ind]
  }
  
  FIM <- FIM[ind,ind]
  if(CANCOR$type == "SM"){
    M = S
  } else {
    M = FIM
  }
  
  val <- list(M   = M,
              S   = S,
              S_all = S_all,
              FIM = FIM,
              FIM_all = FIM_all,
              ind=ind,
              labels=labels[ind],
              labels_all = labels,
              names = names[ind],
              names_all = names,
              zeta = zeta, 
              delta = delta,
              rules.delta = rules.delta,
              K   = K,
              rules.K = !rules.delta,
              call = sys.call(),
              CANCOR = CANCOR)
  class(val) <- "SMC"
  return(val)
}
print.SMC <- function(x,...){
  loginfo("FIM")
  loginfo(paste(x$FIM, sep = " "))
  loginfo("FIM end")
}

info.SMC <- function(x,...){
  loginfo("FIM")
  loginfo(paste(x$FIM, sep = " "))
  loginfo("FIM end")
}

#### Clustering function ####
clusterident.SMC <- function(x,...) {
  # Definition :
  # Non-trival clustering function basing on canonical correlations between set of 
  # sensitivity vectors.
  # Parameters : none
    # n - parameters number
    # m - sensitivity matrix
  # Return     : clustering object (hclust)
    # merge  
    # height 
    # order  
    # labels 
    # representant 
    # width 
    # stability 
    # d_measure 
    # merge_type 
    # rep 
    # redundancy 
    # identifiability 
    # weakness
  # ATTENTION 
  # - if number of sensitive parameters is low code may throw errors, 
  #   it was not completly debugged
  # - there must be at least two sensitive parameters
  
  k <- 0 
  n <- length(x$labels)
  
  
  ## normalization of grammian 
  x$gramm_norm <- diag(x$FIM)#^2/sum(diag(x$FIM)^2)
  
  ## stability and weakness are parameters that were used to calucalte some 
  #  statistcs/heurystics during clustering process
  stability <- 0
  width     <- 0
  
  parameters.activate <- 1:n
  parameters.weakness <- rep(x = 0, times = n)
  
  cluster.size      <- c(rep(x=1, times=n), rep(x=0, times=n))
  cluster.sizerep   <- c(rep(x=1, times=n), rep(x=0, times=n))
  
  ## cluster.active - tells ids of active clusters i.e. clusters that in this step are leafs
  #                 - cluster.correlation[cluster.active, cluster.active] - matrix of correlations
  #                   between clusters in current clustering step
  cluster.active    <- c(rep(x=TRUE, times=n), rep(x=FALSE, times=n))
  
  cluster.representant    <- c(1:n, rep(x=0, times=n))
  cluster.name      <- c(-(1:n), rep(x=0, times=n))
  cluster.hist      <- array(0, c(n-1, 2))
  cluster.level     <- array(0, 2 * n)   # how many merges have been done in cluster
  cluster.merge_type  <- array(0, n - 1) # array to store additional info about merges
  cluster.d_measure   <- array(0, n-1)
  cluster.heights     <- array(0, n-1)
  
  cluster.redundancy  <- array(1, n)
  
  cluster.elements  <- array(0, c(2 * n, n)) # list of elements
  cluster.elements[,1] <- c(1:n, rep(x=0, times=n))
  cluster.elementsrep  <- array(0, c(2 * n, n)) # list of elements
  cluster.elementsrep[,1] <- c(1:n, rep(x=0, times=n))
  
  ## cluster.correlation - correlation between clusters
  ##                     - in the first step cluster == parameters
  ##                     - array of size 2 * n, 2 * n beacuse of number of clustering steps
  cluster.correlation <- array(-Inf, c(2 * n, 2 * n)) 
  cluster.correlation[1:n,1:n] <- sapply(1:n, function(i){
    sapply(1:n, function(j){
      ifelse(i < j, 
             {loginfo(paste(i,j));
             x$CANCOR$get_information_cancor(x$M,
                                    cluster.elementsrep[i, 1:cluster.sizerep[i]],
                                    cluster.elementsrep[j, 1:cluster.sizerep[j]])},
             -Inf)
    })
  })
  
  ## global.correlation - correlation with all parameters
  # global.correlation <- sapply(1:n, 
  #                              function(i){
  #                                x$CANCOR$get_correlations_cancor(x$M,
  #                                                        i,
  #                                                        ((1:n)[cluster.active])[-i])
  #                               })
  

  ## clustering procedure - n - 1 steps
  while(k < n - 1
       # && k < 1
  ) {
    
    correlationMatrix <- cluster.correlation[cluster.active, cluster.active]
    loginfo(paste(k, sum(cluster.active), dim(correlationMatrix)))
    
    ## indices of active clusters
    activate  <- (1:(2*n))[cluster.active]
    
    ## Choosing pair if clusters with the highest correlation between them (best)
    #  - if there are pairs with equal correlation value random pair is chosen
    #bestPair <- which(correlationMatrix == max(correlationMatrix), arr.ind = TRUE)[1,]
    bestPairList <- which(correlationMatrix == max(correlationMatrix), arr.ind = TRUE)
    bestPair <- bestPairList[sample(1:nrow(bestPairList), 1),]
    
    bestInformation  <- correlationMatrix[bestPair[1], bestPair[2]]
    
    ## numbers of clusters with the highest correlation
    i <- activate[bestPair[1]]
    j <- activate[bestPair[2]]
    
    ## best clusters gonna be merged -- they are not leafs (not active) anymore
    cluster.active[i] <- FALSE
    cluster.active[j] <- FALSE
    
    ## statistics of best clusters
    bestCorrelation  <- x$CANCOR$get_correlations_cancor(
      x$M,
      cluster.elementsrep[i, 1:cluster.sizerep[i]],
      cluster.elementsrep[j, 1:cluster.sizerep[j]])
    
    k <- k + 1   
    
    cluster.hist[k,]      <- c(cluster.name[i], cluster.name[j])
    cluster.heights[k]    <- max(cluster.level[i], cluster.level[j]) + mean(1 - bestCorrelation*bestCorrelation)
    cluster.merge_type[k] <- 1#sum(best_cor[1]> delta) #1 #CHANGED
    
    
    ## updating tree's statistics
    width <- width + cluster.size[i] * cluster.size[j]
    #stability <- stability + 1/( max(cluster.heights[i] + 1, cluster.heights[j]+1)) * bestInformation
    
    cluster.size[n + k]    <- cluster.size[i] + cluster.size[j]
    cluster.sizerep[n + k] <- cluster.sizerep[i] + cluster.sizerep[j]
    
    cluster.level[n + k]   <- max(cluster.level[i], cluster.level[j]) + mean(1-bestCorrelation*bestCorrelation)
    
    cluster.name[n + k]      <- k;
    cluster.elements[n + k,] <- c(cluster.elements[i, 1:cluster.size[i]],
                                  cluster.elements[j, 1:cluster.size[j]],
                                  array(0, n - cluster.size[i] - cluster.size[j]))
    cluster.elementsrep[n + k,] <- c(cluster.elementsrep[i, 1:cluster.sizerep[i]],
                                     cluster.elementsrep[j, 1:cluster.sizerep[j]],
                                     array(0, n - cluster.sizerep[i]- cluster.sizerep[j]))
    
    ZL <- c(cluster.elementsrep[i, 1:cluster.sizerep[i]])
    ZR <- c(cluster.elementsrep[j, 1:cluster.sizerep[j]])
    Z  <- c(ZL,ZR)
    
    Zscore.K <- sapply(1:length(Z), function(z){
      (1 - x$CANCOR$get_correlations_cancor(x$M, Z[z], Z[-z]))*x$gramm_norm[Z[z]]
    })
    
    Zscore.delta <- sapply(1:length(Z), function(z){
      x$CANCOR$get_correlations_cancor(x$M, Z[z], Z[-z])
    })
    
    if(x$rules.delta){
      ile <- ifelse(max(Zscore.delta) > x$delta, 1, 0 )
    } else {
      ile <- ifelse(max(Zscore.K) < x$K, 1, 0 )
    }
    
    ## Removing non identifiable parameters from merged cluster !
    #  - maximizing subset of parameters with identifible parameters within this subset
    #  - non-identifible parameters are removed from cluster 
    #  - Removing parameters procedure consist of a couple steps
    #    (1) Zscore reduction :
    #        Choose parameters with the biggest Zscore value according to one of two conditions:
    #        - delta condition
    #        - K condition
    #    (2) Local correlation condition :
    #        Calculate correlation between parameters in cluster (locally).
    #        Choose parameter that have the largest local correlation
    #    (3) Local out of thrown correlation condition : 
    #        Calculate correlation between parameters in cluster (locally), 
    #        that in previous steps were chosen to remove
    #        Choose parameter that have the largest local correlation
    #    (4) Sensitivity condition :
    #        Choose parameter with the lowest sensitivity
    #    (5) Sampling condition :
    #        Sample parameter to remove
    
    if(ile){
  
      while(ile && length(ZL) > 0 && length(ZR) > 0){
        
        ## Removing parameters procedure - (1)
        if(x$rules.delta){
          remove <- which(Zscore.delta == max(Zscore.delta)) ### local
          loginfo(paste(("REDUCE: Zscore reduction Zscore.delta"), Zscore.delta, remove))
        } else {
          remove <- which(Zscore.K == max(Zscore.K)) ### local
          loginfo(paste(("REDUCE: Zscore reduction Zscore.K"), Zscore.K, remove))
        }
        
        ## Removing parameters procedure - (2)
        local.correlation <- sapply(remove, function(r){
          x$CANCOR$get_correlations_cancor(x$M, 
                                  Z[r],
                                  parameters.activate[!parameters.activate %in% Z[r]])
        })
        remove <- remove[which(local.correlation == max(local.correlation))]
        
        loginfo(paste(("REDUCE: Local correlation condition"), local.correlation, remove))
        
        ## Removing parameters procedure - (3)
        if(sum(!parameters.activate %in% Z[remove]) > 0){
          local.correalation.out <-  sapply(remove, function(r){
            x$CANCOR$get_correlations_cancor(x$M, 
                                    Z[r],
                                    parameters.activate[!parameters.activate %in% Z[remove]])
          })
        
          remove <- remove[which(local.correalation.out == max(local.correalation.out))]
          loginfo(paste(("REDUCE: Local out of thrown correlation condition"), local.correalation.out, remove))
      
        }
        
        ## Removing parameters procedure - (4)
        local.sensitivty <- diag(x$FIM)[Z[remove]]
        remove <- remove[which(local.sensitivty == min(local.sensitivty))]
        loginfo(paste(("REDUCE: Sensitivity condition"), local.sensitivty, remove))
        
        ## We check number of parameters chosen for removing using conditions (1-4)
        if(length(remove) > 1){
          parameters.weakness[Z[remove]] <- parameters.weakness[Z[remove]] + 1
        }
        
        ## Removing parameters procedure - (5)
        remove <- remove[sample(1:length(remove), 1)]
        loginfo(paste(("REDUCE: Sampling condition"),  remove))
        
        ## Removing parameters procedure & updating tree
        parameters.activate <- parameters.activate[!parameters.activate %in% Z[remove]]
        Z      <- Z[-remove]
        
        if(remove > length(ZL)){
          ZR <- ZR[-(remove - length(ZL))]    				
        } else {				
          ZL <- ZL[-remove]				
        }
        
        if(length(ZL) == 0 || length(ZR) == 0){		
          ile <- 0
        } else {
          
          Zscore.delta <- sapply(1:length(Z), function(z){
            x$CANCOR$get_correlations_cancor(x$M, Z[z], Z[-z])
          })
          Zscore.K <- sapply(1:length(Z), function(z){
            (1 - x$CANCOR$get_correlations_cancor(x$M, Z[z], Z[-z]))*x$gramm_norm[Z[z]]
          })
          cluster.redundancy[Z] <- Zscore.delta
          
          if(x$rules.delta){
            ile <- ifelse(max(Zscore.delta) > x$delta, 1, 0 )
          } else {
            ile <- ifelse(max(Zscore.K) < x$K, 1, 0 )
          }
        }
      }
      
      cluster.elementsrep[n + k,] <- c(Z, array(0, n - length(Z)))
      cluster.sizerep[n + k] <- length(Z)		
    }
    
    cluster.correlation[cluster.active, n+k] <- 
      sapply((1:(2*n))[cluster.active], function(i){
        x$CANCOR$get_information_cancor(x$M,
                               cluster.elementsrep[i, 1:cluster.sizerep[i]],
                               cluster.elementsrep[n+k, 1:cluster.sizerep[n+k]])
      }) 
    
    cluster.active[n+k] <- TRUE
    
  }
  
  rep <- cluster.elementsrep[2*n-1,1:cluster.sizerep[2*n-1]]
  identifiability <- rep(0, length(x$names_all))
  identifiability[x$ind][rep] <- rep(1, length(rep))
  
  ## hclust object
  cluster <- list(merge  = cluster.hist,
                  height = cluster.heights,
                  order  = cluster.elements[2 * n - 1,], # the order is constructed not to make lines cross in dendrogram
                  labels = x$labels,
                  representant = cluster.representant,
                  width = width,
                  stability  = stability,
                  d_measure = cluster.d_measure,
                  merge_type = cluster.merge_type,
                  rep = rep,
                  redundancy = cluster.redundancy,
                  identifiability = identifiability,
                  weakness = parameters.weakness)
  class(cluster) = "hclust" 
  
  return(cluster)
}

#### plot ####
plot.SMC <- function(x,
                     ...,
                     cluster) {
  # n - parameters number 
  # m - sensitivity matrix
  # labels - parameters names 
  # delta - delta
  
  plotDendogram.SMC(x, cluster=cluster, fig=c(0,0.75,0.4,1))
  
  plotEV.SMC(x, cluster=cluster, fig1=c(0.7,1,0.2,0.6), fig2=c(0.7,1,0,0.4), new=TRUE)
  
  plotFI.SMC(x, cluster=cluster, fig=c(0.55,1,0.45,1), new=TRUE)  
  
  barplot.SMC(x,cluster=cluster, fig=c(0,0.75,0,0.5), new=TRUE)
  
}

#### barplot ####
barplot.SMC <- function(x,
                        ...,
                        cluster,
                        fig      = c(0,1,0,1),
                        new      = FALSE,
                        order.cluster  = TRUE,
                        ident.cluster  = TRUE,
                        labels.args  = list()
                        
                        ){
  par(fig=fig, new=new)
  order <- cluster$order
  if( !order.cluster ){
    order <- 1:length(cluster$order)
  }
  
  labels.args.default <- list(las = 2,
                              cex.names  = 1)
  if(length(labels.args) != 0){
    labels.args.default <- c( labels.args[!(c("height", "names.args") %in% names(labels.args))],
                             labels.args.default[!(names(labels.args.default) %in% names(labels.args))])
  }
  FIM  <- x$FIM
  labels <- x$labels
  if( ! ident.cluster ){
    FIM  <- x$FIM_all
    labels <- x$labels_all
    order  <- 1:length(labels) 
  }
  
  ylim <- log(diag(FIM)[order]) 
  
  MIN  <- 2*min(ylim[ylim != -Inf])
  ylim <- sapply(ylim, function(y){ifelse(y != -Inf, y, MIN)})  
  do.call(barplot, 
          c(list(height = ylim,
                 names.arg = labels[order]),
            labels.args.default))
}

#### Sensitivity analysis ####
sa.SMC <- function(x,
                   ...,
                   cluster,
                   order.cluster = TRUE,
                   ident.cluster = TRUE
){
  
  labels <- x$labels
  names  <- x$names
  FIM  <- x$FIM
  order <- cluster$order
  
  colors <- rep(2,length(labels))  
  colors[cluster$rep] <- 1
  
  if( !ident.cluster ){
    labels <- x$labels_all
    names  <- x$names_all
    FIM  <- x$FIM_all
    order.cluster  <- FALSE 
    colors <- rep(2, length(labels))
    colors[match(x$labels, x$labels_all)][cluster$rep] <- 1
  }
  
  if( !order.cluster ){
    order <- 1:length(labels)
  }
  
  return(data.frame(
             ids = labels[order],
             names = names[order],
             sens = round(log(diag(FIM)[order]), digits = 4),
             identifiability = ifelse(colors[order] == 1, TRUE, FALSE)))
}

#### plotEV ####
plotEV.SMC <- function(x,
                          ...,
                          cluster,
                          fig1 = c(0.7,1,0.2,0.6), 
                          fig2 = c(0.7,1,0,0.4),
                          new  = FALSE
                          ){  
  n <- length(x$labels)
  
  par(fig=fig1, new=new)		
  matplot((1:n),(svd((cov2cor(x$FIM)))$d),xlab='',ylab='cor EV',pch=20,cex=1,lwd=1) #,xaxt="n")
  
  
  par(fig=fig2, new=new)		
  matplot((1:n),log(svd(((x$FIM)))$d),xlab='',ylab='cov EV',pch=20,cex=1,lwd=1) #,xaxt="n")
}

#### plotEV.Cor ####
plotEV.Cor.SMC <- function(x,
                          ...,
                          cluster,
                          fig = c(0.7,1,0.2,0.6), 
                          new  = FALSE
){  
  n <- length(x$labels)

  par(fig=fig, new=new)  	
  matplot((1:n),(svd((cov2cor(x$FIM)))$d),xlab='',ylab='cor EV',pch=20,cex=1,lwd=1) #,xaxt="n")
}

#### plotEV.Cov ####
plotEV.Cov.SMC <- function(x,
                       ...,
                       cluster,
                       fig = c(0.7,1,0,0.4),
                       new  = FALSE
){  
  n <- length(x$labels)
  
  par(fig=fig, new=new)		
  matplot((1:n),log(svd(((x$FIM)))$d),xlab='',ylab='cov EV',pch=20,cex=1,lwd=1) #,xaxt="n")
}

#### plotDendrogram ####
plotDendogram.SMC <- function(x,
                              
                              ...,
                              cluster,
                              fig      = c(0,1,0,1),
                              new      = TRUE,
                              labels.args  = list()
                              ){
  
  change_colors <- function(n, dend, num, merges, repr, merge_type, turn) {
    
    attr(dend,'edgePar')=list(col = 1 , lty=1, lwd= 1);
    
    if(is.leaf(dend)) {    
      dend
    }
    else {# holds num > n  
      s1 = merges[num - n, 1]
      if(s1 < 0) s1 = -s1 # leaf
      else s1 = n + s1 # normal node, its indice in merges array is shifted by n
      s2 = merges[num - n, 2]
      if(s2 < 0) s2 = -s2
      else s2 = n + s2
      
      if(repr[s1] == repr[num]) {
        dend[[1]] = change_colors(n, dend[[1]], s1, merges, repr, merge_type, turn)
        dend[[2]] = change_colors(n, dend[[2]], s2, merges, repr, merge_type, 1 - turn)
      }
      else {
        dend[[1]] = change_colors(n, dend[[1]], s1, merges, repr, merge_type, 1 - turn)
        dend[[2]] = change_colors(n, dend[[2]], s2, merges, repr, merge_type, turn)
      }
      
      
      if(num > 2*n - 1) line_style = 2
      else line_style = 1
      
      if (s1 > 2*n -1) col1 = 8
      else col1 = attr(dend[[1]],'edgePar')$col
      
      if (s2 > 2*n -1) col2 = 8
      else col2 = attr(dend[[2]],'edgePar')$col
      
      attr(dend[[1]],'edgePar')=list(col = col1, lty= line_style, lwd=2 + 2 * merge_type[num - n])
      attr(dend[[2]],'edgePar')=list(col = col2, lty= line_style, lwd=2 + 2 * merge_type[num - n])
      
      dend
    }
  }
  n <- length(x$labels)
  
  dend <- as.dendrogram(cluster)
  dend <- change_colors(n,
                        dend,
                        2*n-1,
                        cluster$merge,
                        cluster$representant,
                        cluster$merge_type)
  
  
  merges <- 1:(n-1)
  
  par(fig=fig, new=new)
  colors <- rep(2,length(x$labels))	
  colors[cluster$rep] <- 1
  plot(dend, leaflab= 'none')
  
  title(paste("identifiable parameters", as.character(sum(colors==1)),sep = " "))
  
  labels.args.default <- list(side = 1,
                              at = 1:length(cluster$order),
                              col = colors[cluster$order],
                              las = 2,
                              cex = 1) 
  if(length(labels.args) != 0){
    labels.args.default <- c(labels.args[!(c("text") %in% names(labels.args))],
                             labels.args.default[!(names(labels.args.default) %in% names(labels.args))])
  }
  do.call(mtext, c(list(text = x$labels[cluster$order]), 
                   labels.args.default
                   ))
  #mtext(x$labels[cluster$order], side = 1, at = 1:length(cluster$order), col = colors[cluster$order], las = 2, cex = 0.5)
  
}  

#### plot Fisher Information  ####
plotFI.SMC <- function(x,
                       ...,
                       cluster,
                       fig      = c(0,1,0,1),
                       new      = FALSE,
                       order.cluster  =  TRUE
                       ){
  par(fig=fig, new=new)
  plot(im((abs(cov2cor(x$FIM)))[cluster$order,cluster$order]), main="Fisher Information Matrix")
}  

#### redundancy ####
redundancy.SMC <- function(x,
                           ...,
                           cluster,
                           order.cluster  =  TRUE
){
#   get_correlations_cancor <- function(m, c1, c2) {
#     ccor = cancor(m[,c1], m[,c2], FALSE, FALSE)
#     ccor$cor
#   }
  order <- cluster$order
  if( !order.cluster ){
    order <- 1:length(x$labels_all) 
  }
  sapply(1:length(order), function(i){
    x$CANCOR$get_correlations_cancor(x$FIM, order[i], order[-i])
  })
}

