quickenrichmentscore_SIDE <-
  function(S,S1,List,thr){
    ## S=UP,S1=DOWN,List=COMP
    
    ## a matrix is build according to the distances that are calculated from the PPI embeddings
    ## regarding the up-regulated factors and the list to be compared with
    distances_up <- matrix(0,nrow=nrow(List),ncol=nrow(S))
    for (c in 1:nrow(S)){
      for (r in 1:nrow(List)){
        distances_up[r,c] <- ppi_dist[as.character(List[r,1]),as.character(S[c,1])]
        if ((distances_up[r,c]<=thr)&(distances_up[r,c]>0)){
          distances_up[r,c] <- 0.8
        } else if (distances_up[r,c]==0){
          distances_up[r,c] <- 1
        } else {
          distances_up[r,c] <- 0
        }
      }
    }
    colnames(distances_up) <- S[,1]
    
    ## a matrix is build according to the distances that are calculated from the PPI embeddings
    ## regarding the down-regulated factors and the list to be compared with
    distances_down <- matrix(0,nrow=nrow(List),ncol=nrow(S1))
    for (c in 1:nrow(S1)){
      for (r in 1:nrow(List)){
        distances_down[r,c] <- ppi_dist[as.character(List[r,1]),as.character(S1[c,1])]
        if ((distances_down[r,c]<=thr)&(distances_down[r,c]>0)){
          distances_down[r,c] <- 0.8
        } else if (distances_down[r,c]==0){
          distances_down[r,c] <- 1
        } else {
          distances_down[r,c] <- 0
        }
      }
    }
    colnames(distances_down) <- S1[,1]
    
    ## calculation of the up-regulated factors' enrichment score
    ES_col <- matrix(0,nrow=nrow(List),ncol=1)
    for (i in 1:nrow(ES_col)){
      if (1%in%distances_up[i,]){
        ES_col[i,1] <- 1
        next()
      } else if (sum(distances_up[i,])>=0.8){
        #t <- sum(distances_up[i,])/0.8
        #ES_col[i,1] <- 0.8+0.2*(t-1)/(t+1)
        ES_col[i,1] <- 1
      }
    }
    
    N <- nrow(List)
    Nh <- nrow(S)
    NR1 <- length(which(ES_col[,1][ES_col[,1]>0]<1)) #numbers between 0.8-1
    NS <- length(which(ES_col[,1]==1)) 
    ES_hit <- ES_col
    for (h in 1:nrow(ES_hit)){
      if (ES_col[h,1]==1){
        ES_hit[h,1] <- ES_col[h,1]/NS
      } else if (ES_col[h]==0){
        ES_hit[h] <- 0
      } #else {
        #ES_hit[h] <- 1/(N/(NR+NS)) #ES_col[h]/(N-NR-NS)
      #}
    }
    
    ## in ES_miss we don't want e.g thr for the similars, so 
    ## it will be 0 as well
    ES_miss <- round(1-ES_col)
    for (h in 1:nrow(ES_miss)){ #NEW 
      if (ES_col[h]==0){#if (ES_miss[h]==1){
        ES_miss[h] <- ES_miss[h]/(N-NS)
      } else if (ES_col[h]==1){#(ES_miss[h]==0){
        ES_miss[h] <- 0
      } #else {
        #ES_miss[h] <- -(NR/N) ##test
      #}
    }
    
    hitc <- cumsum(ES_hit) 
    missc <- cumsum(ES_miss)
    Phit_n <- hitc #NEW
    Pmiss_n <- missc
    
    abshm_n <- abs(Phit_n-Pmiss_n)
    abshm_n <- as.matrix(abshm_n)
    t_n <- apply(abshm_n,2,which.max)
    ESUP_n <- Phit_n[t_n]-Pmiss_n[t_n]
    
    RS <- Phit_n-Pmiss_n
    
    ## calculation of the down-regulated factors' enrichment score
    ES_col <- matrix(0,nrow=nrow(List),ncol=1)
    for (i in 1:nrow(ES_col)){
      if (1%in%distances_down[i,]){
        ES_col[i,1] <- 1
        next()
      } else if (sum(distances_down[i,])>=0.8){
        #t <- sum(distances_down[i,])/0.8
        #ES_col[i,1] <- 0.8+0.2*(t-1)/(t+1)
        ES_col[i,1] <- 1
      }
    }
    
    N <- nrow(List)
    Nh <- nrow(S)
    NR2 <- length(which(ES_col[,1][ES_col[,1]>0]<1)) #numbers between 0.8-1 
    NS <- length(which(ES_col[,1]==1)) 
    ES_hit <- ES_col
    for (h in 1:nrow(ES_hit)){ 
      if (ES_col[h,1]==1){
        ES_hit[h,1] <- ES_col[h,1]/NS
      } else if (ES_col[h]==0){
        ES_hit[h] <- 0
      } #else {
        #ES_hit[h] <- 1/(N/(NR+NS)) #ES_col[h]/(N-NR-NS)
      #}
    }
    
    ## in ES_miss we don't want e.g thr for the similars, so 
    ## it will be 0 as well
    ES_miss <- round(1-ES_col)
    for (h in 1:nrow(ES_miss)){ #NEW 
      if (ES_col[h]==0){#if (ES_miss[h]==1){
        ES_miss[h] <- ES_miss[h]/(N-NS)
      } else if (ES_col[h]==1){#(ES_miss[h]==0){
        ES_miss[h] <- 0
      } #else {
        #ES_miss[h] <- -(NR/N) ##test
      #}
    }
    
    hitc <- cumsum(ES_hit) 
    missc <- cumsum(ES_miss)
    Phit_n <- hitc 
    Pmiss_n <- missc
    
    abshm_n <- abs(Phit_n-Pmiss_n)
    abshm_n <- as.matrix(abshm_n)
    t_n <- apply(abshm_n,2,which.max)
    ESDOWN_n <- Phit_n[t_n]-Pmiss_n[t_n]
    
    ## total enrichment score
    ES <- (ESUP_n-ESDOWN_n)/2
  }
    
