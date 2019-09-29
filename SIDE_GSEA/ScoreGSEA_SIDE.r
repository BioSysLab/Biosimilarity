ScoreGSEA_SIDE <-
  function(MergingSet,SignatureLength,ScoringDistance=c("avg", "max"),thr,p.value=F){
    
    #compute the distance of ranked PRL
    PRL.Distance_SIDE <-function(PRL,SignatureLength,ScoringDistance=c("avg", "max"))
    {
      ScoringDistance=match.arg(ScoringDistance,c("avg","max"))
      ES=matrix(1) 
      if (ncol(PRL)<=1){
        stop("the number of column of argument PRL must be greater than 1")
      }
      KBR=PRL[,1]
      PRLcol=ncol(PRL)
      for (i in 2:PRLcol){
        print(i) # check how it advances
        newPRL=as.matrix(PRL[,i],ncol=1,nrow=nrow(PRL))
        rownames(newPRL) <- rownames(PRL)
        d=integratePRL_SIDE(ES,KBR,newPRL,SignatureLength,ScoringDistance,thr)
        ES=d[[2]]
        KBR=d[[1]]
      }
      return((d[[3]]))
    }
    
    #compute the distance between two random permutations(RP)
    RP.Distance_SIDE <- function(RP1,RP2,SignatureLength,ScoringDistance=c("avg", "max"))
    {
      up1=which(RP1<=SignatureLength)
      down1=which(RP1>=length(RP1)-SignatureLength+1)
      up2=which(RP2<=SignatureLength)
      down2=which(RP2>=length(RP2)-SignatureLength+1)
      ES1=quickenrichmentscore_SIDE(up1,down1,RP2)
      ES2=quickenrichmentscore_SIDE(up2,down2,RP1)
      TES1=1-ES1
      TES2=1-ES2
      if (ScoringDistance=="avg")
        distances = (TES1+TES2)/2
      else
        distances= min(TES1+TES2)/2 
    } 
    
    ScoringDistance=match.arg(ScoringDistance,c("avg","max"))
    PRL=as.matrix(MergingSet)
    #rownames(PRL) <- rownames(MergingSet)
    if (is.null(ncol(PRL))==TRUE){
      stop("the class of argument MergingSet is incorrect")
    }
    if (ncol(PRL)<=1){
      stop("the number of samples of argument MergingSet must be greater than 1")
    }
    dis=PRL.Distance_SIDE(PRL,SignatureLength,ScoringDistance)
    if (p.value){
      p.results=matrix(0,nrow(dis),ncol(dis))
      prb_num=nrow(PRL)
      RPdis=c()
      RP1=sample(prb_num);
      for (i in 1:10000){
        RP2=sample(prb_num);
        RPdis[i]=RP.Distance_SIDE(RP1,RP2,SignatureLength,ScoringDistance)
      }
      RPdis=sort(RPdis)
      for (i in 1:nrow(dis)){
        for (j in 1:ncol(dis)){
          p.results[i,j]=(length(which(RPdis<=dis[i,j]))+1)/10001
        }
      }
      return(list(distances=dis,p.results=p.results))
    }
    return(dis)
  }
