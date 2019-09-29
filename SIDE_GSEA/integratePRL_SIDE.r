integratePRL_SIDE <-
  function(ES,PRL,newPRL,SignatureLength,ScoringDistance=c("avg", "max"),thr){
    
    ScoringDistance=match.arg(ScoringDistance,c("avg","max"))
    origpd=rownames(PRL)
    newpd=rownames(newPRL)
    #PRL=exprs(PRL)
    #prpd=new("AnnotatedDataFrame",data=rbind(as(origpd,"data.frame"),as(newpd,"data.frame")))
    #newPRL=exprs(newPRL)
    PRL=cbind(PRL,as.matrix(newPRL))
    PRL=as.matrix(PRL)
    nelement=ncol(PRL)
    if (length(newPRL)>0) {
      n=length(newPRL);
      DOWN <- matrix(0,ncol=1,nrow=SignatureLength)
      
      for (i in 1:nrow(DOWN)){
        DOWN[i,1] <- rownames(newPRL)[which(newPRL[,1]==i)]
      }
      
      UP <- matrix(0,ncol=1,nrow=SignatureLength)
      
      for (i in 1:nrow(UP)){
        UP[i,1] <- rownames(newPRL)[which(newPRL[,1]==(nrow(newPRL)+1-i))]
      }
      
      ESrow=matrix(0,1,nelement)
      EScol=matrix(0,nelement,1)
      for (e in 1:(nelement-1)){
        COMP <- matrix(0,ncol=1,nrow=nrow(PRL))
        for (k in 1:nrow(COMP)){
          COMP[k,1] <- rownames(PRL)[which(PRL[,e]==(nrow(COMP)+1-k))]
        }
        #COMP <- rownames(PRL[,e])
        ESrow[e]=quickenrichmentscore_SIDE(UP,DOWN,COMP,thr)#PRL[,i])
        brgPRL=PRL[,e] #PRL[,i]
        brgPRL=as.matrix(brgPRL)
        rownames(brgPRL) <- rownames(PRL)
        n=nrow(brgPRL);
        
        down <- matrix(0,ncol=1,nrow=SignatureLength)
        
        for (i in 1:nrow(down)){
          down[i,1] <- rownames(brgPRL)[which(brgPRL==i)]
        }
        
        up <- matrix(0,ncol=1,nrow=SignatureLength)
        
        for (i in 1:nrow(up)){
          up[i,1] <- rownames(brgPRL)[which(brgPRL==(nrow(brgPRL)+1-i))]
        }
        comp <- matrix(0,ncol=1,nrow=nrow(newPRL))
        
        for (k in 1:nrow(comp)){
          comp[k,1] <- rownames(newPRL)[which(newPRL==(nrow(comp)+1-k))]
        }
        
        EScol[e]=quickenrichmentscore_SIDE(up,down,comp,thr)#(up,down,comp)
      }
      EScol[nelement] <- 1
      ESrow[nelement] <- 1
      ES=cbind(ES,EScol[c(1:(length(EScol)-1))])
      ES=rbind(ES,ESrow)#[c(1:(length(ESrow)-1))])
    }  
    ES=as.matrix(ES)
    distances=ES
    distances=as.matrix(distances)
    if (ScoringDistance=="avg"){
      distances = (distances+t(distances))/2
    }else{
      distances = pmax(distances,t(distances))/2
    }
    distances = 1-distances
    #colnames(ES)=rownames(rbind(as(origpd,"data.frame"),as(newpd,"data.frame")))
    colnames(ES)=rownames(origpd)
    colnames(PRL)=colnames(ES)
    colnames(distances)=colnames(PRL)
    rownames(distances)=colnames(ES)
    rownames(ES)=colnames(ES)
    #PRL=new("ExpressionSet",exprs=PRL,phenoData=prpd)
    list(PRL,ES,distances)
  }
