kappa_cluster<-function(x,deg=NULL,useTerm=FALSE,cutoff=0.5,overlap=0.5,minSize=5,escore=3){
    if(isTRUE(useTerm)){
        rownames(x) <- x$Term
    }else{
        rownames(x) <- x$Annot
    }
    if(is.null(deg)){
        deg <- unique(unlist(sapply(x$GeneID,function(x)split_(x))))
    }
    mat<-expand.grid(rownames(x),rownames(x),stringsAsFactors=F)
    mat <- mat[mat$Var1>mat$Var2,]
    mat$kappa <- apply(mat,1,function(y)kappa(x[y[1],"GeneID"],x[y[2],"GeneID"],deg))
    mat<-mat[mat$kappa > cutoff,]
    ml1<-split(mat$Var1,mat$Var2)
    ml2<-split(mat$Var2,mat$Var1)
    ###
    ml<-sapply(union(names(ml1),names(ml2)),function(x)c(ml1[[x]],ml2[[x]]))
    ###
    res<-merge_term(ml,overlap)
    #### remove the class with smaller number 
    res<-res[unlist(lapply(res, function(x)length(x)>=minSize))]
    ###find the smallest p value
    idx<-lapply(res,function(y)which.min(x[y,"Pvalue"]))
    ###use the term with smallest p value as name
    names(res)<-unlist(lapply(names(idx),function(x)res[[x]][idx[[x]]]))
    res<-sapply(unique(names(res)),function(x)unique(unlist(res[names(res)==x])))
    #####
    es <- lapply(res, function(y)calculate_Enrichment_Score(y,x))
    tl <- unlist(lapply(res,length))
    es <- unlist(es)/tl
    es <- es[es>escore]
    dx <- df[names(es),]
    tl <- tl[names(es)]
    dy <-data.frame(AnnotationCluster=1:length(es),EnrichmentScore=es,filteredClusterSize=tl)
    rs<-cbind(dy,dx)
    rs<-rs[order(rs$EnrichmentScore,decreasing = T),]
    return(rs)
}

merge_term<-function(x,overlap){
    ml <- x
    res<-list();
    for(i in names(ml)){
        lhs <- setdiff(names(ml),i)
        for(j in lhs){
            ov<-intersect(ml[[i]],ml[[j]])
            un<-union(ml[[i]],ml[[j]])
            ovl<-length(ov)/length(un)
            if(ovl > overlap){
                res[[i]]<-c(i,un)
                ml <- ml[setdiff(names(ml),j)]
            }else{
                res[[i]]<-c(i,ml[[i]])
            }
        }
    }
    return(res)
}

calculate_Enrichment_Score<-function(x,df){
    pvalue <- df[x,"Pvalue"]
    esp = ifelse(pvalue==0,16,-log10(pvalue))
    es = sum(esp);
}

kappa<-function(x,y,geneall){
    x<-unlist(strsplit(x,","))
    y<-unlist(strsplit(y,","))
    if(length(intersect(x,y))==0){
        kab=0
    }else{
        tmp<-matrix(0,2,2)
        tmp[1,1]<-length(intersect(x,y))
        tmp[2,1]<-length(setdiff(x,y))
        tmp[1,2]<-length(setdiff(y,x))
        tmp[2,2]<-length(setdiff(geneall,union(x,y)))
        oab<-(tmp[1,1]+tmp[2,2])/sum(tmp)
        aab<-((tmp[1,1]+tmp[2,1])*(tmp[1,1]+tmp[1,2])+(tmp[1,2]+tmp[2,2])*(tmp[2,1]+tmp[2,2]))/(sum(tmp)*sum(tmp))
        if(aab==1){
            kab=0
        }else{
            kab<-(oab-aab)/(1-aab)
        }
    }
    return(kab)
}
####
df<-read.delim("TableGO.txt",sep="\t",row.names=1)
deg<-read.delim("DEG.txt",sep="\t")
deg<-deg[,1]
res<-kappa_cluster(df,deg=deg)

