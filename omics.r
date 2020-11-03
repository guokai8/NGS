#################################################################################
# Author:Kai Guo guokai8@gmail.com ##############################################
#
# Last modified: 2019-12-10 11:00######################################
#
# Filename: omics.r
#
##### Functions collection for omics data analysis  #############
##### version 0.1.3 #############################################
##### Kai Guo ###################################################
##################################################################
options(stringsAsFactors = F)
###check and install package needed
check.get.packages<-function(pkg){
  options(warn=-1)
  res<-character()
  need<-as.matrix(sapply(1:length(pkg),function(i)
  {
    
    if(require(pkg[i],character.only = TRUE)==FALSE)
    {
      res<-c(res,pkg[i])
    }
  }))
  need<-as.character(unlist(need[!is.null(need)]))
  if(length(need)>0)
  {
    
    x<-sapply(need,install.packages,dependencies = TRUE)
    lib.fun<-function(need){
      sapply(1:length(need), function(i){
        out<-tryCatch(library(need[i], character.only= TRUE), error=function(e){need[i]})
        if(all(out==need[i])){need[i]}
      })
    }
    out<-as.character(unlist(lib.fun(need)))
    #try bioconductor
    if(length(out)>0){
      cat(paste("Package not found, trying Bioconductor..."),"\n")
      if (!requireNamespace("BiocManager"))
      install.packages("BiocManager")
      lib.fun.bioc<-function(need){
        sapply(1:length(need), function(i){
            tryCatch(BiocManager::install(need[i],ask=FAlSE),
                    error=function(e){need[i]})
        })
      }
      tmp<-lib.fun.bioc(out)
      final<-lib.fun(tmp)
      if(length(final)>0){cat(paste("could not find package: ",paste(as.character(unlist(final)),collapse=", "),sep=""),"\n")}
    }
  }
}
##########################################
check.get.packages(c("mixOmics","FactoMineR","factoextra","tidyr","dplyr",
"broom","pheatmap","ggrepel","ggpubr"))
##########################################
###load all needed library
library(tidyverse)
library(mixOmics)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(pheatmap)
###########################################
###input t.test and pre-normalized methods created by Kai
#####normalized functions
##1. Prepare dataset and normalized
####get CV ##### assume you got Testing Pool
getCV <- function(df){
  suppressMessages(require(tidyverse))
  TestPool<-df%>%dplyr::select(contains("Test"))
  dd<-df[,setdiff(colnames(df),colnames(TestPool))]
  dd$CV<-apply(TestPool, 1,sd)/apply(TestPool, 1, mean)*100
  return(dd)
}
####replace NA with min/2
##
replacezero <- function(x) "[<-"(x, !x | is.na(x), min(x[x > 0],
na.rm = TRUE) / 2)
###get Internal Standard if it exist
getIS<-function(df){
  is.dat<-df[grep('IS',df[,1]),]
  return(is.dat)
}
##remove Internal Standard
removeIS<-function(df){
  dd<-df[grep('IS',df[,1],invert = T),]
  return(dd)
}
###separate Lipid name
separateLipidInfo <- function(df,sep1=";",sep2=" ",sep3=":"){
  dd<-removeIS(df)
  colnames(dd)[1]<-"Sample"
  dx<-dd%>%separate(Sample,c("sam","un"),sep=sep1)%>%
    separate(sam,c("Class","sam1"),sep=sep2)%>%
    separate(sam1,c("Carbon","Double_bond"),sep=sep3)
  dx<-dx[,c(1:3,5:ncol(dx))]
  return(dx)
}
###separate IS lipid name
separateISLipid<-function(df,sep1=";",sep2=" ",sep3=":"){
  colnames(df)[1]<-"Sample"
  dx<-df%>%separate(Sample,c("sam","un"),sep=sep1)%>%
    separate(sam,c("IS","Class","sam1"),sep=sep2)%>%
    separate(sam1,c("Carbon","Double_bond"),sep=sep3)
  dx<-dx[,c(2:4,6:ncol(dx))]
  return(dx)
}
###quantile normalize 
qnormalize<-function(df){
  require(preprocessCore)
  dd<-normalize.quantiles(as.matrix(df[,2:ncol(df)]),copy = F)
  dd<-data.frame('Sample'= df$Sample, dd)
  colnames(dd)<-colnames(df)
  return(dd)
}
#####normalize with Internal Standard if there are IS 
normalizeInternalStandards <- function(df,is){
  is <- is[order(is$CV),] 
  # Internal standard normalization
  # If there's no corresponding IS, the IS having smallest CV used
  NewLipidData <- data.frame()
  #row.names(NewLipidData) <- row.names(new_df)
  for(i in 1:length(df$Class))
  { 
    lipidClass <- as.character(df$Class[i]) # lipid class
    for(j in 1:length(is$Class))
    {
      isClass <- as.character(is$Class[j])# IS class
      # Compare between IS and lipid class
      if(lipidClass == isClass)
      {
        is_data <- is[1, 4:(ncol(df)-1)] ####to decide if you matched class or just use loweset CV
        lipid_data <- df[i,4:(ncol(df)-1)]
        new_data <- lipid_data/is_data
        break
      }
      else
      {
        #if no matched class
        #Divide by the IS having the lowest CV 
        #print(paste("no matching",j, isClass))
        is_data <- is[1, 4:(ncol(df)-1)]
        lipid_data <- df[i,4:(ncol(df)-1)]
        new_data <- lipid_data/is_data
      }
      
    }
    processed_lipid <- cbind(df[i,1:3], new_data)
    NewLipidData <- rbind(NewLipidData,processed_lipid)
    #print(NewLipidData)
  }
  NewLipidData$CV<-df$CV
  return(NewLipidData)
}
###remove CV>=30% if you have Test pool for calculate CV
removeCV<-function(df){
  dx<-df%>%dplyr::filter(CV<=30)
  return(dx)
}
###combine postive and negative lipids
combinePosNeg<-function(ld,posidx,negidx,newnames){
  require(purrr)
  ll<-map2(ld[posidx],ld[negidx],function(x,y)rbind(x,y))
  names(ll)<-newnames
  return(ll)
}
###remove redundant lipids--use mean value or choose low CV
RemoveDuplicates<-function(df,useMean=F){
  if(useMean==FALSE){
  dx<-df%>%mutate(lipid=paste(Class,Carbon,Double_bond,sep=""))%>%
    group_by(lipid)%>%arrange(CV)%>%slice(1)
  dx<-dx[,1:(ncol(dx)-1)]
  }else{
    dx<-df%>%mutate(lipid=paste(Class,Carbon,Double_bond,sep="_"))%>%
      dplyr::select(4:(ncol(df)+1))%>%gather(sample,val,-lipid)%>%group_by(sample,lipid)%>%
      summarise(mu=mean(val))%>%spread(sample,mu)%>%
      separate(lipid,c("Class","Carbon","Double_bond"),sep="_")
  }
  return(dx)
}
###Remove Odd Chain 
RemoveOddChain <- function(df){
    ###some lipid don't have Carbon information
  even_chain_idx <- which((as.numeric(as.matrix(df$Carbon)) %% 2 == 0)|is.na(as.numeric(df$Carbon)))
  new_LipidData <- df[even_chain_idx,]
  cat("Removed: ",nrow(df)-nrow(new_LipidData),"species","\n")
  cat("Final species number is:",nrow(new_LipidData),"\n")
  return(as.data.frame(new_LipidData))
}
##prepare file for following analysis
prepare<-function(df){
  dd<-df[,4:(ncol(df)-1)]
  ###some lipid don't get Double_bond information
  df$Double_bond<-ifelse(is.na(as.numeric(df$Double_bond)),'',df$Double_bond)
  rownames(dd)<-paste(df$Class,paste(df$Carbon,df$Double_bond,sep=":"))
  rownames(dd)<-sub(":$",'',rownames(dd))
  return(dd)
}
quickFilter<-function(lhs){
  #loadpkg()
  #lhs<-lapply(lhs, function(x)qnormalize(x))
  lhs.cv<-lapply(lhs, function(x)getCV(x))#correct
  lhs.is<-lapply(lhs.cv, function(x)removeIS(x))#correct
  is.dat<-lapply(lhs.cv, function(x)getIS(x))#correct
  lhs.sp<-lapply(lhs.is, function(x)separateLipidInfo(x))#correct
  is.sp<-lapply(is.dat, function(x)separateISLipid(x))#correct
  lhs.nor<-map2(lhs.sp,is.sp,function(x,y)normalizeInternalStandards(x,y))#correct
  lhs.rm<-lapply(lhs.nor,function(x)removeCV(x))#correct
  lhs.cb<-combinePosNeg(lhs.rm,c(2,4,6),c(1,3,5),newnames = c("Liver","Plasma","SCN"))#correct
  lhs.me<-lapply(lhs.cb, function(x)RemoveDuplicates(x))
  lhs.od<-lapply(lhs.me, function(x)RemoveOddChain(x))
  rhs<-lapply(lhs.od, function(x)prepare(x))
  return(rhs)
}


################Tissue and Strain comparision###
################################################
###add options to do wilcox
####assume colnames or group were separate by "_"
my.t.test<- function(...,wilcox=FALSE) {
if(isTRUE(wilcox)){
   obj<-try(wilcox.test(...), silent=TRUE)
}else{
    obj<-try(t.test(...), silent=TRUE)
}
    if (is(obj, "try-error")) return(1) else return(obj$p.value)
}
####assume colnames or group were separate by "_"
ttesto<-function(x,wilcox=FALSE){
    te<-unnest(x)
    te<-as.data.frame(te)    
    return(my.t.test(val~group,data=te,wilcox=wilcox))
}
### t test with variance check
### add options to do wilcox
##################
ttestc<-function(x,wilcox=FALSE){
    te<-unnest(x)
    te<-as.data.frame(te)
    e.var<-tryCatch(var.test(te$val~te$group)$p.value, error=function(e){1})
    if(is.nan(e.var)|is.na(e.var)){e.var<-1}
    if(e.var<=0.05){equal.var<-FALSE} else {equal.var <-TRUE}
    #return(e.var)
    if(isTRUE(wilcox)){
        obj<-my.t.test(val~group,data=te,wilcox=TRUE)
    }else{
        obj<-my.t.test(val~group,data=te,var.equal = equal.var,wilcox=FALSE)
    }  
    return(obj)
}
fc.base<-function(x,name){
    tmp<-x%>%unnest()%>%mutate(gg=factor(group,levels=name))%>%group_by(gg)%>%summarise(mu=mean(val,na.rm=T),.groups="drop")%>%dplyr::select(mu)%>%purrr::map(function(x)x[2]/x[1])
    return(unlist(tmp))
}
fc2.base<-function(x){
    tmp<-x%>%unnest()%>%group_by(group)%>%summarise(mu=mean(val,na.rm=T),.groups="drop")%>%dplyr::select(mu)%>%purrr::map(function(x)x[2]/x[1])
    return(unlist(tmp))
}
fc.log<-function(x,name){
  tmp<-x%>%unnest()%>%mutate(gg=factor(group,levels=name))%>%group_by(gg)%>%summarise(mu=mean(val,na.rm=T),.groups="drop")%>%dplyr::select(mu)%>%purrr::map(function(x)x[2]-x[1])
  return(unlist(tmp))
}
##########scale method###
paretoscale <- function(z) {
  rowmean <- apply(z, 1, mean) # row means
  rowsd <- apply(z, 1, sd)  # row standard deviation
  rowsqrtsd <- sqrt(rowsd) # sqrt of sd
  rv <- sweep(z, 1, rowmean,"-")  # mean center
  rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
  return(rv)
}

multigroup.t.test<-function(dd,group,sep="_",FDR="BH",na.omit=TRUE,ref=NULL,t.correct=FALSE,log=FALSE,base=2,wilcox=FALSE){
    cn<-combn(levels(as.factor(group)),2)
    colnames(dd)<-group
    names<-apply(cn,2,function(x)paste(x[2],x[1],sep="_vs_"))
    do.call("cbind",lapply(1:ncol(cn), function(x){
         dt<-dd[,colnames(dd)%in%cn[,x]]
         dg<-group[group%in%cn[,x]]
         dr<-multi.t.test(dt,group=dg,FDR=FDR,na.omit=na.omit,t.correct=t.correct,log=log,base=base,wilcox=wilcox)
         colnames(dr)<-paste(names[x],colnames(dr),sep="_")
         dr
    }))
}
multilevel.t.test<-function(data,group=NULL,FDR="BH",levelsep="_",na.omit=TRUE,t.correct=FALSE,log=FALSE,base=2,wilcox=FALSE,...){
  if(!is.null(group)){
    colnames(data)=paste(group,1:ncol(data),sep=levelsep)
  }
  tnest<-data%>%rownames_to_column(var="spe")%>%gather(sample,val,-spe)%>%
    separate(sample,c("str","group","id"),sep=levelsep)%>%dplyr::select(spe,group,str,val)%>%
    group_by(str,spe)%>%nest()
    if(t.correct==TRUE){
        ttest=ttestc
    }else{
        ttest=ttesto
    }
    if(log==TRUE){
        fc=fc.log
    }else{
        fc=fc.base
    }
  res<-tnest%>%mutate(pvalue=purrr::map(data,function(x)ttest(x,wilcox=wilcox)))%>%
    dplyr::select(str,spe,pvalue)%>%unnest()%>%spread(str,pvalue)%>%
    as.data.frame()%>%column_to_rownames(var="spe")
  respadj<-purrr::map_df(res,function(x)p.adjust(x,method="fdr"))
  colnames(respadj)<-paste(colnames(res),"padj",sep="_")
  colnames(res)<-paste(colnames(res),"pvalue",sep="_")
  res_fc<-tnest%>%mutate(FC=purrr::map(data,function(x)fc(x)))%>%dplyr::select(str,spe,FC)%>%unnest()%>%
    spread(str,FC)%>%as.data.frame()%>%column_to_rownames(var="spe")
    if(log==TRUE){
        res_fc2=log2(base^(res_fc))
        res_fc=base^(res_fc)
    }else{
        res_fc2=log2(res_fc)
    }
  #colnames(res_fc)<-paste(colnames(res_fc),c("FC","log2FC"),sep="_")
  colnames(res_fc2)<-paste(colnames(res_fc),"log2FC",sep="_")
  colnames(res_fc)<-paste(colnames(res_fc),"FC",sep="_")
  dd<-data.frame(res,respadj,res_fc,res_fc2)
  if(na.omit==TRUE){
      dd<-na.omit(dd)
  }
 return(dd)
}
multi.t.test<-function(data,group=NULL,FDR="BH",levelsep="_",na.omit=TRUE,t.correct=FALSE,log=FALSE,base=2,wilcox=FALSE,...){
    if(!is.null(group)){
        colnames(data)=paste(group,1:ncol(data),sep=levelsep)
    }
    name<-levels(as.factor(group))
    tnest<-data%>%rownames_to_column(var="spe")%>%gather(sample,val,-spe)%>%
    separate(sample,c("group","id"),sep="_")%>%dplyr::select(spe,group,val)%>%
    group_by(spe)%>%nest()
    if(t.correct==TRUE){
        ttest=ttestc
    }else{
        ttest=ttesto
    }
    if(log==TRUE){
        fc=fc.log
    }else{
        fc=fc.base
    }
    res<-tnest%>%mutate(pvalue=purrr::map(data,function(x)ttest(x,wilcox=wilcox)))%>%
    dplyr::select(spe,pvalue)%>%unnest()%>%as.data.frame()%>%
    column_to_rownames(var="spe")
    respadj<-purrr::map_df(res,function(x)p.adjust(x,method=FDR))
    colnames(respadj)<-"padj"
    #colnames(res)<-paste(colnames(res),"pvalue",sep="_")
    res_fc<-tnest%>%mutate(FC=purrr::map(data,function(x)fc(x,name)))%>%dplyr::select(spe,FC)%>%unnest()%>%
    as.data.frame()%>%column_to_rownames(var="spe")
    if(log==TRUE){
        res_fc[,2]=log2(base^(res_fc[,1]))
        res_fc[,1]=base^(res_fc[,1])
    }else{
        res_fc[,2]=log2(res_fc[,1])
    }
    colnames(res_fc)[2]<-"log2FC"
    #colnames(res_fc)<-paste(colnames(res_fc),c("FC","log2FC"),sep="_")
    dd<-data.frame(res,respadj,res_fc)
    if(na.omit==TRUE){
    dd<-na.omit(dd)
    }
    return(dd)
}
multi.fc<-function(data,group=NULL,log=TRUE,base=2,ref=NULL,...){
  if(!is.null(group)){
    colnames(data)=paste(group,1:ncol(data),sep="@")
  }
      if(log==TRUE){
        fc=fc.log
    }else{
        fc=fc.base
    }
    name<-levels(as.factor(group))
    if(!is.null(ref)){
        name <- c(ref,setdiff(name,ref))
    }
  tnest<-data%>%rownames_to_column(var="spe")%>%gather(sample,val,-spe)%>%
    separate(sample,c("group","id"),sep="@")%>%dplyr::select(spe,group,val)%>%
    group_by(spe)%>%nest()
    
  res_fc<-tnest%>%mutate(FC=purrr::map(data,function(x)fc(x,name)))%>%dplyr::select(spe,FC)%>%unnest()%>%
    as.data.frame()%>%column_to_rownames(var="spe")
        if(log==TRUE){
        res_fc[,2]=log2(base^(res_fc[,1]))
        res_fc[,1]=base^(res_fc[,1])
    }else{
        res_fc[,2]=log2(res_fc[,1])
    }
    colnames(res_fc)[2]<-"log2FC"
  return(res_fc)
}

##make figures
###########################################
###volcano plot for differential results  ###modified by Kai
volcano.plot<-function(FC=FC,p.value=p.value,labels=1:length(FC),
                       FC.lim=3,p.value.lim=.05,xsize=10,ysize=10,tsize=20,x.lab="Log2 Fold change",
                       y.lab="-log10 Pvalue",size=2,alpha=.5,high="red",mid="grey",low="skyblue",
                       topn=20,usepadj=FALSE,returnData=FALSE){
    #convert to log (non-log makes less sense and is harder to implement
    require(ggrepel)
    t.p.lim<--log10(p.value.lim)
    t.FC.lim<-log2(FC.lim)
    if(usepadj==TRUE){
        y.lab="-log10 Padj"
    }else{
        y.lab<-y.lab
    }
    x.lab<-x.lab
    FC<-log(FC,base=2)
    p.value<--log10(p.value)
    fact<-factor(FC>=t.FC.lim&p.value>t.p.lim|FC<(t.FC.lim*-1)&p.value>t.p.lim)
    labels[fact!=TRUE]<-""
    tmp<-data.frame(FC,p.value,labels,fact)
    tmp<-tmp[order(tmp$p.value,decreasing = T),]
    tmp$color<-ifelse(tmp$FC>0,"Up","Down")
    tmp$color[tmp$fact==FALSE]<-"No"
    if(isTRUE(returnData)){
      return(tmp%>%filter(fact==TRUE))
    }
    if(sum(na.omit(fact==TRUE))>=topn){
        tmp$labels[cumsum(tmp$fact==TRUE)>topn]<-""
        # tmp$labels[(topn+1):nrow(tmp)]<-""
    }
    # darkcyan
    cols=c("Up"=high,"No"=mid,"Down"=low)
    p<-ggplot(tmp, aes(x=FC,y=p.value,color=color)) + geom_point(alpha=alpha,size=size)+ geom_text_repel(aes(label=labels),color="black",size=size+1,show_guides=FALSE)+
        geom_vline(xintercept = t.FC.lim*c(-1,1),linetype=2,size=.75,color="red") +
        geom_hline(yintercept = t.p.lim,linetype=2,size=.75,color="red")+
        scale_color_manual(name="Significant", values=cols) +
        ylab(y.lab) + xlab(x.lab) +theme_light(base_size = 15) + guides(color = guide_legend(override.aes = list(size = 8)))+
        theme(axis.text.x=element_text(size=xsize),axis.text.y=element_text(size=ysize))+
        ggtitle(paste0(sum(na.omit(fact==TRUE))," metabolites selected"))+theme(plot.title = element_text(size = tsize, face = "bold"))
    p
}
##############################################
##Venn Diagram plot 
###venn plot  ###created by Kai
venn.plot<-function(rlist,n,...){
dev.new()
check.get.packages("VennDiagram")
mycol=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
ven.plot<-venn.diagram(x,filename = NULL,
col = "transparent",
fill = mycol[1:n],
alpha = 0.50,
cex=2,
cat.col = mycol[1:n],
cat.cex = 1.5,
cat.fontface = "bold",
margin = 0.05
)
grid.draw(ven.plot)
}
venn.plot2<-function(rlist,n,...){
    dev.new()
    check.get.packages("VennDiagram")
    mycol=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
    ven.plot<-venn.diagram(rlist,filename = NULL,
    col = "black",
    fill = mycol[1:n],
    alpha = 0.50,
    cex=2,
    cat.col = mycol[1:n],
    cat.cex = 1.5,
    cat.fontface = "bold",
    margin = 0.05
    )
    grid.draw(ven.plot)
}
### draw VennDiagram
###install by devtools::install_github("guokai8/VennDetail")
library(VennDetail)
ven.plot<-function(rlist,return=TRUE,cex=2,cat.cex = 1.5){
    venlist<-venndetail(rlist)
    plot(venlist,cex=cex,cat.cex=cat.cex)
    if(return==TRUE){
        return(venlist)
    }
}
###
pheat<-function(df,sel,group,fontsize=3,border="white"){
    annotation_col = data.frame(
    Group = group
    )
    rownames(annotation_col) = colnames(df)
    pheatmap(df[sel,],scale="row",color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
    annotation_col = annotation_col,fontsize_row = fontsize,border_color=border )
}
##############
###draw heatmap ###modified by Kai
draw.heatmap<-function(data, class.factor=NULL, class.color=NULL, heatmap.color = NULL, border.color=NULL, match.dim=2, type=c("none","z.scale", "spearman", "pearson","biweight"),
                       cluster.method = c("none","ward", "single", "complete", "average", "mcquitty", "median" , "centroid"),
                       distance.method = c("none","euclidean", "maximum", "manhattan", "canberra", "binary" ,"minkowski"),
                       alpha = NULL,font.size = 12,show.names=F, ncolors = 100, level.limit=10){
  check.get.packages("pheatmap")
  type<-match.arg(type)
  cluster.method<-match.arg(cluster.method)
  distance.method<-match.arg(distance.method)
  #prepare data object
  if(match.dim==1){tmp.data<-data.frame(t(data.frame(data)))} else { tmp.data<-data.frame(data)}
  
  if(type == "z.scale") {tmp.data<-scale(tmp.data, center=TRUE, scale = TRUE)}		
  # calculate correlations
  if(!type=="none"& !type == "z.scale" ){
    tmp<-calculate.correlations(tmp.data,type=type,results="matrix")
    tmp.data<-tmp$cor
    tmp.data.pvalue<-tmp$p.value
    
    if(is.numeric(alpha)){ # make discrete
      tmp.data[tmp.data>0]<-1
      tmp.data[tmp.data<0]<--1
      tmp.data[tmp.data.pvalue>alpha]<-0
      ncolors<-3 # limit heat map to 3 colors
    }
  }
  
  if(!cluster.method=="none"){
    cluster_rows<-cluster_cols<-T
    # calculate distances
    if(!distance.method=="none"){
      if(type=="none"|type == "z.scale" ){
        # if(match.dim==2){
        # clustering_distance_cols<-dist(tmp.data, method= distance.method)
        # clustering_distance_rows<-dist(data.frame(t(tmp.data)), method= distance.method)
        # } else {
        clustering_distance_rows<-dist(tmp.data, method= distance.method)
        clustering_distance_cols<-dist(data.frame(t(tmp.data)), method= distance.method)
        # }
      } else {                         				
        tmp.data.dist<-dist(tmp.data, method= distance.method) # for correlations ignore sign focus on magnitude
        clustering_distance_cols<-clustering_distance_rows<-tmp.data.dist
      }
    } 
  } else {
    cluster_rows<-cluster_cols<-F
    cluster.method<-NULL			
  }
  # set colors for top of the heatmap annotation
  if(!is.null(class.factor)){
    annotation<- data.frame(sapply(class.factor,as.factor))
    dimnames(annotation)<-dimnames(class.factor)
  } else {
    annotation<-NA
  }
  
  if(is.na(annotation)){
    annotation.color<-NA
  } else {
    if(is.null(class.color)){
      choices<-colors()[-agrep(c("gray","grey"),colors())]
      #test if factor is continuous or discrete
      fct.type<-sapply(1:ncol(annotation),function(i){
        obj<-as.factor(annotation[,i])
        nlevels(obj)>level.limit
      })
      
      annotation.color<-lapply(1:ncol(class.factor), function(i){
        if(fct.type[i]){
          tmp<-choices[sample(1:length(choices),1)]
          tmp<-colorRampPalette(c("white",tmp))(nlevels(annotation[,i]))
          names(tmp)<-levels(annotation[,i])
        } else {
          tmp<-choices[sample(1:length(choices),nlevels(annotation[,i]))]
          names(tmp)<-levels(annotation[,i])
        }
        tmp
      })
      names(annotation.color)<-colnames(annotation)
    } else {
      annotation.color<- lapply(1:ncol(class.color), function(i){
        fct<-as.factor(annotation[,i])
        sapply(1:nlevels(fct), function(j){
          name<-levels(fct)[j]
          tmp<-fixlc(class.color[which(fct==name)[1], i])
          names(tmp)<-name
          tmp
        })
      })
      names(annotation.color)<-colnames(annotation)
    }
  }
  if(show.names==FALSE){
    show_rownames<-show_colnames <-F
  } else { 
    show_rownames<-show_colnames <-T
  }
  #set colors
  heat.col<-function(col=1,n){ # color for heat map
    # cat("Choices  \n","1 = orange-white-blue",
    # "\n", "2 = white-red \n","3 = white-black \n",
    # "4 = green-black-red \n",
    # "5 = yellow-red \n",
    # "6 = white-yellow-red \n",
    # "7 = white-yellow-orange-red \n")
    col= switch(col,
                "1" = colorRampPalette(c("navy", "white", "orange"))(n),
                "2" = colorRampPalette(c( "white", "red"))(n),
                "3" = colorRampPalette(c( "white", "black"))(n),
                "4" = colorRampPalette(c( "green","black", "red"))(n),
                "5" = colorRampPalette(c( "yellow", "red"))(n),
                "6" = colorRampPalette(c( "white","yellow", "red"))(n),
                "7" = colorRampPalette(c( "white","yellow","orange", "red"))(n),
                "8" = colorRampPalette(c("navy", "white", "firebrick3"))(n))
    return(col)
  }
  if(is.null(heatmap.color)){ heat.color<-heat.col(8,ncolors) } else { heat.color<-colorRampPalette(heatmap.color)(ncolors)}
  #plot heat map
  pheatmap(tmp.data,col= heat.color,
           show_rownames = show_rownames,show_colnames = show_colnames, # show labels
           cluster_rows = cluster_rows, cluster_cols=cluster_cols, # cluster
           clustering_method=cluster.method,
           clustering_distance_rows = clustering_distance_rows,
           clustering_distance_cols = clustering_distance_cols,
           border_color = border.color,
           annotation = annotation, # annotation factor
           annotation_colors = annotation.color, # colors for each column of annotation factor
           annotation_legend= TRUE,
           fontsize = font.size)
}
#########################################################
##PCA part
####PCA part # created by Kai
make.pca<-function(data,scale=F){
  check.get.packages(c("FactoMineR","factoextra"))
  res.pca<-PCA(data,scale.unit = scale,graph=F)
  return(res.pca)
}
##plot the PCA results created by Kai
plot.pca<-function(pca.obj,type=c("scree","scores","contrib","biplot"),top=50,group=group,axes=1,choice="var",repel=TRUE){
  if(length(type)==4){
    type="scores"
  }
  if(length(type)>1){
    stop("you must select one type of plot in scree,scores,contrib,biplot \n")
  }
  if(!is.factor(group)){
   group<-factor(group)
  }
  switch(type,
         scree=.local<-function(pca.obj){
           fviz_screeplot(pca.obj)+theme(axis.text =element_text(size=15),axis.title = element_text(size=15))
         },
         scores=.local<-function(pca.obj){
           fviz_pca_ind(pca.obj,repel = repel,habillage = group,addEllipses = T,ellipse.level=0.95, palette = "Dark2")+
             theme(axis.text =element_text(size=15),axis.title = element_text(size=15))
         },
         contrib=.local<-function(pca.obj){
           .theme<-theme_bw(base_size = 15) %+replace% 
             theme(legend.background = element_blank(), 
                   legend.key = element_blank(), panel.background = element_blank(), 
                   panel.border = element_blank(), strip.background = element_blank(), 
                   plot.background = element_blank(), complete = TRUE,axis.text.x=element_text(angle=90,size=10)
                   )
           fviz_contrib(pca.obj,choice=choice,top = top,axes=axes)+.theme
         },
         biplot=.local<-function(pca.obj){
           fviz_pca_biplot(pca.obj,habillage = group,addEllipses = T,ellipse.level=0.95,palette = "Dark2",select.var=list(contrib = top))
         }
         )
  .local(pca.obj)
}
##plot the heatmap based on selected species created by Kai
##top should be how many species you want to display on the heatmpa based on contribution
pca.heatmap<-function(data,pca.obj,top=50,group=group,contrib=1,fontsize_row=6,fontsize_col=6,mycol=colorRampPalette(c("navy", "white", "firebrick3"))(78),border="lightgrey"){
  agroup<-data.frame(condition=group)
  rownames(agroup)<-rownames(data)
  pheatmap(t(data)[names(sort(pca.obj$var$contrib[,contrib],decreasing = T)[1:top]),],scale="row",
           annotation_col = agroup,cluster_cols = F,fontsize_row = fontsize_row,fontsize_col= fontsize_col,color=mycol,border=border)
}
####plot loadings value for PCA #created by Kai
pca.loadings.plot<-function(pca.obj,top=50,choice="var"){
  require(tidyr)
  if(choice=="var"){
    pt<-as.data.frame(pca.obj$var$coord)%>%rownames_to_column(var="sample")%>%arrange(desc(abs(Dim.1)))
    ploadings<-pt[1:top,1:3]
  }else{
    pt<-as.data.frame(pca.obj$ind$coord)%>%rownames_to_column(var="sample")%>%arrange(desc(abs(Dim.1)))
    ploadings<-pt[,1:3]
  }
  colnames(ploadings)[2:3]=c("PC1","PC2")
  ploadings<-as.data.frame(ploadings)
  ploadings%>%gather(PC,value,-sample)%>%ggplot(aes(sample,value))+geom_bar(stat="identity",position="dodge",aes(fill=PC))+
    theme(axis.text.x=element_text(angle=90))+xlab("Sample")+ylab("Loadings value")
}
#####################
###PLS-DA part
###make pls-da  ###created by Kai
make.PLSDA.model<-function(y,data,ncomp=NULL,scale=T,cv=TRUE, ref=NULL,
cv.method="Mfold",folds=10,nrepeat=10,log=TRUE,base=2,...){
  rr<-list()
  require(mixOmics)
  if(is.null(ncomp)){
      ncomp=10
    }
  res<-plsda(data,y,ncomp=ncomp,scale =scale)
  if(isTRUE(cv)){
     res.perf <- perf(res, validation = cv.method, folds = folds,
                    progressBar = FALSE, nrepeat = nrepeat)
     ncomp <- res.perf$choice.ncomp[2,1]
     res<-plsda(data,y,ncomp=ncomp,scale =scale)
  }else{
      res.perf=NULL
    }
      name<-levels(as.factor(y))
      if(!is.null(ref)){
          name<-c(ref,setdiff(name,ref))
        }
      if(length(name)==2){
      fc<-multi.fc(as.data.frame(t(data)),y,log=log,base=base,ref=ref)
      fc$dir<-ifelse(fc$log2FC>0,name[2],name[1])
      fc<-fc[rownames(res$loadings$X),]
}else{
     fc<-NULL
}
  rr$fc<-fc
  rr$plsda<-res
  rr$vip <- vip(res);
  rr$perf <- res.perf
  return(rr)
}
### plot PLS-DA VIP value created by Kai
plot.vip<-function(data,cutoff=1,top=50){
    data<-as.data.frame(data)
    data$species=rownames(data)
    dd<-data%>%filter(comp1>=cutoff)%>%arrange(desc(comp1))%>%head(top)
    dd$species=factor(dd$species,levels =dd$species )
    ggplot(dd,aes(species,comp1))+geom_bar(stat="identity")+
    theme_light(base_size = 15)+ylab("VIP")+xlab("")+
    theme(axis.text.x = element_text(angle=90))
    
}
#### plot PLS-DA value #created by Kai
plot.vip2<-function(data,cutoff=1,top=50,fc=NULL,label.color="white",label.size=8,dot.size=8,rotate=TRUE){
    require(ggpubr)
    data<-as.data.frame(data)
    data$species=rownames(data)
    if(!is.null(fc)){
        data$comp1=ifelse(fc$log2FC>0,data$comp1,-data$comp1)
        data$group<-fc$dir
    }
    dd<-data%>%filter(abs(comp1)>=cutoff)%>%arrange(desc(abs(comp1)))%>%head(top)
    dd$species=factor(dd$species,levels =dd$species )
    colnames(dd)[1]<-"VIP"
    #I("#00AFBB")
    #"#00AFBB""
    if(!is.null(fc)){
    ggdotchart(dd,x="species",y="VIP",add="segments",color="group",rotate=rotate,dot.size=dot.size,add.params = list(color = "group", size = 2),
    label=round(abs(dd$VIP),2),font.label = list(color = label.color, size = label.size,vjust=0.2))+xlab("")
    }else{
      ggdotchart(dd,x="species",y="VIP",add="segments",color=I("#00AFBB"),rotate=rotate,dot.size=dot.size,add.params = list(color = "#00AFBB", size = 2),
    label=round(abs(dd$VIP),2),font.label = list(color = label.color, size = label.size,vjust=0.2))+xlab("")   
    }
}

#### write PLS-DA VIP  out # created by Kai
write.vip<-function(data,name,cutoff=1){
    data<-as.data.frame(data)
    data$species=rownames(data)
    dd<-data%>%filter(comp1>=cutoff)%>%arrange(desc(comp1))
    write.table(dd,file=paste(name,"_VIP_cut.txt",sep=""),row.names=F)
}
###plot pls-da scores created by Kai
plot.pls.scores<-function(obj,label=NA,group=group,repel=F,ellipse=T){
  obj<-obj$plsda
  scores<-as.data.frame(obj$variates$X)
  colnames(scores)<-paste("LV",1:ncol(scores),sep="")
  scores$group=group
  if(is.na(label)){
    scores$label=rownames(scores)
  }else{
    scores$label=label
  }
  exp<-obj$explained_variance$X
  if(repel==TRUE){
    require(ggrepel)
    p<-ggplot(scores,aes(LV1,LV2,color=group))+geom_point()+xlab(paste("LV1 (",round(exp[1]*100,2),"%)",sep=""))+
      ylab(paste("LV2 (",round(exp[2]*100,digits = 2),"%)",sep=""))+
      geom_text_repel(aes(label=label),show.legend = F)+theme_minimal(base_size=15)
  }else{
    p<-ggplot(scores,aes(LV1,LV2,color=group))+geom_point()+xlab(paste("LV1 (",round(exp[1]*100,2),"%)",sep=""))+theme_minimal(base_size=15)+
      ylab(paste("LV2 (",round(exp[2]*100,digits = 2),"%)",sep=""))+
      geom_text(aes(label=label),show.legend = F)
  }
  if(ellipse==TRUE){
      print(p+stat_ellipse(aes(LV1,LV2,color=group),type="norm"))
  }else{
      print(p)
  }
}
##################### created by Kai
## plot cross-validation results
plot.perf<-function(res){
require(mixOmics)
plot(res$perf, overlay = 'measure', sd = TRUE,legend.position = "horizontal")
}

####auc plot for plsda
plot.auc <- function(res){
require(mixOmics)
auc.plsda<-auroc(res$plsda)
}

##################################################################
###OPLS-DA part
#
###OPLS-DA  ##modified by Kai ####pls.y need to be one column factor
make.OSC.PLS.model<-function(pls.y,pls.data,comp=2,OSC.comp=1,validation = "LOO",progress=TRUE,cv.scale=FALSE,return.obj="stats",train.test.index=NULL,OPLSDA=FALSE,...){
    
    check.get.packages("pls")
    #initialize
    OSC.results<-list()
    OSC.results$data[[1]]<-pls.data # may need to
    OSC.results$y[[1]]<-pls.y<-as.matrix(pls.y)
    if(!is.null(train.test.index)){ # objects fo predictions
        OSC.results$test.data[[1]]<-OSC.results$data[[1]][train.test.index=="test",]
        # if(cv.scale==TRUE){ # the same?
        # OSC.results$test.y<-test.y<-as.matrix(OSC.results$y[[1]][train.test.index=="test",])
        # } else {
        OSC.results$test.y<-test.y<-as.matrix(OSC.results$y[[1]][train.test.index=="test",])
        # }
        OSC.results$data[[1]]<-OSC.results$data[[1]][train.test.index=="train",] # data may need to be a data.frame
        OSC.results$y[[1]]<-as.matrix(OSC.results$y[[1]][train.test.index=="train",])
    }
    
    if (progress == TRUE){ pb <- txtProgressBar(min = 0, max = (OSC.comp+1), style = 3)}
    
    #need to iteratively fit models for each OSC
    for(i in 1:(OSC.comp+1)){
        data<-OSC.results$data[[i]]
        tmp.model<-plsr(OSC.results$y[[1]]~., data = data, ncomp = comp, validation = validation ,scale=cv.scale,...)
        ww<-tmp.model$loading.weights[,1] #
        pp<-tmp.model$loadings[,1]
        w.ortho<- pp - crossprod(ww,pp)/crossprod(ww)*ww
        t.ortho<- as.matrix(data) %*% w.ortho
        p.ortho<- crossprod(as.matrix(data),t.ortho)/ c(crossprod(t.ortho))
        Xcorr<- data - tcrossprod(t.ortho,p.ortho)
        
        
        #stats for classifiers, currently for two group comparisons only
        # for training data only
        if(OPLSDA==TRUE){
            pred.val<-as.data.frame(tmp.model$fitted.values)[,comp]
            OSC.results$OPLSDA.stats[[i]]<-O.PLS.DA.stats(pred=pred.val,truth=unlist(OSC.results$y[[1]])[,1])
        }
        
        #prediction objects
        if(!is.null(train.test.index)){
            test.data<-OSC.results$test.data[[i]]
            predicted.mod<-    predict(tmp.model,newdata=test.data, ncomp=1:comp, comps=1:comp, type="response")
            OSC.results$predicted.Y[[i]]<-predicted.mod
            # predicted.RMSEP
            OSC.results$predicted.RMSEP[[i]]<-sapply(1:ncol(test.y), function(i){
                (sum((predicted.mod[,i]-test.y[,i])^2)/nrow(predicted.mod))^.5
            })
            #stats for classifiers, currently for two group comparisons only
            if(OPLSDA==TRUE){
                OSC.results$OPLSDA.stats[[i]]<-O.PLS.DA.stats(pred=predicted.mod,truth=test.y)
            }
            
            t.tst<-as.matrix(test.data)%*%w.ortho
            p.tst <- crossprod(as.matrix(test.data), t.tst) / c(crossprod(t.tst))
            OSC.test.data <- as.matrix(test.data) - tcrossprod(t.tst, p.tst)
            OSC.results$test.data[[i+1]]<-OSC.test.data # for next round
        }
        
        #store results
        OSC.results$RMSEP[[i]]<-matrix(t(RMSEP(tmp.model)$val[dim(RMSEP(tmp.model)$val)[1],,]),,ncol=ncol(pls.y)) #
        OSC.results$rmsep[[i]]<- RMSEP(tmp.model)$val[dim(RMSEP(tmp.model)$val)[1],,comp+1]# CV adjusted rmsep for each y by column
        OSC.results$Q2[[i]]<-matrix(pls::R2(tmp.model)$val,ncol=ncol(pls.y),byrow=TRUE)
        OSC.results$Xvar[[i]]<-drop(tmp.model$Xvar/tmp.model$Xtotvar)#matrix(drop(tmp.model$Xvar/tmp.model$Xtotvar),ncol=1)
        OSC.results$fitted.values[[i]]<-tmp.model$fitted.values
        OSC.results$scores[[i]]<-tmp.model$scores
        OSC.results$loadings[[i]]<-tmp.model$loadings
        OSC.results$loading.weights[[i]]<-tmp.model$loading.weights
        OSC.results$total.LVs[[i]]<-comp
        OSC.results$OSC.LVs[[i]]<-i-1 # account for first model not having any OSC LVs
        #coefficients
        OSC.results$coefficients[[i]]<-matrix(coefficients(tmp.model),ncol=1)
        #initialize data for next round
        OSC.results$data[[i+1]]<-as.data.frame(Xcorr)
        OSC.results$model.description<-as.list( sys.call() )#as.list(environment())
        #update timer
        if (progress == TRUE){setTxtProgressBar(pb, i)}
        
        #exit loop if no more orthogonal dimensions can be calculated
        if(any(is.na(Xcorr))){break}
    }
  
  if (progress == TRUE){close(pb)}
  #calculate VIP if single Y and oscorespls
  #based on http://mevik.net/work/software/VIP.R 
  if(ncol(pls.y)==1){
    object<-tmp.model
    VIP<-function(object){
      SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
      Wnorm2 <- colSums(object$loading.weights^2)
      SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
      t(sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS)))[,1:comp,drop=FALSE]
    }
    OSC.results$VIP<-tryCatch(VIP(object),error=function(e){data.frame(VIP=matrix(1,nrow(tmp.model$loadings[,]),1))})	
  } else { OSC.results$VIP<-data.frame(VIP=matrix(1,nrow(tmp.model$loadings[,]),ncol(pls.y)))}
  
  
  if (return.obj=="model"){return(tmp.model)} else {	return(OSC.results)	}
}
### plot OSC results ###modified by Kai

#fit many OPLS models to overview optimal LV and OLV
optimize.OPLS<-function(max.LV=4,tolerance =0.01,pls.y,pls.data,validation = "LOO",method="oscorespls",cv.scale=T,...){
    
    #iterate and fit OSC models for each possible LV > 1
    out<-lapply(1:max.LV, function(i){
        mod<-OSC.correction(pls.y=pls.y,pls.data=pls.data,comp=i,OSC.comp=i,validation = validation,cv.scale=cv.scale,...)
        tmp<-data.frame(RMSEP=do.call("rbind",mod$RMSEP))
        tmp$LV<-i
        tmp$OLV<-rep(mod$OSC.LVs,each=(ncol(pls.y)*(i+1)))
        tmp$pls.y<-rep(1:ncol(pls.y), each=(i+1))
        #do not report partials
        get<-matrix(rep(0:i),nrow=nrow(tmp))
        tmp[get==max(get),]
    })
    obj<-do.call("rbind",out)
    
    #choose optimal combination of LV/OLV for all Ys
    choose.opt.OPLS.comp(obj=obj,pls.y=pls.y,tolerance=0.01)
}

#choose optimal model LV and OLV component number
choose.opt.OPLS.comp<-function(obj,pls.y,tolerance=0.01){
    
    
    tmp.list<-split(obj,obj$pls.y)
    
    results<-lapply(1:length(tmp.list), function(i){
        x<-tmp.list[[i]]
        RMSEP<-x[,1:2]
        comp<-x$LV
        ocomp<-x$OLV
        
        even<-1:ncol(RMSEP)%%2==0 # CV RMSEP currently assuming this was used in modeling
        tmp<-RMSEP[,even]# CV RMSEP
        is.min<-which.min(tmp)
        min.RMSEP<-tmp[is.min]
        #look for smaller model with in tolerance
        # not worse than this, accept smaller
        delta<-tmp-min.RMSEP
        tmp.min<-which(delta<=tolerance)
        data.frame(x[c(tmp.min),], delta[tmp.min])
        
    })
    
    #choose smallest model within tolerance for both
    tmp<-do.call("rbind",results)
    x<-split(tmp$LV, tmp$pls.y )
    LV<-unlist(Reduce(intersect, x))
    x<-split(tmp$OLV, tmp$pls.y )
    OLV<-unlist(Reduce(intersect, x))
    
    #if there is an intersection
    if(length(LV)>0&length(OLV)>0){
        list(best=tmp[tmp$LV==min(LV)&tmp$OLV==min(OLV), ], LV=min(LV), OLV=min(OLV))
    } else {
        list(best=tmp, LV=tmp$LV[which.min(tmp$delta.tmp.min.)], OLV = tmp$OLV[which.min(tmp$delta.tmp.min.)])
    }
}
OSC.correction<-function(pls.y,pls.data,comp=5,OSC.comp=4,validation = "LOO",...){ # later open to  all plsr options, ...
    
    require(pls)
    
    #initialize
    OSC.results<-list()
    OSC.results$data[[1]]<-pls.data
    OSC.results$y[[1]]<-pls.y # also add a place to store plsr options for record keeping
    #  add a place to store plsr options for record keeping
    
    #need to iteratively fit models for each OSC
    for(i in 1:(OSC.comp+1)){
        data<-OSC.results$data[[i]]
        tmp.model<-plsr(OSC.results$y[[1]]~., data = data, ncomp = comp, validation = validation)#,...
        ww<-tmp.model$loading.weights[,1]
        pp<-tmp.model$loadings[,1]
        w.ortho<- pp - crossprod(ww,pp)/crossprod(ww)*ww
        t.ortho<- as.matrix(pls.data) %*% w.ortho
        p.ortho<- crossprod(as.matrix(data),t.ortho)/ c(crossprod(t.ortho))
        Xcorr<- data - tcrossprod(t.ortho,p.ortho)
        
        #store results
        OSC.results$RMSEP[[i]]<-matrix(RMSEP(tmp.model)$val,ncol=2,byrow=TRUE)
        OSC.results$scores[[i]]<-tmp.model$scores
        OSC.results$loadings[[i]]<-tmp.model$loadings
        OSC.results$loading.weights[[i]]<-tmp.model$loading.weights
        OSC.results$total.LVs[[i]]<-comp
        OSC.results$OSC.LVs[[i]]<-i-1 # account for first model not having any OSC LVs
        #initialize data for next round
        OSC.results$data[[i+1]]<-as.data.frame(Xcorr)
    }
    
    return(OSC.results)
}

### plot OSC results ###modified by Kai
plot.OSC.results<-function(obj,plot="RMSEP",groups=NULL){
    check.get.packages("ggplot2")
    # obj is from make.OSC.PLS.model
    #plot = one of: c("RMSEP","scores","loadings","delta.weights")
    #groups is a factor to show group visualization in scores plot
    switch(plot,
    RMSEP             =  .local<-function(obj){
        #bind info and RMSEP
        comps<-obj$total.LVs
        ocomps<-obj$OSC.LVs
        plot.obj<-obj$RMSEP
        bound<-do.call("rbind",lapply(1:length(comps),function(i)
        {
            out<-as.data.frame(cbind(plot.obj[[i]][,1],c(0:comps[i]),paste(comps[i]," LVs and ",ocomps[i]," OSC LVs",sep="")))
            colnames(out)<-c("RMSEP","component","model")
            out
        }))
        bound[,1:2]<-as.numeric(as.matrix(bound[,1:2]))
        
        #custom theme
        .theme<- theme(
        axis.line = element_line(colour = 'gray', size = .75),
        panel.background = element_blank(),
        plot.background = element_blank()
        )
        #plot
        p<-ggplot(data=bound, aes(x=component, y=RMSEP,color=model)) + geom_line(size=1,alpha=.5) + geom_point(size=2)+.theme
        print(p)
    },
    scores             =    .local<-function(obj){
        comps<-obj$total.LVs
        ocomps<-obj$OSC.LVs
        plot.obj<-obj$scores
        if(is.null(groups)){groups<-rep("gray",nrow(plot.obj[[1]][,]))}
        bound<-do.call("rbind",lapply(1:length(comps),function(i)
        {
            out<-as.data.frame(cbind(plot.obj[[i]][,1:2],unlist(groups),paste(comps[i]," LVs and ",ocomps[i]," OSC LVs",sep="")))
            colnames(out)<-c("Comp1","Comp2","groups","model")
            out
        }))
        bound[,1:2]<-as.numeric(as.matrix(bound[,1:2]))
        
        #calculate convex hull for polygons for each group
        data.obj <- split(bound, bound$model)
        tmp.obj <- lapply(1:length(data.obj), function(i){
            obj<-data.obj[[i]]
            s2<-split(obj,obj[,3])
            do.call(rbind,lapply(1:length(s2),function(j){
                tmp<-s2[[j]]
                tmp[chull(tmp[,1:2]),]
            }))
        })
        chull.boundaries <- do.call("rbind", tmp.obj)
        
        #custom theme
        .theme<- theme(
        axis.line = element_line(colour = 'gray', size = .75),
        panel.background = element_blank(),
        panel.border = element_rect(colour="gray",fill=NA),
        plot.background = element_blank()
        )
        
        #make plot
        p<-ggplot(data=bound, aes(x=Comp1, y=Comp2, group=groups,color=groups)) + #geom_density2d(aes(group=groups))+
        geom_hline(aes(yintercept=0),color="gray60",linetype="dashed")+geom_vline(aes(xintercept=0),color=I("gray60"),linetype=2)+facet_grid(. ~ model)
        p<-p+geom_polygon(data=chull.boundaries,aes(x=Comp1,y=Comp2,fill=groups),alpha=.5) +geom_point(size=2)+.theme
        print(p)
    },
    loadings         =     .local<-function(obj){ # will only plot first component for each model
        comps<-obj$total.LVs
        ocomps<-obj$OSC.LVs
        plot.obj<-obj$loadings
        bound<-do.call("rbind",lapply(1:length(comps),function(i)
        {
            out<-as.data.frame(cbind(plot.obj[[i]][,1:2],rownames(plot.obj[[i]]),paste(comps[i]," LVs and ",ocomps[i]," OSC LVs",sep="")))
            colnames(out)<-c("Comp1","Comp2","variable","model")
            out
        }))
        bound[,1:2]<-as.numeric(as.matrix(bound[,1:2]))
        
        #custom theme
        .theme<- theme(
        axis.line = element_line(colour = 'gray', size = .75),
        panel.background = element_blank(),
        legend.position = "none",
        plot.background = element_blank()
        )
        
        #make plot
        p<-ggplot(data=bound, aes(x=variable,y=Comp1, fill=variable)) + geom_bar(stat = "identity") + coord_flip() + #geom_density2d(aes(group=groups))+
        facet_grid(. ~ model) +.theme
        print(p)
    },
    delta.weights     =     .local<-function(obj){ # will only plot first component for each model
        comps<-obj$total.LVs
        ocomps<-obj$OSC.LVs
        plot.obj<-obj$loading.weights
        bound<-do.call("rbind",lapply(2:(length(ocomps)),function(i)
        {
            out<-as.data.frame(cbind(plot.obj[[1]][,1]-plot.obj[[i]][,1],names(plot.obj[[i]][,1]),paste(comps[i]," LVs and ",ocomps[i]," OSC LVs",sep="")))
            colnames(out)<-c("delta_weight","variable","model")
            out
        }))
        bound[,1]<-signif(as.numeric(as.matrix(bound[,1])),3)
        
        #theme
        .theme<- theme(
        axis.line = element_line(colour = 'gray', size = .75),
        panel.background = element_blank(),
        legend.position = "none",
        plot.background = element_blank()
        )
        #make plot
        p<-ggplot(data=bound, aes(x=variable,y=delta_weight, fill=variable)) + geom_bar(stat = "identity") + coord_flip() + #geom_density2d(aes(group=groups))+
        facet_grid(. ~ model) +.theme
        print(p)
    }
    )
    .local(obj)
}


#recreating plots based on plot.PCA options with slight modifications (good example of a place to use oob, need to have helper function to switch top level inputs based on class and use generic plotter)
#modified by Kai

plot.PLS<-function(obj, results = c("screeplot","scores","loadings","biplot"),xaxis=1,yaxis=2,size=3,color=NULL, shape=NULL, label=TRUE, legend.name =  NULL, font.size=5,group.bounds="ellipse",alpha=.5,g.alpha=.2,print.plot=TRUE,extra=NULL,repel=T,...){
    require(ggplot2)
    require(grid)
    require(ggrepel)
    require(reshape2)
    #obj is the results of type get.OSC.model
    #plot = one of: c("screeplot","scores","loadings","biplot","multi")
    #color is a factor to show group visualization of scores based on color
    
    
    local<-switch(results, # only tested for single Y models!
    RMSEP             =  function(obj,...){
        
        .theme<- theme(
        axis.line = element_line(colour = 'gray', size = .75),
        panel.background = element_blank(),
        plot.background = element_blank()
        )
        # RMSEP<-obj$RMSEP[length(obj$RMSEP)]][,ncol(obj$RMSEP[[1]])] # get for all LVs and optionally CV version
        # Q2<-obj$Q2[[length(obj$Q2)]][,ncol(obj$Q2[[1]])]
        # Xvar<-c(0,obj$Xvar[[length(obj$Xvar)]]) # 0 is for intercept only model
        
        RMSEP<-obj$RMSEP[,ncol(obj$RMSEP)] # get for all LVs and optionally CV version
        Q2<-obj$Q2[,ncol(obj$Q2)]
        Xvar<-cumsum(c(0,obj$Xvar)) #
        
        
        LV<-paste0("",0:(length(RMSEP)-1))
        tmp<-melt(data.frame(LV,RMSEP,Q2,Xvar))
        
        
        # RMSEP<-obj$RMSEP[,ncol(obj$RMSEP)] # get for all LVs and optionally CV version
        # Q2<-obj$Q2[,ncol(obj$Q2)]
        # Xvar<-c(0,obj$Xvar) #
        
        # LV<-as.character(c(0:(length(RMSEP)-1))) # account for intercept only model
        # tmp<-melt(data.frame(LV,RMSEP,Q2,Xvar),id=LV)
        
        p<-ggplot(data=tmp ,aes(y=value,x=LV,fill=variable))+
        geom_bar(stat="identity",position=position_dodge())+.theme +ylab("value")+xlab("LV") +extra
        if(print.plot){
            print(p)
        }     else {
            return(p)
        }
        
    },
    scores             =    function(obj,color,size,alpha,shape,...){
        comps<-obj$total.LVs[1]
        tmp.obj<-tryCatch(obj$scores[[comps]][,c(xaxis,yaxis)],error=function(e){obj$scores[,c(xaxis,yaxis)]}) # not sure how to simply unclass and coerce to data.frame
        
        tmp<-data.frame(tmp.obj,id = rownames(tmp.obj))
        #plot
        .theme2<- theme(
        axis.line = element_line(colour = 'gray', size = .75),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.background=element_rect(fill='white'),
        legend.key = element_blank()
        )
        
        if(is.null(color)){
            tmp$color<-"gray"
        }else{
            tmp$color<-as.factor(color[,])
            if(is.null(legend.name)){legend.name<-colnames(color)}
        }
        
        if(is.null(shape)){
            tmp$shape<-21
        }else{
            tmp$shape<-as.factor(shape[,])
        }
        
        points<-if(all(tmp$color=="gray")) {
            geom_point(color="gray",size=size,alpha=alpha,show_guide = FALSE)
        } else {
            geom_point(aes(color=color,),size=size,alpha=alpha)  # shape disabled
        }
        #labels
        tmp$lab.offset<-tmp[,2]-abs(range(tmp.obj[,2])[1]-range(tmp.obj[,2])[2])/50
        labels<-if(label==TRUE){
            geom_text_repel(size=font.size,aes_string(x=colnames(tmp)[1], y="lab.offset",label="id"),color="black",show_guide = FALSE)
            
        } else {
            NULL
            
        }
        
        #group visualizations
        #Hoettellings T2 ellipse
        polygons<-NULL
        if(group.bounds=="ellipse"){
            ell<-get.ellipse.coords(cbind(tmp.obj[,1],tmp.obj[,2]), group=tmp$color)# group visualization via
            polygons<-if(is.null(color)){
                geom_polygon(data=data.frame(ell$coords),aes(x=x,y=y), fill="gray", color="gray",linetype=2,alpha=g.alpha, show_guide = FALSE)
            } else {
                geom_polygon(data=data.frame(ell$coords),aes(x=x,y=y, fill=group),linetype=2,alpha=g.alpha, show_guide = FALSE)
            }
        }
        
        if(group.bounds=="polygon"){
            ell<-get.polygon.coords(data.frame(tmp.obj),tmp$color)# group visualization via
            polygons<-if(is.null(color)){
                geom_polygon(data=data.frame(ell),aes(x=x,y=y), fill="gray", color="gray",linetype=2,alpha=g.alpha, show_guide = FALSE)
            } else {
                geom_polygon(data=data.frame(ell),aes(x=x,y=y, fill=group),linetype=2,alpha=g.alpha, show_guide = FALSE)
            }
        }
        #making the actual plot
        p<-ggplot(data=tmp,aes_string(x=colnames(tmp)[1], y=colnames(tmp)[2])) +
        geom_vline(xintercept = 0,linetype=2, size=.5, alpha=.5) +
        geom_hline(yintercept = 0,linetype=2, size=.5, alpha=.5) +
        points +
        .theme2 +
        labels +
        polygons +
        scale_x_continuous(paste(colnames(tmp)[1],sprintf("(%s%%)", round(obj$Xvar[xaxis],digits=2)*100),sep=" "))+
        scale_y_continuous(paste(colnames(tmp)[2],sprintf("(%s%%)", round(obj$Xvar[yaxis],digits=2)*100),sep=" "))
        if(!is.null(legend.name)) {p<-p+scale_colour_discrete(name = legend.name)}
        p<-p+extra
        if(print.plot){
            print(p)
        }     else {
            return(p)
        }
    },
    "loadings"        = function(obj,color,size,alpha,...){
        comps<-obj$total.LVs[1]
        tmp.obj<-tryCatch(obj$loadings[[comps]][,c(xaxis,yaxis)],error=function(e){obj$loadings[,c(xaxis,yaxis)]}) # not sure how to simply unclass and coerce to data.frame
        tmp<-data.frame(tmp.obj,id = rownames(tmp.obj))
        #plot
        .theme2<- theme(
        axis.line = element_line(colour = 'gray', size = .75),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.background=element_rect(fill='white'),
        legend.key = element_blank()
        )
        #check to make sure color length matches dim[1]
        if(is.null(color)){
            tmp$color<-"gray"
        }else{
            if(!length(color[,])==nrow(tmp)){tmp$color<-"gray"# reset if doesn't match
            } else {
                tmp$color<-as.factor(color[,])
                if(is.null(legend.name)){legend.name<-colnames(color)}
            }
        }
        
        points<-if(all(tmp$color=="gray")) {
            geom_point(color="gray",size=size,alpha=alpha,show_guide = FALSE)
        } else {
            geom_point(aes(color=color),size=size,alpha=alpha)
        }
        #labels
        tmp$lab.offset<-tmp[,2]-abs(range(tmp.obj[,2])[1]-range(tmp.obj[,2])[2])/50
        labels<-if(label==TRUE){geom_text_repel(size=font.size,aes_string(x=colnames(tmp)[1], y="lab.offset",label="id"),color="black",show_guide = FALSE)} else { NULL }
        
        #group visualizations
        #Hoettellings T2 ellipse
        polygons<-NULL
        if(group.bounds=="ellipse"){
            ell<-get.ellipse.coords(cbind(tmp.obj[,1],tmp.obj[,2]), group=tmp$color)# group visualization via
            polygons<-if(all(tmp$color=="gray")){
                geom_polygon(data=data.frame(ell$coords),aes(x=x,y=y), fill="gray", color="gray",linetype=2,alpha=g.alpha, show_guide = FALSE)
            } else {
                geom_polygon(data=data.frame(ell$coords),aes(x=x,y=y, fill=group),linetype=2,alpha=g.alpha, show_guide = FALSE)
            }
        }
        
        if(group.bounds=="polygon"){
            ell<-get.polygon.coords(data.frame(tmp.obj),tmp$color)# group visualization via
            polygons<-if(all(tmp$color=="gray")){
                geom_polygon(data=data.frame(ell),aes(x=x,y=y), fill="gray", color="gray",linetype=2,alpha=g.alpha, show_guide = FALSE)
            } else {
                geom_polygon(data=data.frame(ell),aes(x=x,y=y, fill=group),linetype=2,alpha=g.alpha, show_guide = FALSE)
            }
        }
        
        #making the actual plot
        p<-ggplot(data=tmp,aes_string(x=colnames(tmp)[1], y=colnames(tmp)[2])) +
        geom_vline(xintercept = 0,linetype=2, size=.5, alpha=.5) +
        geom_hline(yintercept = 0,linetype=2, size=.5, alpha=.5) +
        points +
        .theme2 +
        labels +
        polygons +
        scale_x_continuous(paste(colnames(tmp)[1],sprintf("(%s%%)", round(obj$Xvar[xaxis],digits=2)*100),sep=" "))+
        scale_y_continuous(paste(colnames(tmp)[2],sprintf("(%s%%)", round(obj$Xvar[yaxis],digits=2)*100),sep=" "))
        if(!is.null(legend.name)) {p<-p+scale_colour_discrete(name = legend.name)}
        p<-p+extra
        if(print.plot){
            print(p)
        }     else {
            return(p)
        }
    },
    "biplot"        = function(obj,color,size,alpha,...){
        require(scales)
        comps<-obj$total.LVs[1]
        loadings<-tmp.loadings<-tryCatch(obj$loadings[[comps]][,c(xaxis,yaxis)],error=function(e){obj$loadings[,c(xaxis,yaxis)]}) # not sure how to simply unclass and coerce
        scores<-tmp.obj<-data.frame(tryCatch(obj$scores[[comps]][,c(xaxis,yaxis)],error=function(e){obj$scores[,c(xaxis,yaxis)]})) # not sure how to simply unclass and coerce to data.frame
        .theme2<- theme(
        axis.line = element_line(colour = 'gray', size = .75),
        panel.background = element_blank(),
        plot.background = element_blank()
        )
        #based on https://groups.google.com/forum/#!topic/ggplot2/X-o2VXjDkQ8
        tmp.loadings[,1]<-rescale(loadings[,1], range(scores[,1]))
        tmp.loadings[,2]<-rescale(loadings[,2], range(scores[,2]))
        tmp.loadings<-data.frame(tmp.loadings,label=rownames(loadings))
        
        #using tmp.obj because no need for labels, and started badly need to rewrite
        #Adding Hoettellings T2 ellipse
        if(is.null(color)){
            tmp.obj$color<-"gray"
        }else{
            tmp.obj$color<-as.factor(color[,])
            if(is.null(legend.name)){legend.name<-colnames(color)}
        }
        
        #group visualizations
        #Hoettellings T2 ellipse
        polygons<-NULL
        if(group.bounds=="ellipse"){
            ell<-get.ellipse.coords(cbind(tmp.obj[,1],tmp.obj[,2]), group=tmp.obj$color)# group visualization via
            polygons<-if(all(tmp.obj$color=="gray")){
                geom_polygon(data=data.frame(ell$coords),aes(x=x,y=y), fill="gray", color="gray",linetype=2,alpha=g.alpha, show_guide = FALSE)
            } else {
                geom_polygon(data=data.frame(ell$coords),aes(x=x,y=y, fill=group),linetype=2,alpha=g.alpha, show_guide = FALSE)
            }
        }
        
        if(group.bounds=="polygon"){
            ell<-get.polygon.coords(data.frame(tmp.obj[,1:2]),tmp.obj$color)# group visualization via
            polygons<-if(all(tmp.obj$color=="gray")){
                geom_polygon(data=data.frame(ell),aes(x=x,y=y), fill="gray", color="gray",linetype=2,alpha=g.alpha, show_guide = FALSE)
            } else {
                geom_polygon(data=data.frame(ell),aes(x=x,y=y, fill=group),linetype=2,alpha=g.alpha, show_guide = FALSE)
            }
        }
        
        points<-if(all(tmp.obj$color=="gray")) {
            geom_point(data=data.frame(tmp.obj),aes_string(x=colnames(tmp.obj)[1], y=colnames(tmp.obj)[2]),color="gray",size=size,alpha=alpha,show_guide = FALSE)
        } else {
            geom_point(data=data.frame(tmp.obj), aes_string(x=colnames(tmp.obj)[1], y=colnames(tmp.obj)[2],color="color"),size=size,alpha=alpha)
        }
        #plot
        p<-ggplot()+
        points +
        polygons+
        geom_segment(data=tmp.loadings, aes_string(x=0, y=0, xend=colnames(tmp.loadings)[1], yend=colnames(tmp.loadings)[2]), arrow=NULL, alpha=0.25)+
        geom_text_repel(data=tmp.loadings, aes_string(x=colnames(tmp.loadings)[1], y=colnames(tmp.loadings)[2], label="label"), alpha=0.5, size=font.size)+
        scale_colour_discrete("Variety")+
        scale_x_continuous(paste(colnames(tmp.obj)[1],sprintf("(%s%%)", round(obj$Xvar[xaxis],digits=2)*100),sep=" "))+
        scale_y_continuous(paste(colnames(tmp.obj)[2],sprintf("(%s%%)", round(obj$Xvar[yaxis],digits=2)*100),sep=" ")) +
        .theme2
        if(!is.null(legend.name)) {p<-p+scale_colour_discrete(name = legend.name)}
        p<-p+extra
        if(print.plot){
            print(p)
        }     else {
            return(p)
        }
    }
    )
    
    local(obj,color=color,size=size,alpha=alpha,shape,...)
}

#create PLS model
make.PLS.model<-function(y,data,pls.method="simpls",
ncomp=2, CV="LOO",CV.segments=NULL,segment.type=NULL, CV.scale=FALSE, opt.comp=FALSE){
    #use opt.comp=TRUE to dynamically optimize the # of latent variables
    #minimum of 2 components
    #need to switch based on CV specifications
    
    #make sure the number of supplied components won't error plsr
    if(ncomp>nrow(data)-1)
    {
        ncomp<-nrow(data)-1
        #cat("The number of components was changed to", ncomp, "to accommodate sample number","\n")
    }
    if(CV=="LOO")
    
    {    if(opt.comp==TRUE)
        {
            mod1 <- plsr(y~ as.matrix(data), ncomp=ncomp,
            data=as.data.frame(data) ,method=pls.method,validation=CV,scale=CV.scale)
            new.comp<-c(2:ncomp)[which.max(R2(mod1)$val[-c(1:2)])]
            #cat("PCs were changed from",PCs,"to", new.comp)
            if(dim(as.data.frame(new.comp))[1]==0){new.comp<-2}
            mod1 <- plsr(y~ as.matrix(data), ncomp=new.comp,
            data=as.data.frame(data) ,method=pls.method,validation=CV,scale=CV.scale)
        }else{
            mod1 <- plsr(y~ as.matrix(data), ncomp=ncomp,
            data=as.data.frame(data) ,method=pls.method,validation=CV,scale=CV.scale)
        }
    }else{
        if(opt.comp==TRUE)
        {
            mod1 <- plsr(y~ as.matrix(data), ncomp=ncomp,
            data=as.data.frame(data) ,method=pls.method,validation=CV,segments=CV.segments,segment.type=segment.type,scale=CV.scale)
            new.comp<-c(2:ncomp)[which.max(R2(mod1)$val[-c(1:2)])]
            if(dim(as.data.frame(new.comp))[1]==0){new.comp<-2}
            mod1 <- plsr(y~ as.matrix(data), ncomp=ncomp,
            data=as.data.frame(data) ,method=pls.method,validation=CV,segments=CV.segments,segment.type=segment.type,scale=CV.scale)
        }else{
            mod1 <- plsr(y~ as.matrix(data), ncomp=ncomp,
            data=as.data.frame(data) ,method=pls.method,validation=CV,segments=CV.segments,segment.type=segment.type,scale=CV.scale)
        }
    }
    
    mod1
}

#extract OSC submodel from OSC results object
## modified by Kai
get.OSC.model<-function(obj,OSC.comp){
    #obj = results from OSC.correction()
    #OSC.comp = number of orthogonally corrected components
    
    index<-c(1:length(obj$OSC.LVs))
    id<-index[obj$OSC.LVs==OSC.comp]
    
    #extract and return
    out<-list()
    out$data<-obj$data[[id]]
    out$y<-obj$y
    out$fitted.values<-data.frame(obj$fitted.values[[id]][,,dim(obj$fitted.values[[id]])[3],drop=FALSE])
    colnames(out$fitted.values)<-paste0("Y_",c(1:ncol(out$y[[1]])),"_fitted.values")
    out$residuals<-data.frame(out$fitted.values-out$y[[1]])
    colnames(out$residuals)<-paste0("Y_",c(1:ncol(out$y[[1]])),"_residuals")
    out$RMSEP<-obj$RMSEP[[id]]
    out$predicted.RMSEP<-obj$predicted.RMSEP[[id]]
    out$predicted.Y<-obj$predicted.Y[[id]]
    out$Q2<-obj$Q2[[id]]
    out$Xvar<-obj$Xvar[[id]]
    out$scores<-obj$scores[[id]]
    out$loadings<-obj$loadings[[id]]
    out$loading.weights<-obj$loading.weights[[id]]
    out$total.LVs<-obj$total.LVs[[id]]
    out$OSC.LVs<-obj$OSC.LVs[[id]]
    out$VIP<-obj$VIP
    out$OPLSDA.stats<-obj$OPLSDA.stats[[id]]
    out$coefficients<-data.frame(obj$coefficients[[id]])
    out$model<-data.frame("Xvar"=c(0,round(cumsum(out$Xvar)*100,2)),"Q2"=out$Q2,"RMSEP"= out$RMSEP)
    colnames(out$coefficients)<-"coefficients"
    return(out)
}

#calculate root mean squared error
RMSE<-function(values,predictions){sqrt(sum((values-predictions)^2)/length(values))}

#model function
model.fxn<-function(data,inds,algorithm="pcr",y,ncomp=2,return="pred.error",...){
    mod<-do.call(algorithm,list(formula=y~.,data=data,ncomp=ncomp,subset=inds,...=...))
    switch(return,
    "pred.error" = .local<-function(){
        y-predict(mod, newdata=data,ncomp=ncomp,...)
    },
    "coef"         =     .local<-function(){
        c(coef(mod))
    }
    )
    .local()
}

#bootstrap function
boot.fxn<-function(algorithm="pcr",data=tmp.data,y,ncomp=2,return="pred.error", R=499,...){
    library(boot);library(pls)
    boot(data,statistic=model.fxn,R=R,algorithm=algorithm,ncomp=ncomp,y=y,return=return,...=...)
}

#generate boostrapped parameters for model (only works for single Y models)
boot.model<-function(algorithm="pcr",data=tmp.data,y,ncomp=2,return="pred.error",R=499,...){
    # function currently tailored to pls and tested with pcr and plsr
    # return can be one of c("pred.error","coef")
    # pred.error = out of bag (sample) error of prediction (RMSEP)
    # coef = bootstrapped coefficent weights
    
    x<-boot.fxn(algorithm,data,y,ncomp,return,R,...)
    in.bag<-boot.array(x)
    out.bag<-in.bag==0
    switch(return,
    "pred.error"     = .local<-function(){
        in.bag<-boot.array(x)
        oob.error<-mean((x$t^2)[in.bag==0])
        oob.error.sd<-sd((x$t^2)[in.bag==0])
        app.error<-MSEP(do.call(algorithm,list(formula=y~.,data=data,ncomp=ncomp,...=...)),ncomp=ncomp,intercept=FALSE)
        est.error<-sqrt(0.368*c(app.error$val) + 0.632 * oob.error)
        sd.error<-sqrt(0.368*c(app.error$val) + 0.632 * oob.error.sd)
        return(list(RMSEP=data.frame(bootstrapped.0.632_RMSEP = est.error,mean.error = mean(abs(x$t0)),sd.error=sd(abs(x$t0))),boot.obj = x))
    },
    "coef"             = .local<-function(){
        
        return(list(coef=data.frame(bootstrapped.coef=x$t0,CI.95.percent=t(apply(x$t,2,quantile, c(0.025,.975)))),boot.obj = x))
    }
    )
    .local()
}

#function to select upper boundaries of a vector based on quantile or number
#returns row id of positions (or logical)
feature.cut2<-function(obj,type="quantile",thresh=.9,separate=FALSE){
    # type can be one of c("number", "quantile")
    # thresh = select values above quantile or number of largest values
    # separate = value tested separately based on sign
    
    obj<-fixln(obj) # make a vector of type numeric
    #get position of selections
    tmp<-data.frame(id=1:length(obj),obj=obj)
    if(separate==TRUE){
        tmp.l<-split(abs(tmp),sign(tmp$obj))
    } else {
        tmp.l<-list(abs(tmp))
    }
    
    if(type=="quantile"){
        filter<-lapply(1:length(tmp.l),function(i){
            cut.point<-quantile(tmp.l[[i]][,2],prob=thresh)
            tmp.l[[i]][tmp.l[[i]][,2]>=cut.point,1]
        })
    }
    
    if(type=="number"){
        filter<-lapply(1:length(tmp.l),function(i){
            cut.point<-thresh
            x<-tmp.l[[i]][order(tmp.l[[i]][,2],decreasing=TRUE),]
            x[1:nrow(x)<=cut.point,1]
        })
    }
    
    #a logical vector
    # return(c(1:length(obj))%in%unlist(filter))
    
    #return row
    return(unlist(filter))
}

#function to bootstrap many models
multi.boot.model<-function(algorithm="pcr",data=data,y,
feature.subset=NULL,ncomp=4,return="pred.error",R=10,
parallel=FALSE,plot=TRUE,...){
    
    .local<-function(var.id,...){
        boot.model(algorithm=algorithm,data=data[,var.id],y=y,ncomp=ncomp,return=return,R=R,...)[[1]] # no boot obj returned
    }
    
    if (parallel == TRUE){
        check.get.packages(c("snow","doSNOW","foreach"))
        #start cluster
        cl.tmp = makeCluster(rep("localhost",Sys.getenv('NUMBER_OF_PROCESSORS')), type="SOCK")
        registerDoSNOW(cl.tmp)
        
        #
        out<-list()
        out<-foreach(i=c(1:length(feature.subset)),.combine="rbind") %dopar% .local(var.id=feature.subset[[i]])
        stopCluster(cl.tmp)
    } else {
        pb <- txtProgressBar(min = 0, max = length(feature.subset), style = 3)
        out<-do.call("rbind",lapply(1:length(feature.subset),
        function(i){
            setTxtProgressBar(pb, i)
            .local(var.id=feature.subset[[i]])})) #,...
        dimnames(out)<-list(c(1:length(feature.subset)),c("bootstrapped.0.632_RMSEP", "mean.error", "sd.error"))
        close(pb)
    }
    
    #output
    val<-data.frame(feature.set=c(1:nrow(out)),number.of.features=sapply(feature.subset,length),RMSEP_0.632=as.numeric(out[,1]),
    mean.error = out[,2] )
    
    #calculate loess mins
    get.lo<-function(){
        lo <- loess(RMSEP_0.632 ~ number.of.features, data=val)
        which.min(predict(lo,new.data=val))}
    lo.min<-tryCatch(get.lo(),error=function(e){NULL})
    
    out<-list()
    out$results<-val
    out$loess.min<-lo.min
    out$RMSEP_0.632.min<-val$number.of.features[which.min(val$RMSEP_0.632)]
    
    if(plot == TRUE){
        check.get.packages(c("ggplot2","reshape"))
        
        plot.obj<-data.frame(melt(val[,-c(1:2)]),features=rep(val$number.of.features,(ncol(val)-2)))
        print(plot.obj)
        #make plot
        p<-ggplot(data=plot.obj, aes(y=value, x=features,group=variable,color=variable, fill=variable)) +
        xlab("number of features")
        
        if(length(lo.min)>0){
            p<-p+stat_smooth(level = 0.95,size=.75,alpha=.15,legend=FALSE) +
            geom_vline(xintercept = out$RMSEP_0.632.min,lty=2,col="red") +
            ggtitle(paste("minimum at ",out$RMSEP_0.632.min," features"))+
            geom_point(size=3,alpha=.75)
        } else {
            p<-p+geom_point(size=3,alpha=.75)
        }
        print(p)
    }
    return(out)
}

get.ellipse.coords<-function(obj,group=NULL, ellipse.level=.95){
    check.get.packages(c("ellipse","splancs"))
    
    fct<-if(is.null(group)) as.factor(rep(1,nrow(obj))) else factor(unlist(group))
    .obj<-split(as.data.frame(obj),fct)
    names<-c(colnames(obj),"")
    #calculate points for ellipse
    #for level of group
    ellipse.var<-lapply(1:nlevels(fct),function(i)
    {
        pts<-.obj[[i]]
        m<-colMeans(pts)
        tmp<-cbind(tryCatch(ellipse::ellipse(as.matrix(cov(pts)),centre=c(m[1],m[2]),level=ellipse.level),
        error=function(e){matrix( c(NA,NA),nrow=1)}),rep(levels(fct)[i],nrow(pts)))
        colnames(tmp)<-NULL
        tmp
    })
    
    #format for ggplot 2
    tmp<-do.call("rbind",ellipse.var)
    #remove errors
    tmp<-tmp[!is.na(tmp[,1]),]
    
    ellipse.size<-sapply(1:length(ellipse.var),function(i)
    {
        tryCatch(areapl(ellipse.var[[i]]),error=function(e){NA})
    })
    #avoiding issues with x/y as factors
    obj<-data.frame(matrix(as.numeric(as.matrix(tmp[,1:2])),ncol=2),tmp[,3])
    colnames(obj)<-c("x","y","group")
    #may need to maintain ordered factor levels
    obj$group<-factor(obj$group,levels=levels(fct),ordered=is.ordered(fct))
    return(list(coords=data.frame(obj), area=ellipse.size))
}

#get polygon coordinates for each group
get.polygon.coords<-function(obj,group){
    require(plyr)
    fct<-if(is.null(group)) data.frame(fct=as.factor(rep(1,nrow(obj)))) else factor(unlist(group))
    obj<-data.frame(obj,fct)
    .chull<-function(x){tryCatch(chull(x),error=function(e){NA})} #
    chull.bounds <- data.frame(ddply(obj, .(fct), function(x) data.frame(x[.chull(as.matrix(x)),])))
    colnames(chull.bounds)<-c("x","y","group")
    #may need to maintain ordered factor levels
    chull.bounds$group<-factor(chull.bounds$group,levels=levels(fct),ordered=is.ordered(fct))
    return(chull.bounds)
}
### modified by Kai
feature.bar.plot<-function(feature.set,weights.set,extra.plot=NULL){
    #feature.set = index for selected features
    #all weights to plot as a bar graph
    #show.all determines if only selected or all features are ploted
    
    #Bargraphs
    
    #plotting data
    bound<-data.frame(weights=weights.set,
    variable=c(1:length(weights.set)),
    show = c(1:length(weights.set))%in%feature.set)
    
    #cut offs
    cuts<-range(bound$weights[!bound$show])
    plot.title<- paste ("upper/lower bounds = ", signif(cuts[2],4), " / " ,signif(cuts[1],4))
    #plot.colors<-scale_fill_brewer(palette="Blues")
    
    #make plot of variable and weight
    p<-ggplot(data=bound, aes(x=variable,y=weights, fill=show)) +
    geom_bar(stat = "identity") + #geom_density2d(aes(group=groups))+
    theme_minimal(base_size = 15) +geom_hline(yintercept = cuts,lty=2,col="red") +
    labs(title = plot.title, fill= "Selected") #+
    #plot.colors
    
    
    # sorted weight
    sorted.bound<-bound[order(bound$weights),]
    sorted.bound$variable<-c(1:length(bound$weights))
    
    
    p2<-ggplot(data=sorted.bound, aes(x=variable,y=weights, fill=show)) +
    geom_bar(stat = "identity") + xlab(" ") + #geom_density2d(aes(group=groups))+
    theme_minimal(base_size = 20) + geom_hline(yintercept = cuts,lty=2,col="red")# +
    #plot.colors
    #print(p2)
    
    #print plots
    if(is.null(extra.plot)){
        multiplot (p,p2, plotlist=NULL, cols=1)
    } else {
        multiplot (extra.plot,p,p2, plotlist=NULL, cols=1)
    }
}

###make the S plot ###modified by Kai
make.S.plot<-function(pls.data,pls.scores,pls.loadings, cut.off=0.05, FDR=TRUE,plot=TRUE,...){
    
    check.get.packages("ggplot2")
    
    #pls.data     = data used for model
    #scores     = scores for selected (1st) component
    #loadings     = loadings for selected (1st) component
    #plot         = make S-Plot
    #cut.off     = select optimal features based on a test of the significance of the correlation (pearsons) between variable and scores
    #FDR         = use q-value as the cut off
    #...         = can specify correlation type = c("pearson","spearman","biweight")
    
    # calculate p(corr) or correlation between scores and the original variable
    cor.mat<-calculate.correlations(cbind(pls.scores,pls.data),results="matrix",...) #
    corrs<-cor.mat$cor[-1,1]
    p.vals<-cor.mat$p.value[-1,1]
    
    #false discovery rate correction
    if(FDR==TRUE){
        #p.vals<-FDR.adjust(p.vals,type="pvalue",return.all=TRUE)$qval #
        p.vals<-p.adjust(p.vals, method="BH")
    }
    
    #index to draw visualization
    show<-p.vals
    show[]<-1
    if(is.numeric(cut.off)){
        show[p.vals>cut.off]<-0
    }
    
    #make plot
    plot.obj<-data.frame(pcorr=corrs,loadings=pls.loadings,value=p.vals, significant=as.logical(show))
    
    if(plot==TRUE){
        #theme
        # .theme<- theme(
        #    axis.line = element_line(colour = 'black', size = .75),
        #    panel.background = element_blank(),
        #     #legend.position = "none",
        #    plot.background = element_blank()
        #  )
        
        #cut offs
        selected<-plot.obj$significant==1
        plot.title<- paste (sum(selected)," selected features / ",round(sum(selected)/length(pls.loadings)*100,0),"%",sep="")
        
        
        #make plot of variable and weight
        p<-ggplot(data=plot.obj, aes(x=loadings,y=pcorr, color=significant)) +
        geom_point(stat = "identity",alpha=.75,size=2) + #geom_density2d(aes(group=groups))+
        labs(title = plot.title, fill= "Selected")+theme_minimal(base_size = 20)
        print(p)
    } else { p<-"NULL"}
    
    return(list(plot=p,selected.data=pls.data[,plot.obj$significant==1],feature.info=plot.obj))
}

##### modified by Kai
plot.S.plot<-function(obj,names=NULL,return=c("all","splot","barplot","top"),extra=NULL, plot=TRUE){
    
    check.get.packages(c("ggplot2","gridExtra"))
    
    #model.weight             = x cut
    #pcorr                     = correlation between scores and variables (y cut)
    #combined.selection     = selected features
    
    
    #make plot
    plot.obj<-data.frame(names=rownames(obj),pcorr=obj$pcorr,loadings=obj$model.weight, significant=obj$combined.selection)
    
    
    #cut offs
    selected<-plot.obj$significant==1
    plot.title<- paste (sum(selected)," selected features ( ",round(sum(selected)/length(plot.obj$loadings)*100,0),"%)",sep="")
    
    
    #make S plot of variable and weights
    p1<-ggplot(data=plot.obj, aes(x=loadings,y=pcorr, color=significant)) +
    geom_point(stat = "identity",alpha=.75,show_guide=FALSE,size=2) + #geom_density2d(aes(group=groups))+
    theme_minimal(base_size=15) + labs(title = plot.title, fill= "Selected") +extra
    
    #
    #feature.set = index for selected features
    #all weights to plot as a bar graph
    #show.all determines if only selected or all features are ploted
    
    #plotting data (redundant object!, need to stay consistent with names to avoid this idiocy)
    bound<-data.frame(name=plot.obj$name,weights=plot.obj$loadings,
    variable=c(1:length(plot.obj$loadings)),
    show = plot.obj$significant)
    
    #cut offs
    tmp<-bound$weights[bound$show]
    tmp<-split(tmp,sign(tmp))
    
    cuts<-c(max(tmp[['1']]),min(tmp[['-1']]))
    plot.title<- paste ("upper/lower bounds = ", signif(cuts[2],4), " / " ,signif(cuts[1],4))
    #plot.colors<-scale_fill_brewer(palette="Blues")
    
    #make plot of variable and weight
    p2<-ggplot(data=bound, aes(x=variable,y=weights, fill=show)) +
    geom_bar(stat = "identity") + #geom_density2d(aes(group=groups))+
    theme_minimal(base_size=15) +geom_hline(yintercept = cuts,lty=2,col="red") +
    labs(title = plot.title, fill= "Selected") +extra #+
    #plot.colors
    
    .theme<- theme(
    axis.line = element_line(colour = 'black'),
    panel.background = element_blank(),
    legend.position = "none",
    plot.background = element_blank(),
    axis.title.x=element_text(size=15),
    axis.text.y=element_text(size=5)
    )
    # sorted weight and show names in vertical bar plot (selected only)
    sorted.bound<-bound[order(bound$weights,decreasing=FALSE),]
    sorted.bound$index<-1:nrow(sorted.bound)
    #remove unused factor levels
    sorted.bound<-sorted.bound[sorted.bound$show==TRUE,,drop=FALSE]
    if(nrow(sorted.bound)>0){
        sorted.bound$name<-fixlc(sorted.bound$name)
        sorted.bound$index<-1:nrow(sorted.bound)
        p3<-ggplot(sorted.bound, aes(x = index, y = weights, fill = show))+
        geom_bar(stat = "identity",show_guide=FALSE,fill="gray") + xlab(" ") + #geom_density2d(aes(group=groups))+
        .theme+ geom_hline(yintercept = cuts,lty=2,col="red") + coord_flip() +
        scale_x_continuous(breaks=c(1:length(sorted.bound$name)),labels=as.character(sorted.bound$name)) + extra
    } else {p3<-ggplot(sorted.bound, aes(x = index, y = weights, fill = show))+geom_bar(stat = "identity",show_guide=FALSE,fill="gray")}
    
    
    #plot
    if(plot){
        switch(return,
        "all" = print(grid.arrange(p1,p2,p3, ncol = 1)),
        "splot" = print(p1),
        "barplot" = print(p2) ,
        "top" = print(p3)
        )
    } else {
        switch(return,
        "all" = list(tmp.list),
        "splot" = p1,
        "barplot" = p2 ,
        "top" = p3
        )
    }
    
}
####modified by Kai
#feature select using a combination of analyte correlation to scores (S-plot) and feature weights
PLS.feature.select<-function(pls.data,pls.scores,pls.loadings,pls.weight,plot=TRUE,p.value=0.05, FDR=TRUE,
cut.type="quantile",top=0.95,separate=TRUE,...){
    #combined args from
    #feature.cut() & make.S.plot()
    #cuts is a single value which is a propability for type = quantile or integer for number
    
    #first selection criteria based on magnitude of model weight
    # weight.cut<-feature.cut(obj=pls.weight,type=cut.type,cuts=top,separate=separate,plot=FALSE)
    weight.cut<-feature.cut2(obj=pls.weight,type=cut.type,thresh=top,separate=separate)
    weight.cut.selected<-data.frame(selected.weights=rep(0,length(pls.weight)))
    weight.cut.selected[unlist(weight.cut),]<-1
    
    #second selection criteria based on variable correlation with scores
    cor.cut<-make.S.plot(pls.data=pls.data,pls.scores=pls.scores,pls.loadings=unname(pls.loadings),cut.off=p.value, FDR=FDR,plot=FALSE,...)
    
    #combine and plot
    combo.cut<-data.frame(model.weight=unname(pls.weight),weight.cut.selected, cor.cut$feature.info)
    combo.cut$combined.selection<-combo.cut$significant&combo.cut$selected.weights==1
    
    #return results!!! check
    # invisible(as.data.frame(combo.cut))
    
    if(plot){
        #create updated S-plot
        plot.obj<-combo.cut
        selected<-plot.obj$combined.selection==1
        plot.title<- paste (sum(selected)," selected features or ",round(sum(selected)/length(pls.loadings)*100,0),"%",sep="")
        
        #theme
        #.theme<- theme(
        #  axis.line = element_line(colour = 'gray', size = .75),
        #  panel.background = element_blank(),
        #  legend.position = "none",
        #  plot.background = element_blank()
        #)
        
        
        #make plot of variable and weight
        p<-ggplot(data=plot.obj, aes(x=loadings,y=pcorr, color=combined.selection)) +
        geom_point(stat = "identity",alpha=.75) + #geom_density2d(aes(group=groups))+
        theme_minimal(base_size=15) + labs(title = plot.title, fill= "Selected")
        
        
        #plot results
        feature.bar.plot(feature.set=c(1:length(combo.cut$model.weight))[combo.cut$combined.selection],weights.set=combo.cut$model.weight, extra.plot=p)
    }
    
    #return results
    # return(as.data.frame(combo.cut))
    return(invisible(as.data.frame(combo.cut)))
}

#permute PLS model
permute.PLS<-function(data,y,n=10,ncomp,...){# could be made parallel
    #permuted Y
    perm.y<-lapply(1:n,function(i)
    {
        apply(y,2,gtools::permute)
    })
    
    #generate permuted models
    model<-lapply(1:n,function(i)
    {
        # cat("permuting model",i,"\n")
        model<-make.PLS.model(y=perm.y[[i]],data,ncomp=ncomp,...)
        #get stats
        q2<-R2(model)$val[,,ncomp+1]
        rx2<-drop(model$Xvar/model$Xtotvar)[ncomp]
        pred.val<-model$fitted.values[,,ncomp]
        rmsep<-pls::RMSEP(model)$val[2,,ncomp] # take CV adjuste RMSEP
        list(Q2=q2,RX2=rx2,RMSEP=rmsep)#,predicted=pred.val,actual=perm.y[[i]])
    })
    
    tmp<-matrix(unlist(do.call("rbind",model)),ncol=3)
    colnames(tmp)<-c("Q2","Xvar","RMSEP")
    means<-apply(tmp,2,mean)
    sds<-apply(tmp,2,sd)
    summary<-paste(signif(means,4)," ± ", signif(sds,3))
    return(list(permuted.values=tmp, mean = means, standard.deviations = sds, summary = summary))
}

#permute OSC-PLS model
permute.OSC.PLS<-function(data,y,n=10,ncomp,OSC.comp=1,train.test.index=NULL,...){ # should be made parallel
    
    #only using first Y
    # message("Only first Y permuted")
    #permuted Y
    perm.y<-lapply(1:n,function(i)
    {
        apply(y[,1,drop=FALSE],2,gtools::permute)
    })
    
    #collect correlation between y and permuted y
    # disabled for multi y
    if(ncol(y)==1){
        cor.with.y<-data.frame(correlation=abs(cor(cbind(y,do.call("cbind",perm.y))))[-1,1])
    } else {
        cor.with.y<-NULL
    }
    
    #generate permuted models
    model<-sapply(1:n,function(i)
    {
        # cat("permuting model",i,"\n")
        if(!is.null(train.test.index)) {tmp.train.test.index<-train.test.index[,i,drop=FALSE]} else {tmp.train.test.index<-train.test.index}
        #if train data contains all the same value then this will cause an error for OPLS
        model<-make.OSC.PLS.model(pls.y=as.matrix(perm.y[[i]]),pls.data=data,comp=ncomp,OSC.comp=OSC.comp,train.test.index=tmp.train.test.index,...) #,...
        tmp.OSC.comp<-max(model$OSC.LVs)
        #get stats
        q2<-model$Q2[[tmp.OSC.comp+1]][ncomp+1,,drop=FALSE]# cv adjusted
        rx2<-round(sum(model$Xvar[[tmp.OSC.comp+1]])*100,1)
        pred.val<-as.matrix(model$fitted.values[[tmp.OSC.comp+1]][,,ncomp])
        rmsep<-model$rmsep[[tmp.OSC.comp+1]]# take CV adjusted internal RMSEP (see true RMSEP below)
        if(!is.null(train.test.index)) {
            rmsep<-model$predicted.RMSEP[[tmp.OSC.comp+1]]
        }
        if(!is.null(model$OPLSDA.stats)){oplsda.stats<-data.frame(model$OPLSDA.stats[[tmp.OSC.comp+1]])} else {oplsda.stats<-data.frame(empty=NA)}
        
        data.frame(RX2=rep(rx2,length(q2)),Q2=q2,RMSEP=rmsep,oplsda.stats)#,predicted=pred.val,actual=perm.y[[i]])
    })
    #return results
    tmp<-as.matrix(t(model[!rownames(model)=="empty",]))
    names<-colnames(tmp)
    tmp<-matrix(unlist(tmp),ncol=ncol(tmp))#needs to be atomic has to be a better way
    colnames(tmp)<-names
    means<-apply(tmp,2,mean, na.rm=TRUE)
    sds<-apply(tmp,2,sd,na.rm=TRUE)
    summary<-matrix(paste(signif(means,4),"+/-", signif(sds,3)),ncol=length(sds))
    colnames(summary)<-colnames(tmp)
    return(list(permuted.values=cbind(tmp,cor.with.y), mean = means, standard.deviations = sds, summary = summary))
}

#IMPROVED version of permute.OSC.PLS, using a modification of the test/train OSC.PLS.train.test to return full prediction results
permute.OSC.PLS.train.test<-function(pls.data,pls.y,perm.n=10,train.test.index,comp,OSC.comp,...) {
    pls.y<-as.matrix(pls.y)
    #permute the Y
    perm.y<-lapply(1:perm.n,function(i)
    {
        apply(pls.y[,1,drop=FALSE],2,gtools::permute)
    })
    
    results<-lapply(1:ncol(train.test.index), function(i){
        pls.y<-as.matrix(perm.y[[i]])
        pls.train.index<-as.matrix(train.test.index[,i])
        #order for merging test with train stats in one object
        new.order<-c(c(1:nrow(pls.data))[pls.train.index=="train"],c(1:nrow(pls.data))[pls.train.index=="test"])
        back.sort<-order(new.order)
        
        train.y<-train.real<-pls.y[pls.train.index=="train",]
        train.data<-pls.data[pls.train.index=="train",]
        test.real<-pls.y[pls.train.index=="test",]
        #all arguments gave been preset elsewhere
        test.pls.results<-make.OSC.PLS.model(pls.y=pls.y,pls.data=pls.data,comp=comp,OSC.comp=OSC.comp, train.test.index=pls.train.index,...)
        tmp.OSC.comp<-max(test.pls.results$OSC.LVs) # control when limited with Orthogonal dimensions
        Q2<-data.frame(Q2=test.pls.results$Q2[[tmp.OSC.comp+1]][comp,])
        
        #fitted values
        train.pred<-test.pls.results$fitted.values[[tmp.OSC.comp+1]][,,comp]
        test.pred<-test.pls.results$predicted.Y[[tmp.OSC.comp+1]]
        
        RMSEP<-data.frame(RMSEP=test.pls.results$predicted.RMSEP[[tmp.OSC.comp+1]])
        Xvar<-data.frame(Xvar=round(sum(test.pls.results$Xvar[[tmp.OSC.comp+1]])*100,1))
        if(!is.null(test.pls.results$OPLSDA.stats)){oplsda.stats<-data.frame(test.pls.results$OPLSDA.stats[[tmp.OSC.comp+1]])} else {oplsda.stats<-NULL}
        
        #results
        predicted.y<-rbind(as.matrix(train.pred),as.matrix(test.pred))
        actual.y<-rbind(as.matrix(train.real),as.matrix(test.real))
        test.index<-pls.train.index
        if(is.null(oplsda.stats)){
            res<-list(data.frame(predicted = predicted.y[back.sort,], actual = actual.y[back.sort,],train.test.id=test.index), data.frame(Xvar,Q2,RMSEP))
        } else {
            res<-list(data.frame(predicted = predicted.y[back.sort,], actual = actual.y[back.sort,],train.test.id=test.index), data.frame(Xvar,Q2,RMSEP,oplsda.stats))
        }
        return(res)
    })
    
    #need to summarize results
    id<-c(1:length(results))[c(1:length(results))%%2==0]
    aggregated<-do.call("rbind",lapply(1:length(results),function(i){data.frame(results[[i]][2])}))
    
    aggregated.summary<-matrix(paste(signif(apply(aggregated,2,mean,na.rm=TRUE),4),"+/-",signif(apply(aggregated,2,sd,na.rm=TRUE),3)),nrow=1)
    colnames(aggregated.summary)<-colnames(aggregated)
    list(full.results=results, performance=aggregated, summary=aggregated.summary)
}

#permute OSC-PLS model (in progress version for multiple Ys)
permute.OSC.PLS2<-function(data,y,n=10,ncomp,OSC.comp=1,train.test.index=NULL,...){ # should be made parallel
    
    
    #permuted Y
    perm.y<-lapply(1:n,function(i)
    {
        apply(y,2,gtools::permute)
    })
    
    #collect correlation between y and permuted y\
    # disa bled for multi y until main fxn can handle these correctly
    if(ncol(y)==1){
        cor.with.y<-data.frame(correlation=abs(cor(cbind(y,do.call("cbind",perm.y))))[-1,1])
    } else {
        cor.with.y<-NULL
    }
    
    #generate permuted models
    model<-lapply(1:n,function(i)
    {
        # cat("permuting model",i,"\n")
        if(!is.null(train.test.index)) {tmp.train.test.index<-train.test.index[,i,drop=FALSE]} else {tmp.train.test.index<-train.test.index}
        model<-make.OSC.PLS.model(pls.y=as.matrix(perm.y[[i]]),pls.data=data,comp=ncomp,OSC.comp=OSC.comp,train.test.index=tmp.train.test.index,...) #
        #get stats
        q2<-model$Q2[[OSC.comp+1]][ncomp+1,,drop=FALSE]# cv adjusted
        rx2<-round(sum(model$Xvar[[OSC.comp+1]])*100,1)
        pred.val<-as.matrix(model$fitted.values[[OSC.comp+1]][,,ncomp])
        rmsep<-model$rmsep[[OSC.comp+1]]# take CV adjusted internal RMSEP (see true RMSEP below)
        if(!is.null(train.test.index)) {
            rmsep<-model$predicted.RMSEP[[OSC.comp+1]]
        }
        list(RX2=rep(rx2,length(q2)),Q2=q2,RMSEP=rmsep)#,predicted=pred.val,actual=perm.y[[i]])
    })
    
    tmp<-matrix(unlist(do.call("rbind",model)),ncol=ncol(y)*3 )
    colnames(tmp)<-paste0("Y.",1:ncol(y),"_",rep(c("Xvar","Q2","RMSEP"),each=ncol(y)))
    means<-apply(tmp,2,mean, na.rm=TRUE)
    sds<-apply(tmp,2,sd,na.rm=TRUE)
    summary<-paste(signif(means,4),"±", signif(sds,3))
    return(list(permuted.values=cbind(tmp,cor.with.y), mean = means, standard.deviations = sds, summary = summary))
}

#statistical test to compare permuted distribution to model performance
OSC.validate.model<-function(model, perm, train= NULL, test="t.test",...) {
    #model is an object generated with get.OSC.model
    #perm must be object generated with permute.OSC.PLS
    # if train = NULL  and test = "t.test" perform a one-sample t-test to test if model stat comes from permuted distribution
    # else perform a two sample t-test to compare train/test to permuted stats
    # if test= "perm.test" calculate p-value based on http://www.ncbi.nlm.nih.gov/pubmed/21044043
    # number tests different than the permuted +1/ number of permutations + 1
    
    #match model and perm objects
    perm.vals<-perm$permuted.values
    if(is.null(perm.vals)) perm.vals<-perm$performance # ugly OOP, switching between two versions to get perm
    
    comp<-length(model$Q2)
    if(is.null(train)){
        if(any(colnames(perm.vals)=="AUC")){
            mod.vals<-data.frame(RX2 = cumsum(model$Xvar*100)[comp-1],
            Q2 = model$Q2[comp],
            RMSEP = model$RMSEP[comp],model$OPLSDA.stats)[1,]
        } else {
            mod.vals<-data.frame(RX2 = cumsum(model$Xvar*100)[comp-1],
            Q2 = model$Q2[comp],
            RMSEP = model$RMSEP[comp])[1,]
        }
    } else {
        mod.vals<-data.frame(train$performance) # from train object
    }
    
    #test single model value (should test log distributions to avoid the effect of outliers?)
    single.test<-function(mod,perm){
        data.frame(matrix(tryCatch(t.test(perm,mu=unlist(mod))$p.value, error=function(e) {1}),ncol=1))
    }
    
    #two group test
    group.test<-function(mod,perm){
        data.frame(matrix(tryCatch(t.test(mod,perm)$p.value, error=function(e) {1}),ncol=1))
    }
    
    if(is.null(train)){
        p.vals<-lapply(1:ncol(mod.vals),function(i){
            if(names(mod.vals[i])%in%c("RMSEP","ER","FPR")) {dir<-">"} else {dir<-"<"} # used for permutation tests
            switch(test,
            "t.test"        = single.test(mod=mod.vals[i],perm=perm.vals[,i]),
            "perm.test" = perm.test(mod=mod.vals[i],perm=perm.vals[,i],dir,type=1),
            "perm.test2" = perm.test(mod=mod.vals[i],perm=perm.vals[,i],dir,type=2))
        })
        
        if(is.null(perm$summary)){perm$summary<-"not permuted"}
        #make output in table form
        res<-data.frame(rbind(signif(mod.vals,4),perm$summary,unname(signif(unlist(p.vals),4))))
        rownames(res)<-c("model","permuted model","p-value")
    } else {
        p.vals<-lapply(1:ncol(mod.vals),function(i){
            if(names(mod.vals[i])%in%c("RMSEP","ER","FPR")) {dir<-">"} else {dir<-"<"} # used for permutation tests
            switch(test,
            "t.test"        = group.test(mod=mod.vals[,i],perm=perm.vals[,i]),
            "perm.test" = perm.test(mod=mod.vals[i],perm=perm.vals[,i],dir),
            "perm.test2" = perm.test(mod=mod.vals[i],perm=perm.vals[,i],dir,type=2))
            
        })
        
        if(is.null(perm$summary)){perm$summary<-"not permuted"}
        #make output in table form
        res<-data.frame(rbind(train$summary,perm$summary,unlist(signif(unlist(p.vals),4))))
        rownames(res)<-c("model","permuted model","p-value")
    }
    return(res)
}

#conservative p-value based on permutation tests  http://www.ncbi.nlm.nih.gov/pubmed/21044043 (could go elsewhere)
perm.test<-function(mod,perm,compare="<",type=1){
    f<-function(mod,perm,compare){tryCatch((sum(do.call(compare,list(na.omit(mod),na.omit(perm))))+1)/(mean(length(na.omit(perm)),length(na.omit(mod)))+1),error=function(e) {NA})}
    if(type==1){f(mod,perm,compare)} else {f(mod,perm,compare=compare)-1/(mean(length(na.omit(perm)),length(na.omit(mod)))+1)}
}

#compare train stats between two models
OSC.PLS.model.compare<-function(model1, model2,test="t.test",...){
    #models must be object generated with OSC.PLS.train.test
    
    p.vals<-do.call("cbind",lapply(1:ncol(model1$performance), function(i) {
        if(colnames(model1$performance)[i]%in%c("RMSEP","ER","FPR")) {dir<-">"} else {dir<-"<"}
        switch(test,
        
        "t.test"     = tryCatch(t.test(model1$performance[,i],model2$performance[,i])$p.value, error=function(e) {1}), #force error = insiginificant
        "perm.test" = perm.test(model1$performance[,i],model2$performance[,i],compare=dir),
        "perm.test2" = perm.test(model1$performance[,i],model2$performance[,i],compare=dir,type=2)
        )
    })
    )
    
    res<-data.frame(rbind(model1$summary,model2$summary,signif(p.vals,4))) # don't include Xvar
    dimnames(res)<-list(c("model1","model2","p-value"),colnames(model1$summary))
    return(res)
}

#conduct train/test validations on PLS model
PLS.train.test<-function(pls.data,pls.y,pls.train.index,comp,...) {
    #only test first Y
    
    pls.y<-as.matrix(pls.y)
    #order for merging test with train stats in one object
    new.order<-c(c(1:nrow(pls.data))[pls.train.index=="train"],c(1:nrow(pls.data))[pls.train.index=="test"])
    back.sort<-order(new.order)
    
    train.y<-train.real<-pls.y[pls.train.index=="train",]
    train.data<-pls.data[pls.train.index=="train",]
    #all arguments gave been preset elsewhere
    test.pls.results<-make.PLS.model(train.y,train.data,ncomp=comp,...)
    Q2<-unlist(R2(test.pls.results)$val[,,max(comp)+1])
    
    train.pred<-test.pls.results$fitted.values[,,comp]
    
    #use model to predict test values
    test.data<-as.matrix(pls.data[pls.train.index=="test",])
    test.pred<- as.data.frame(predict(object=test.pls.results, newdata=test.data, ncomp = c(1:comp), comps=c(1:comp),type ="response"))
    
    test.real<-pls.y[pls.train.index=="test",]
    RMSEP<-sapply(1:ncol(pls.y), function(i){
        (sum((test.pred[,i]-test.real[,i])^2)/nrow(test.pred))^.5
    })
    
    #results
    predicted.y<-rbind(train.pred,test.pred)
    actual.y<-rbind(train.real,test.real)
    test.index<-pls.train.index
    res<-list(predicted.y[back.sort,], actual.y[back.sort,], test.index,RMSEP,Q2,LVs=PCs)
    names(res)<-c("predicted.y","actual.y=","pls.train.index","RMSEP","Q2","LVs")
    return(res)
}

#conduct train/test validations on O-PLS model
OSC.PLS.train.test<-function(pls.data,pls.y,train.test.index,comp,OSC.comp,...) {
    pls.y<-as.matrix(pls.y)
    results<-lapply(1:ncol(train.test.index), function(i){
        pls.train.index<-as.matrix(train.test.index[,i])
        #order for merging test with train stats in one object
        new.order<-c(c(1:nrow(pls.data))[pls.train.index=="train"],c(1:nrow(pls.data))[pls.train.index=="test"])
        back.sort<-order(new.order)
        
        train.y<-train.real<-pls.y[pls.train.index=="train",]
        train.data<-pls.data[pls.train.index=="train",]
        test.real<-pls.y[pls.train.index=="test",]
        #all arguments gave been preset elsewhere
        test.pls.results<-make.OSC.PLS.model(pls.y=pls.y,pls.data=pls.data,comp=comp,OSC.comp=OSC.comp, train.test.index=pls.train.index,...)
        tmp.OSC.comp<-max(test.pls.results$OSC.LVs) # control when limited with Orthogonal dimensions
        Q2<-data.frame(Q2=test.pls.results$Q2[[tmp.OSC.comp+1]][comp,])
        
        #fitted values
        train.pred<-test.pls.results$fitted.values[[tmp.OSC.comp+1]][,,comp]
        test.pred<-test.pls.results$predicted.Y[[tmp.OSC.comp+1]]
        
        RMSEP<-data.frame(RMSEP=test.pls.results$predicted.RMSEP[[tmp.OSC.comp+1]])
        Xvar<-data.frame(Xvar=round(sum(test.pls.results$Xvar[[tmp.OSC.comp+1]])*100,1))
        if(!is.null(test.pls.results$OPLSDA.stats)){oplsda.stats<-data.frame(test.pls.results$OPLSDA.stats[[tmp.OSC.comp+1]])} else {oplsda.stats<-NULL}
        
        #results
        predicted.y<-rbind(as.matrix(train.pred),as.matrix(test.pred))
        actual.y<-rbind(as.matrix(train.real),as.matrix(test.real))
        test.index<-pls.train.index
        if(is.null(oplsda.stats)){
            res<-list(data.frame(predicted = predicted.y[back.sort,], actual = actual.y[back.sort,],train.test.id=test.index), data.frame(Xvar,Q2,RMSEP))
        } else {
            res<-list(data.frame(predicted = predicted.y[back.sort,], actual = actual.y[back.sort,],train.test.id=test.index), data.frame(Xvar,Q2,RMSEP,oplsda.stats))
        }
        return(res)
    })
    
    #need to summarize results
    id<-c(1:length(results))[c(1:length(results))%%2==0]
    aggregated<-do.call("rbind",lapply(1:length(results),function(i){data.frame(results[[i]][2])}))
    
    aggregated.summary<-matrix(paste(signif(apply(aggregated,2,mean,na.rm=TRUE),4),"+/-",signif(apply(aggregated,2,sd,na.rm=TRUE),3)),nrow=1)
    colnames(aggregated.summary)<-colnames(aggregated)
    list(full.results=results, performance=aggregated, summary=aggregated.summary)
}

# function for splitting dataset in to test and trainning sets
test.train.split<-function(nsamples, n=1, strata=NULL, prop.train=2/3, split.type="random",data=NULL){
    #nsamples the number of samples in the data set to split among test and trainnig data sets
    #n the number or test/training splits to return
    #strata factor within whose levels the test/trainning sets will be derived
    #prop.train the proportion of samples in the trainning set
    #split.type how sample assignment to trainning/test splits is determined,  the options are "random", "duplex"
    #data needed for duplex method
    res<-lapply(1:n,function(i){
        if(is.null(strata)){
            if(split.type=="random"){
                t.num<-ceiling(nsamples*prop.train)
                test.train.id<-rep("test",nsamples)
                test.train.id[sample(c(1:length(test.train.id)),t.num)]<-"train"
            }
            
            if(split.type=="duplex"){ # takes the floor, if the sample number is not even the proportion of trainning samples may be one less than specified
                object<-ken.sto2(data, per = "TRUE",per.n= (1-prop.train),va = "TRUE",num=1)
                object<-duplex.select(data,object,percent.in.test=(1-prop.train))
                test.train.id<-rep("train",nrow(data))
                test.train.id[object$`Chosen validation row number`]<-"test"
            }
        } else {
            if(split.type=="random"){
                tmp<-split(1:nsamples,strata)
                train.id<-unlist(sapply(1:length(tmp),function(i){
                    t.num<-ceiling(length(tmp[[i]])*prop.train)
                    tmp[[i]][sample(c(1:length(tmp[[i]])),t.num)]
                }))
                test.train.id<-rep("test",nsamples)
                test.train.id[train.id]<-"train"
            }
            
            if(split.type=="duplex"){ # trainning/test assignments are underestimated for odd number of samples (due rounding down)
                tmp<-split(data,strata)
                test.id<-fixlc(sapply(1:length(tmp),function(i){
                    object<-ken.sto2(tmp[[i]], per = "TRUE",per.n= (1-prop.train),va = "TRUE",num=1)
                    object<-duplex.select(tmp[[i]],object,percent.in.test=(1-prop.train))
                    object$`Chosen validation sample names`
                }))
                test.train.id<-rep("train",nrow(data))
                test.train.id[rownames(data)%in%test.id]<-"test"
            }
            
        }
        test.train.id
    })
    as.data.frame(do.call("cbind",res))
}

#function for carrying out test/trainning split based on duplex or kennard-stone method
ken.sto2<-function(inp, per = "TRUE", per.n = 0.3, num = 7, va = "TRUE"){
    #based on  ken.sto in package "soil.spec"
    #changes: altered number of PCs selection
    #took out saving, plotting
    #opened slot for custom PCA analysis options
    #fixed some bugs
    
    if (class(inp) != "data.frame" & class(inp) != "matrix")
    {
        stop("Invalid argument: 'inp' has to be of class 'data.frame' or 'matrix'.")
    }
    if (per != "TRUE" & per != "FALSE")
    {
        stop("Invalid argument: 'per' has to be either 'TRUE' or 'FALSE'.")
    }
    
    if (per == "TRUE")
    {
        if (class(per.n) != "numeric")
        {
            stop("Invalid argument: 'per' has to be of class 'numeric'.")
        }
        
        if (per.n < 0 | per.n > 1)
        {
            stop("Invalid argument: 'per' has to be between 0 and 1.")
        }
        n <- round(per.n * nrow(inp), 0)
    }
    
    if (per == "FALSE")
    {
        if (class(as.integer(num)) != "integer")
        {
            stop("Invalid argument: 'num' has to be of class 'integer'.")
        }
        
        if (num <= 0)
        {
            stop("Invalid argument: 'num' has to be between 1 and the number of samples minus one.")
        }
        
        if (num >= nrow(inp))
        {
            stop("Invalid argument: 'num' has to be between 1 and the number of samples minus one.")
        }
        
        n <- num
    }
    if (va != "TRUE" & va != "FALSE")
    {
        stop("Invalid argument: 'va' has to be either 'TRUE' or 'FALSE'.")
    }
    
    #allow to use specified PCA analysis results
    pca <- prcomp(inp, scale = T)
    prco <- as.data.frame(pca$x)
    cpv <- summary(pca)[[6]][3,]
    zzz <- matrix(nrow = 1, ncol = (length(cpv)-2))
    for (i in 1:(ncol(zzz)-2))
    {
        e <- (cpv[i] + 0.04) < cpv[i + 3]
        zzz[i] <- e
    }
    pc <- (which(zzz == FALSE) - 1)[1]
    if (pc == 1|is.na(pc))
    {
        pc <- 2
    }
    prco <- prco[, 1:pc]
    min <- c(rep(1, ncol(prco)))
    max <- c(rep(1, ncol(prco)))
    for (i in 1:ncol(prco))
    {
        blub <- which(prco[, i] == min(prco[, i]))
        min[i] <- blub[1]
        bla <- which(prco[, i] == max(prco[, i]))
        max[i] <- bla[1]
    }
    min <- rownames(prco)[min]
    max <- rownames(prco)[max]
    start <- unique(c(min, max))
    start.n <- match(start, rownames(inp))
    
    if (va == "FALSE")
    {
        euc <- as.data.frame(as.matrix(dist(prco)))
        inp.start <- rownames(prco)[-start.n]
        inp.start.b <- inp.start
        cal <- start
        stop<-min(c(n,length(start)))
        for (k in 1:(stop))
        {
            test <- apply(euc[inp.start.b, cal], 1, min)
            bla <- names(which(test == max(test)))
            cal <- c(cal, bla)
            inp.start.b <- inp.start.b[-(which(match(inp.start.b,
            bla) != "NA"))]
        }
        cal.n <- match(cal, rownames(inp))
        
        output <- list(`Calibration and validation set` = va,
        `Number important PC` = pc, `PC space important PC` = prco,
        `Chosen sample names` = unique(cal), `Chosen row number` = unique(cal.n),
        `Chosen calibration sample names` = "NULL", `Chosen calibration row number` = "NULL",
        `Chosen validation sample names` = "NULL", `Chosen validation row number` = "NULL")
    }
    
    if (va == "TRUE")
    {
        n<-ceiling(per.n*nrow(inp))
        cal.start <- start
        cal.start.n <- start.n
        val.min <- c(rep(1, ncol(prco)))
        val.max <- c(rep(1, ncol(prco)))
        for (i in 1:ncol(prco))
        {
            blub <- which(prco[-cal.start.n, i] == min(prco[-cal.start.n,
            i]))
            val.min[i] <- blub[sample(length(blub), 1)]
            bla <- which(prco[-cal.start.n, i] == max(prco[-cal.start.n,
            i]))
            val.max[i] <- bla[sample(length(bla), 1)]
        }
        val.min <- rownames(prco[-cal.start.n, ])[val.min]
        val.max <- rownames(prco[-cal.start.n, ])[val.max]
        val.start <- unique(c(val.min, val.max))
        val.start.n <- match(val.start, rownames(inp))
        cal.val <- c(cal.start, val.start)
        cal.val.start <- match(c(cal.start, val.start), rownames(inp))
        euc <- as.data.frame(as.matrix(dist(prco)))
        inp.start <- rownames(prco)[-cal.val.start]
        inp.start.b <- inp.start
        val <- val.start
        stop<-n#min(c(n,length(val.start)))
        k<-1
        for (k in 1:(stop))
        {
            test <- apply(euc[inp.start.b, val], 1, min)
            bla <- names(which(test == max(test)))
            val <- c(val, bla)
            inp.start.b <- inp.start.b[-(which(match(inp.start.b,
            bla) != "NA"))]
        }
        val.n <- match(val, rownames(inp))
        cal.n <- c(1:nrow(inp))[-val.n]
        cal <- rownames(inp)[cal.n]
        
        #tmp fix for problem in function
        n<-ceiling(per.n*nrow(inp))
        if(n<1){n<-1}
        tst.id<-unique(val.n)
        if(n>length(tst.id)){n=length(tst.id)}
        val.n<-sample(tst.id,n)
        cal.n<-c(cal.n,tst.id[!tst.id%in%val.n])
        
        cal <- rownames(inp)[cal.n]
        val<-rownames(inp)[val.n]
        
        
        
        output <- list(`Calibration and validation set` = va,
        `Number important PC` = pc, `PC space important PC` = prco,
        `Chosen sample names` = "NULL", `Chosen row number` = "NULL",
        `Chosen calibration sample names` = unique(cal), `Chosen calibration row number` = unique(cal.n),
        `Chosen validation sample names` = unique(val), `Chosen validation row number` = unique(val.n))
    }
    class(output) <- "ken.sto"
    return(output)
}

#wrapper to iterate ken.sto2
duplex.select<-function(data,ken.sto2.obj,percent.in.test){
    #determine how many more are needed
    start.have<-ken.sto2.obj$`Chosen validation sample names`
    need<-percent.in.test*nrow(data)-length(start.have)
    
    #don't do anything if there are enough
    if(need>0)
    {
        #extract from remainning data
        have<-start.have
        while(need>0)
        {
            tmp.data<-data[!rownames(data)%in%have,]
            more<-ken.sto2(tmp.data, per = "TRUE", per.n = percent.in.test, num = 7, va = "TRUE")
            now.have<-more$`Chosen validation sample names`
            need<-percent.in.test*nrow(data)-(length(now.have)+length(have))
            have<-c(have,now.have)
        }
        
        #adjust for too many selected
        drop<-NA
        if(need<0)
        {
            drop<-now.have[sample(length(now.have),abs(need))]
        }
        
        new.obj<-have[!have%in%drop]
        
        #objects to return
        `Chosen validation sample names`=c(new.obj)
        `Chosen validation row number`= c(1:nrow(data))[rownames(data)%in%new.obj]
        `Chosen calibration sample names`= rownames(data)[!rownames(data)%in%`Chosen validation sample names`]
        `Chosen calibration row number` =c(1:nrow(data))[rownames(data)%in%`Chosen calibration sample names`]
    }else{
        `Chosen validation sample names`=ken.sto2.obj$`Chosen validation sample names`
        `Chosen validation row number`= ken.sto2.obj$`Chosen validation row number`
        `Chosen calibration sample names`= ken.sto2.obj$`Chosen calibration sample names`
        `Chosen calibration row number` =ken.sto2.obj$`Chosen calibration row number`
    }
    output<-list(`Chosen validation row number`= `Chosen validation row number`,
    `Chosen validation sample names`=`Chosen validation sample names`,
    `Chosen calibration sample names` = `Chosen calibration sample names`,
    `Chosen calibration row number` = `Chosen calibration row number`)
}

#function to calculate included/excluded feature model stats
optimize.OPLS.feature.select<-function(model,feature.subset,permute=TRUE,train.test.index,progress=TRUE,test="perm.test",...){
    #need to know OPLS model args (*later store in model and use this a reference)
    #not sure why this is missing
    ntests<-ncol(train.test.index)
    #selected model stats (TODO: avoid recalculating the model!)
    data<-model$data[[1]][,feature.subset,drop=F]
    model1<-make.OSC.PLS.model(pls.y=model$y[[1]],pls.data=data,comp=model$total.LVs[1],OSC.comp=max(model$OSC.LVs), validation = model$model.description$validation,method=model$model.description$method, cv.scale=model$model.description$cv.scale,return.obj="stats",progress = progress,...)
    
    if(permute==TRUE){
        #permutation
        sel.permuted.stats <-permute.OSC.PLS.train.test(data,pls.y = model$y[[1]],perm.n = ntests,train.test.index = train.test.index,comp =  model$total.LVs[1],OSC.comp = max(model$OSC.LVs),progress = progress,...)
        # sel.permuted.stats <- permute.OSC.PLS(data = data, y = model$y[[1]], n = ntests, ncomp = model$total.LVs[1], osc.comp=max(model$OSC.LVs), progress = progress, train.test.index = train.test.index,...) #...
    } else {
        sel.permuted.stats<-NULL
    }
    
    #training/testing to get robust model stats
    sel.OPLS.train.stats <- OSC.PLS.train.test(pls.data = data, pls.y = model$y[[1]], train.test.index, comp = model$total.LVs[1], OSC.comp = max(model$OSC.LVs), cv.scale = model$model.description$cv.scale, progress = progress,...) # ...
    sel.OPLS.model<-OSC.validate.model(model = model1, perm = sel.permuted.stats, train = sel.OPLS.train.stats,test,...)
    
    #excluded model stats
    data<-model$data[[1]][,!feature.subset]
    model2<-make.OSC.PLS.model(pls.y=model$y[[1]],pls.data=data,comp=model$total.LVs[1],OSC.comp=max(model$OSC.LVs), validation = model$model.description$validation,method=model$model.description$method, cv.scale=model$model.description$cv.scale,return.obj="stats",progress = progress,...)
    if(permute==TRUE){
        #permutation
        ex.permuted.stats <-permute.OSC.PLS.train.test(data,pls.y = model$y[[1]],perm.n = ntests,train.test.index = train.test.index,comp =  model$total.LVs[1],OSC.comp = max(model$OSC.LVs),progress = progress,...)
        # ex.permuted.stats <- permute.OSC.PLS(data = data, y = model$y[[1]], n = ntests, ncomp = model$total.LVs[1], osc.comp=max(model$OSC.LVs), progress = progress, train.test.index = train.test.index,...) #...
    } else {
        ex.permuted.stats<-NULL
    }
    
    #training/testing to get robust model stats
    ex.OPLS.train.stats <- OSC.PLS.train.test(pls.data = data, pls.y = model$y[[1]], train.test.index, comp = model$total.LVs[1], OSC.comp = max(model$OSC.LVs), cv.scale = model$model.description$cv.scale, progress = progress,...) # ...
    ex.OPLS.model<-OSC.validate.model(model = model2, perm = ex.permuted.stats, train = ex.OPLS.train.stats,test,...)
    
    full.sel.model.comparison<-OSC.PLS.model.compare(model1=sel.OPLS.train.stats, model2=ex.OPLS.train.stats,test,...)
    #create final table
    out<-data.frame(cbind(model=c(rep("selected",3),rep("excluded",3),"comparison"),rbind(as.matrix(sel.OPLS.model),as.matrix(ex.OPLS.model),as.matrix(full.sel.model.comparison)[3,,drop=F])))
    
    list(selected.train=sel.OPLS.train.stats,selected.permuted=sel.permuted.stats,excluded.train=ex.OPLS.train.stats,excluded.permuted=ex.permuted.stats,summary=out)
    
}

#get classification performance statistics
O.PLS.DA.stats<-function(truth,pred){
    
    #
    check.get.packages("ROCR")
    # library(ROCR)
    # # library(caret) #need e1071
    # library(hmeasure)
    
    
    y.range<-range(as.numeric(truth))
    mid<-mean(y.range)
    binned.pred<-pred
    binned.pred[binned.pred<mid]<-y.range[1]
    binned.pred[binned.pred>=mid]<-y.range[2]
    # scaled.pred<-rescale(as.numeric(pred),y.range)
    # scaled.pred[scaled.pred<mid]<-y.range[1]
    # scaled.pred[scaled.pred>=mid]<-y.range[2] # not sure what to do with a prediction == the mid point
    #get AUC
    mod.AUC<-function(pred,truth){
        # pred1 <- prediction(pred, truth)
        #perf <- performance(pred1, measure="tpr", x.measure="fpr")
        # plot(perf,lty=1,lwd=4,col="#9400D350") # plot not interesting with so  few measurements
        # add precision recall http://stackoverflow.com/questions/8499361/easy-way-of-counting-precision-recall-and-f1-score-in-r
        unlist(performance(prediction(pred, truth),measure= "auc")@y.values)
    }
    
    
    #misclassCounts(binned.pred,truth)  # library(hmeasure)
    
    AUC<-tryCatch(mod.AUC(pred=binned.pred,truth=truth),error=function(e){NA}) # protect errors due to !=2 groups
    # AUC<-mod.AUC(pred=binned.pred,truth=truth) # get NA when groups !=2
    # get other metrics
    # library(hmeasure) # using modified fxn which accepts inputs other than 1 and 0
    results<-tryCatch(misclassCounts2(pred=binned.pred,truth=truth)$metrics ,error=function(e){NULL}) # protect errors due to >2 groups or use caret::confusionMatrix
    #happens when there are perfect predictions
    # library(caret)
    # results<-tryCatch(confusionMatrix(binned.pred,as.numeric(truth)),error=function(e){"error"}) # protect errors due to >2 groups
    if(is.null(results)){results<-list();results$byClass[1:2]<-NA}
    # res<-data.frame(AUC=AUC,sensitivity=results$byClass[1], specificity=results$byClass[2])
    res<-data.frame(AUC=AUC,results)
    
    rownames(res)<-"model"
    return(res)
}

#modified hmeasure::misclassCounts to accept class values besides 0 and 1
# limited to only 2 classes
misclassCounts2<-function (pred,truth){
    
    vals<-sort(unique(truth),decreasing=TRUE) # smaller value is not the class
    TP <- sum(pred == vals[1] & truth == vals[1])
    FP <- sum(pred == vals[1] & truth == vals[2])
    TN <- sum(pred == vals[2] & truth == vals[2])
    FN <- sum(pred == vals[2] & truth == vals[1])
    
    conf.matrix <- data.frame(pred.1 = c(TP, FP), pred.0 = c(FN,
    TN))
    row.names(conf.matrix) <- c("actual.1", "actual.0")
    ER <- (FP + FN)/(TP + FP + TN + FN)
    Sens <- TP/(TP + FN)
    Spec <- TN/(TN + FP)
    Precision <- TP/(TP + FP)
    Recall <- Sens
    TPR <- Recall
    FPR <- 1 - Spec
    F <- 2/(1/Precision + 1/Sens)
    Youden <- Sens + Spec - 1
    metrics <- data.frame(ER = ER, Sens = Sens, Spec = Spec,
    Precision = Precision, Recall = Recall, TPR = TPR, FPR = FPR,
    F = F, Youden = Youden)
    return(list(conf.matrix = conf.matrix, metrics = metrics))
}

# generic mean squared error of prediction
.MSEP<-function(actual,pred){
    if(is.null(dim(actual))){
        mean((actual-pred)^2)
    } else {
        sapply(1:ncol(actual),function(i){mean((actual[,i]-pred[,i])^2)})
    }
}

#wrapper to carry out multiple feature selection and model validation runs
multi.OPLS.feature.select<-function(model,filter,ocomp=max(model$OSC.LVs),extra=NULL,train.test.index=NULL,progress=TRUE,...){
    
    
    library(plyr)
    
    .model<-get.OSC.model(obj=model,OSC.comp=ocomp)
    
    optimal.feature.selection<-lapply(1:length(filter),function(i){
        top<-filter[i]
        #this step can be done only once
        opts<-PLS.feature.select(model$data[[1]],pls.scores=.model$scores[,][,1,drop=F],pls.loadings=.model$loadings[,][,1,drop=F],pls.weight=.model$loadings[,][,1,drop=F],plot=FALSE,p.value=1,FDR=FALSE,cut.type="number",top=top,separate=FALSE)
        optim<-optimize.OPLS.feature.select(model=model,feature.subset=opts$combined.selection,permute=TRUE,train.test.index=train.test.index,progress=progress,...) # check variance explained in X
        if(progress==TRUE){print(paste0(round(i/length(filter)*100,0)," % complete"))}
        # return model and permuted results
        rbind(cbind(type="model",rbind(data.frame(vars=top,model="included",optim$selected.train$performance),data.frame(vars=top,model="excluded",optim$excluded.train$performance))),
        cbind(type="permuted",rbind(data.frame(vars=top,model="included",optim$selected.permuted$performance),data.frame(vars=top,model="excluded",optim$selected.permuted$performance))))
    })
    
    #summary
    obj<-do.call("rbind",optimal.feature.selection)
    
    #get summary
    library(plyr)
    means<-ddply(obj,.(type,vars,model),colwise(mean,na.rm=TRUE))
    sds<-ddply(obj,.(type,vars,model),colwise(sd,na.rm=TRUE))
    
    #calculate p-values for the comparison
    
    #fxn
    compare.mods<-function(mod.vals,perm.vals,test="perm.test",...){
        lapply(1:ncol(mod.vals),function(i){
            if(names(mod.vals[i])%in%c("RMSEP","ER","FPR")) {dir<-">"} else {dir<-"<"} # used for permutation tests
            switch(test,
            "t.test"        = group.test(mod=mod.vals[,i],perm=perm.vals[,i]),
            "perm.test" = perm.test(mod=mod.vals[i],perm=perm.vals[,i],dir),
            "perm.test2" = perm.test(mod=mod.vals[i],perm=perm.vals[,i],dir,type=2))
            
        })
    }
    
    #hack objects to fit
    #included to excluded
    tmp.l<-split(obj,obj$type)
    
    #included vs. excluded
    p.vals<-ddply(tmp.l[["model"]],.(vars),function(tmp,...){
        tmp.obj<-tmp[,-c(1:3)]
        id<-tmp$model=="included"
        res<-data.frame(compare.mods(tmp.obj[id,],tmp.obj[!id,],...))
        colnames(res)<-colnames(tmp)[-c(1:3)]
        return(res)
    })
    
    #model vs. permuted
    tmp.l<-split(obj,obj$model)
    #included vs. permuted
    p.vals2<-ddply(tmp.l[["included"]],.(vars),function(tmp,...){
        tmp.obj<-tmp[,-c(1:3)]
        id<-tmp$type=="model"
        res<-data.frame(compare.mods(tmp.obj[id,],tmp.obj[!id,],...))
        colnames(res)<-colnames(tmp)[-c(1:3)]
        return(res)
    })
    
    #excluded vs. permuted
    #included vs. permuted
    p.vals3<-ddply(tmp.l[["excluded"]],.(vars),function(tmp,...){
        tmp.obj<-tmp[,-c(1:3)]
        id<-tmp$type=="model"
        res<-data.frame(compare.mods(tmp.obj[id,],tmp.obj[!id,],...))
        colnames(res)<-colnames(tmp)[-c(1:3)]
        return(res)
    })
    #p.values
    p.values<-data.frame(rbind(cbind(type="included vs. excluded",p.vals),cbind(type="included vs. permuted",p.vals2),cbind(type="excluded vs. permuted",p.vals3)))
    
    list(all.results=obj,summary=list(mean=means,sd=sds,p.values=p.values))
}

#plot the results of multi.OPLS.feature.select ###modified by Kai
plot.multi.OPLS.feature.select<-function(object,extra=NULL,objects=c("RMSEP","Q2","AUC","F","Sens","Spec","Precision","Recall")){
    check.get.packages(c("gridExtra","ggplot2"))
    
    means<-object$summary$mean[,-c(1:3)]
    sds<-object$summary$sd[,-c(1:3)]
    p.vals<-object$summary$p.values[,-c(1:2)]
    model<-object$summary$mean$model
    vars<-factor(object$summary$mean$vars) # make this factor for equal spacin
    model.type<-paste(model,object$summary$mean$type,sep="_")
    
    #filter error columns
    f<-function(x){sum(is.na(x))}
    id<-!apply(means,2,f)==nrow(means)&colnames(means)%in%objects #exclude all NA columns
    means<-means[,id,drop=FALSE]
    sds<-sds[,id,drop=FALSE]
    p.vals<-p.vals[,id,drop=FALSE] ###fixed by Kai
    plot.list<-list()
    k<-1
    for(i in 1:ncol(means)){
        value<-colnames(means)[i]
        tmp<-data.frame(value=means[,value],error=sds[,value],model=model.type,vars=vars)
        vlines<-p.vals$var[p.vals[,value]<=0.05]
        if(length(vlines)>0){show.sig<-geom_vline(xintercept=vlines,linetype="dotted")}    else {show.sig<-NULL}
        plot.list[[k]]<-ggplot(data = tmp, aes(x = vars, y = value,group=model) ) +
        geom_errorbar(aes(ymin = value - error, ymax = value + error, color=model), size=1, width=.1) + # add error bars (do so before geom_point so the points are on top of the error bars)
        geom_line(aes(color=model),size=1) + # join points with lines (specify this before geom_point, or the lines will be drawn over the shapes)
        geom_point(aes(shape=model,color=model), size=5) +ylab(value) +
        xlab("cutoff") + theme_light() + show.sig + extra # scale_x_continuous("cutoff", breaks=unique(tmp$vars))
        k<-k+1
    }
    #plot
    if(length(plot.list)>1){
        do.call("grid.arrange",plot.list)
    } else {print(plot.list[[1]])}
}


#extract best model
best.OPLS.features<-function(obj=res,decreasing = c("RMSEP"),measures=c("AUC","Sens","Spec","Precision","Recall","TPR","F","Youden")){
    
    #rank
    tmp<-obj$summary$mean[obj$summary$mean$model=="included"&obj$summary$mean$type=="model",,drop=FALSE]
    .tmp<-tmp[,colnames(tmp)%in%measures,drop=FALSE]
    # variables which need small is best rank need to be inverted
    .tmp[,colnames(.tmp)%in%decreasing]<-1/.tmp[,colnames(.tmp)%in%decreasing]
    .rank<-apply(.tmp,2,rank) #rank is decreasing
    rmat<-matrix(.rank,,length(measures))
    rowid<-matrix(1:nrow(rmat),nrow(rmat),ncol(rmat))[which.max(rmat)]
    tmp[rowid,,drop=FALSE]
}
### sort create by Kai
sortchar<-function(x,order=FALSE){
    tmp<-gsub('[a-zA-Z]','',x)
    tmp<-gsub('[,\\_@\\*\\&\\^\\$\\#\\~]','',tmp)
    tmp<-as.numeric(tmp)
    names(tmp)<-x
    if(order==TRUE){
        return(order(tmp))
    }else{
        return(sort(tmp))
    }
}
#iteratively split a vector into fractional units taking the floor
split.vector<-function(n,step=.5){
    
    #adapted from randomForest::rfcv
    
    k <- floor(log(n, base = 1/step))
    n.var <- round(n * step^(0:(k - 1)))
    same <- diff(n.var) == 0
    if (any(same))
    n.var <- n.var[-which(same)]
    if (!1 %in% n.var)
    n.var <- c(n.var, 1)
    return(n.var)
}
#list to character
fixlc<-function(obj){as.character(unlist(obj))}
#factor to numeric
fixlf<-function(obj){as.numeric(as.factor(unlist(obj)))} # should encode text as factors first
#######
### calcluate the correlation modified by Kai
calculate.correlations<-function(data,type="spearman", results = "edge list"){
  check.get.packages(c("impute","WGCNA","Hmisc"))
  #data will be coerced to a matrix
  # type includes pearson (WGCA), biweight(WGCA), spearman
  switch(type,
         pearson 	= .local<-function(data){
           obj<-corAndPvalue(as.matrix(data), use = "pairwise.complete.obs", alternative = "two.sided")
           list(cor=obj$cor,p.value=obj$p)
         },		
         biweight 	= .local<-function(data){
           obj<-bicorAndPvalue(as.matrix(data),use = "pairwise.complete.obs", alternative = "two.sided")
           list(cor=obj$bicor,p.value=obj$p)
         },
         spearman    = .local<-function(data){
           obj<-rcorr(as.matrix(data),type="spearman")
           list(cor=obj$r,p.value=obj$P)
         })
  
  res<-.local(data)	
  
  #add option for matrix or edge list	out put # there is a faster way to get edge list!
  if (results == "edge list"){
    cor<-gen.mat.to.edge.list(res$cor)
    p.vals<-gen.mat.to.edge.list(res$p.value)
    fdr<-p.adjust(p.vals$value,method="BH")
    res<-data.frame(cor,p.values = p.vals[,3],fdr.p.value=fdr)
  }		
  
  return(res)			
}

#test only specific correlations between x and y matrices
# and calculate FDR for only x to y comparisons
xy.correlations<-function(x,y,method="spearman",fdr.method="BH",...){
  #where x and y are data frames
  res<-do.call("rbind",lapply(1:ncol(x),function(i){
    tmp<-do.call("rbind",lapply(1:ncol(y),function(j){
      res<-cor.test(x[,i],y[,j],method=method,...)
      data.frame(source=colnames(x)[i],target=colnames(y)[j],cor=res$estimate,p.value=res$p.value)
    }))
  }))
  res$fdr.p.value<-p.adjust(res$p.value,method=fdr.method)
  return(res)
}	

#accesory function to return position of first instance of unique object 
unique.id<-function(obj)
		{
			tmp<-as.factor(obj)
			id<-seq(along=obj)
			sapply(1:nlevels(tmp),function(i)
				{
					id[tmp==levels(tmp)[i]][1]
				})
		}
#################d
####matrix to edge list ###created by Kai
gen.mat.to.edge.list<-function(mat,symmetric=TRUE,diagonal=FALSE,text=FALSE){
    require(tidyr)
    mat<-as.matrix(mat)
    id<-is.na(mat) # used to allow missing
    mat[id]<-"nna"
    if(symmetric){mat[lower.tri(mat)]<-"na"} # use to allow missing values
    if(!diagonal){diag(mat)<-"na"}
    obj<-as.data.frame(mat)%>%rownames_to_column(var="source")%>%gather(target,value,-source)
    obj<-as.data.frame(obj)
    obj<-obj[!obj$value=="na",]
    obj$value[obj$value=="nna"]<-NA
    if(!text){obj$value<-as.numeric(as.character(obj$value))}
    return(obj)
}
#############################################################################
###construct network
#look up KEGG reactant pairs 
####modified original get.KEGG.pairs to new version ###modified by Kai
get.KEGG.pairs<-function(local=TRUE){
  require(tidyverse)
  if(local){
    DB<-tryCatch(suppressWarnings(read_table2("data/KEGG_reaction_pairs.txt",col_names = F)) ,error=function(e){NULL})
  }else{
    DB<-tryCatch(suppressWarnings(read_table2("https://gist.githubusercontent.com/dgrapov/4964564/raw/121c0924f9e0a8d6b804724a166fc7ee6ea6c38a/KEGG%2520reaction%2520pairs",col_names = F)),error=function(e){NULL})
    }
  return(as.data.frame(DB))
}

#look up CID to KEGG translation this function is replaced with the more general get.Reaction.pairs
####use CTSgetR to get more information compare with get.CID.KEGG.pairs ###modified by Kai
get.CID.KEGG.pairs<-function(cids){
  require(CTSgetR)
  return(CTSgetR(cids,from="PubChem CID",to="KEGG"))
}
#### calculate the reaction relationship ####created by Kai
CID.to.KEGG.pairs<-function(cids){
  require(dplyr)
  database=get.KEGG.pairs()
  lookup=get.CID.KEGG.pairs(cids)
  cc<-as.data.frame(database)%>%dplyr::filter(X1%in%as.vector(lookup[,4]),X2%in%as.vector(lookup[,4]))
  vv<-as.vector(lookup$searchTerm)
  names(vv)<-lookup$value
  cc$X1<-vv[cc$X1]
  cc$X2<-vv[cc$X2]
  colnames(cc)<-c("source","target")
  return(cc)
}

####replace by CID.to.KEGG.pairs ###modified by Kai
#CID.to.KEGG.pairs<-function(cid,database=get.KEGG.pairs(),lookup=get.CID.KEGG.pairs()){
#  
#  matched<-lookup[c(1:nrow(lookup))[fixln(lookup[,1])%in%cid],]
#  ids<-sapply(1:nrow(matched),function(i)
#  {
#    #
#    c(which(as.character(matched[i,2])==as.character(database[,1])),which(as.character(matched[i,2])==as.character(database[,2])))				
#  })
#  
#  names(ids)<-fixln(matched[,1])	# cid of all paired by cid
#  
#  #construct symmetric matrix then extract unique edge list
#  mat<-do.call("rbind",lapply(1:length(ids),function(i)
#  {
#    obj<-ids[[i]]
#    match<-sapply(1:length(ids), function(j)
#    {
#      tmp<-ids[[j]]
#      sum(tmp%in%obj)
#    })
#  }))
#  dimnames(mat)<-list(names(ids),names(ids))
#  elist<-gen.mat.to.edge.list(mat)
 # as.data.frame(elist[fixln(elist[,3])>0,1:2])	#cid source to cid target based on kegg pairs	
#}
##
get.tanimoto.from.SDF<-function(cmpd.DB,type="list",cut.off=0,...){	
  #convert to SDFstr
  #depends on ChemmineR 
  require(ChemmineR)
  cmpd.sdf.list<-new("SDFstr", a = cmpd.DB)
  sd.list<-as(cmpd.sdf.list, "SDFset")
  cid(sd.list) <- sdfid(sd.list)
  # Convert base 64 encoded fingerprints to character vector, matrix or FPset object
  fpset <- fp2bit(sd.list, type=2)
  out<-sapply(rownames(fpset), function(x) ChemmineR::fpSim(x=fpset[x,], fpset,sorted=FALSE)) 
  obj<-as.matrix(out)
  
  if(type=="list"){
    e.list<-gen.mat.to.edge.list(obj)
    final<-edge.list.trim(e.list,index=as.numeric(as.character(unlist(e.list[,3]))),cut=cut.off,less.than=FALSE) # probably overkill
  }else{
    obj[obj<cut.off]<-0
    final<-obj
  }  
  return(final)
}

#get metabolite structure encoding (SDF) from local database or using PubChem webservices
get.SDF.from.CID<-function(cids,DB=NULL,query.limit=25,update.DB=TRUE,save.as="CID.SDF.DB",progress=TRUE,...){
  
  #retrieve metabolite SDF from DB
  #DB should be a list with cids as names
  #for all missing in DB, look up using PubChem PUG
  #if update then update DB with cid entries
  #return list of SDF files for each cid
  
  #remove and print to screen error cids
  #removal ids
  objc<-as.character(unlist(cids))
  objn<-as.numeric(objc)
  dup.id<-duplicated(objn)
  na.id<-is.na(objn)
  if(sum(c(dup.id,na.id))>0){
    
    objc<-as.character(unlist(cids))
    objn<-as.numeric(unlist(cids))
    
    #remove duplicated
    message(cat(paste("The following duplicates were removed:","\n")))
    message(cat(paste(unique(objc[dup.id])),sep="\n"))
    # remove NA
    message(cat(paste("Bad inputs were removed:","\n")))
    message(cat(paste(unique(objc[na.id])),sep="\n"))
    
    cid.objects<-objn[!(na.id | dup.id)]
  } else { 
    cid.objects<-objn 
  }
  
  #check for cid object in local DB
  DB.ids<-names(DB)
  need.id<-!cid.objects%in%DB.ids
  have.id<-DB.ids%in%cid.objects
  cmpd.DB<-DB[have.id]
  
  #use PUG webservices to get missing 
  if(sum(need.id)>0){
    tmp.cids<-unique(cid.objects[need.id])
    
    if(progress){
      message(cat("Using PubChem Power User Gateway (PUG) to get molecular fingerprint(s).\nThis may take a moment...","\n"))
    }
    
    #translate sdf file, modified from ChemmineR
    read.sdf<-function (sdfstr) {
      
      #number of queries controlled in url
      if (length(sdfstr) > 1) {
        mysdf <- sdfstr
      } else {
        mysdf <- readLines(sdfstr)
      }
      
      y <- regexpr("^\\${4,4}", mysdf, perl = TRUE)
      index <- which(y != -1)
      indexDF <- data.frame(start = c(1, index[-length(index)] + 
                                        1), end = index)
      mysdf_list <- lapply(seq(along = indexDF[, 1]), function(x) mysdf[seq(indexDF[x, 
                                                                                    1], indexDF[x, 2])])
      if (class(mysdf_list) != "list") {
        mysdf_list <- list(as.vector(mysdf_list))
      }
      names(mysdf_list) <- 1:length(mysdf_list)
      #mysdf_list <- new("SDFstr", a = mysdf_list)
      return(mysdf_list)
    }
    
    #due to url string size limit query query.limit sdf obj at a time
    blocks<-c(seq(1,length(tmp.cids),by=query.limit),length(tmp.cids))
    # should use cv fold generating fxn to avoid boundary overlaps 
    compounds<-list() #breaks=ceiling(length(cids)/25),include.lowest = TRUE)
    for(i in 1:(length(blocks)-1)){
      url<-paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",paste(tmp.cids[blocks[i]:blocks[(i+1)]],collapse=","),"/SDF")
      compounds[[i]]<-read.sdf(url) 
    }
    
    #create cmpd.list holding sdf strings
    cmpd.list<-list()
    names<-tmp.cids#paste0("CMP",1:length(cid.objects)) # name doesn't matter? should be cid
    for(i in 1:length(compounds)){
      tmp<-unclass(compounds[[i]])
      names(tmp)<-names[blocks[i]:blocks[(i+1)]]
      cmpd.list<-c(cmpd.list,tmp)
    }
    #make sure all are unique
    cmpd.list<-cmpd.list[fixlc(names[!duplicated(names)])] #could get duplicated based on web query
    #combine with DB
    cmpd.DB<-c(cmpd.DB,cmpd.list)
    
    #add new records to DB and save
    if(update.DB) {
      new.DB<-c(DB,cmpd.list)
      assign(save.as,new.DB)
      save(save.as,list=save.as,file=save.as)
    }
  } 
  
  return(cmpd.DB)
  
}

#wrapper to get tanimoto from cid
CID.to.tanimoto<-function(cids,...){
  #wrapper for get sdf from cid 
  #convert sdf to tanimoto similarity
  cmpd.DB<-get.SDF.from.CID(cids,...)
  get.tanimoto.from.SDF(cmpd.DB,...)
}
#querry chemical translation service (CTS) to get tanimoto from inchis (very slow)
CID.to.tanimoto.CTS<-function(cid,lookup=get.CID.INCHIcode.pairs()){
  check.get.packages("XML")
  matched<-lookup[c(1:nrow(lookup))[lookup[,1]%in%cid],]
  matched<-matched[!matched[,2]=="",]
  codes<-gsub("=", "%3D", matched[,2])
  #use webservice to get tanimoto score between cids based in inchi key
  #do in parallel
  check.get.packages(c("snow","doSNOW","foreach"))
  
  i<-1
  out<-list()
  for(i in 1:length(codes))
  {
    cl.tmp = makeCluster(rep("localhost",Sys.getenv('NUMBER_OF_PROCESSORS')), type="SOCK") 
    registerDoSNOW(cl.tmp) 
    .local<-function(i,j,codes)
    {
      url=paste("http://vulcan.fiehnlab.ucdavis.edu:8080/tanimoto-service-1.2/rest/xml/calc.xml?from=",codes[i],"&to=",codes[j],sep="")
      text<-tryCatch(XML::xmlTreeParse(url),error=function(e){NULL})
      if(is.null(text))
      {
        return()
      }else{
        as.numeric(strsplit(unlist(text$doc$children$result[[3]]),">")[3])
      }
    }
    out[[i]]<-foreach(j=c(1:length(codes)),.combine="c") %dopar% .local(i,j,codes=codes)#length(codes)
    stopCluster(cl.tmp)	
  }
  
  names(ids)<-matched[,1]	# cid of all paired by cid
  
  #construct symmetric matrix then extract unique edge list
  mat<-do.call("rbind",lapply(1:length(ids),function(i)
  {
    obj<-ids[[i]]
    match<-sapply(1:length(ids), function(j)
    {
      tmp<-ids[[j]]
      sum(tmp%in%obj)
    })
  }))
  dimnames(mat)<-list(names(ids),names(ids))
  elist<-gen.mat.to.edge.list(mat)
  as.data.frame(elist[elist[,3]==1,1:2])	#cid source to cid target based on kegg pairs	
}
#use chemical resolver to get inchi keys from smiles
get.inchikey.from.smiles<-function(smiles,progress=TRUE){
  # smiles are coerced to a 1 column data frame
  obj<-data.frame(matrix(unlist(smiles),ncol=1))
  if(require(RCurl)==FALSE){install.packages("RCurl");library(RCurl)} else { library(RCurl)} # need RCurl for web querry
  if (progress == TRUE){ pb <- txtProgressBar(min = 0, max = nrow(obj), style = 3)} # show progress bar
  
  start<-"http://cactus.nci.nih.gov/chemical/structure/"
  end<-"/stdinchikey"
  out<-sapply(1:nrow(obj),function(i)
  {
    if (progress == TRUE){setTxtProgressBar(pb, i)}
    
    close(pb)
    url<-paste(start,as.character(unlist(obj[i,])),end,sep="")
    url<-gsub("\\ ","%20",url) # fix spaces 
    tryCatch( getURL(url,ssl.verifypeer=FALSE) ,error=function(e){"error"})
    
  })
  
  if (progress == TRUE){close(pb)}
  
  #format output to only return InchI
  bad<-is.na(smiles)
  out<-as.character(unlist(out))
  out[bad]<-"InChIKey=error"
  #results<-matrix(as.character(unlist(as.data.frame(strsplit(out,"="))[2,])),ncol=1)
  results<-matrix(out,ncol=1)
  colnames(results)<-"InchI Key"
  return(results)
}
get.inchikey.RPAIRS<-function(type="main",url="https://gist.github.com/dgrapov/5674494/raw/9faff56b5f0fe89b554a508bd954605e26b492fc/InchI+Key+Reaction+Pairs"){ 
  #more types should be added based on third column levels
  df<-read.delim("https://gist.githubusercontent.com/dgrapov/5674494/raw/9faff56b5f0fe89b554a508bd954605e26b492fc/InchI+Key+Reaction+Pairs%22",sep="\t",header=F)
  if(require(RCurl)==FALSE){install.packages("RCurl");library(RCurl)} else { library(RCurl)}
  text<-tryCatch( getURL(url,ssl.verifypeer=FALSE) ,error=function(e){NULL})
  tmp<-strsplit(text,"\\n")
  tmp2<-strsplit(as.character(unlist(tmp)), "\t")
  #fix header errors
  tmp2[[1]]<-strsplit(tmp2[[1]],"\\  ")
  full<-out<-matrix(unlist(tmp2),ncol=4, byrow=TRUE)
  
  if(type =="main"){
    out<-full[full[,3]=="main",1:2]
  } 
  
  if(type =="all"){
    out<-full[,1:2]
  }
  
  if(type =="full"){
    out<-full
  }	
  return(out)
}
## transform from one to many other identifiers
## created by Kai Guo
get.meta.information<-function(keys="",from="InChIKey",to=NULL,verbose=TRUE){
    library(tidyverse)
    tmp<-NULL
    keys<-unique(keys)
    keys<-na.omit(keys)
    if(require("CTSgetR",character.only = TRUE)==FALSE)
    {
      library(devtools)
      install_github("dgrapov/CTSgetR")
      library(CTSgetR)
    }else{
        require(CTSgetR)
    }
    if(is.null(to)){
    opt=c("InChIKey","PubChem CID","ChemDB","BioCyc","ChEBI","CAS","Human Metabolome Database","KEGG")
    lef<-setdiff(opt,from)
    }else{
        lef=to
    }
    obj<-do.call(rbind,lapply(1:length(keys),function(id){
        if(isTRUE(verbose)){
            cat("Translation:",id,".......\n")
        }
        t(sapply(1:length(lef),function(i){
        id=keys[id]
        from=from
        to=lef[i]
        tryCatch(
        {
            tmp<-CTSgetR(id,from,to,progress=FALSE)
            colnames(tmp)<-c("from","to")
            tmp
        },
        error=function(e){
            tmp<-data.frame(from=id,to=NA)
            tmp
                    })}))}
    ))
    colnames(obj)<-c(paste0("From:",from),"To")
    obj<-as.data.frame(obj)
    obj$database<-rep(lef,times=length(keys))
    obj<-obj%>%spread(database,To)
    return(obj)
}
