ggcorr.test<-function(data,method ="pearson",use="complete.obs",cor_matrix = NULL,
                      nbreaks = NULL, digits = 2, name = "", low = "#3B9AB2", mid = "#EEEEEE", 
                      high = "#F21A00", midpoint = 0, palette = NULL, 
                      label = FALSE,label.sep="\n",pvalue=FALSE,star=TRUE,label_color = "black", 
                      label_size = 2.5,legend.position = "right", legend.size = 9,upper=FALSE,xangle=0,...)
{
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      data = as.data.frame(data)
    }
    x = which(!sapply(data, is.numeric))
    if (length(x) > 0) {
      warning(paste("data in column(s)", paste0(paste0("'", names(data)[x], "'"), collapse = ", "), "are not numeric and were ignored"))
      data = data[, -x]
    }
  }
  name<-combn(colnames(data),2,simplify=F)
  if (is.null(cor_matrix)) {
    cor_name=do.call(rbind,name)
    colnames(cor_name)<-c("V1","V2")
    cor_matrix=cbind(cor_name,do.call(rbind,lapply(name,function(x).cor_test(data[,x[1]],data[,x[2]],method=method,use=use))))
    if(upper==TRUE){
      cor_matrix$V1<-factor(cor_matrix$V1,levels=unique(cor_matrix$V1))
      cor_matrix$V2<-factor(cor_matrix$V2,levels=unique(cor_matrix$V2))
    }else{
      cor_matrix$V1<-factor(cor_matrix$V1,levels=rev(unique(cor_matrix$V1)))
      cor_matrix$V2<-factor(cor_matrix$V2,levels=rev(unique(cor_matrix$V2)))
    }
  }
  tmp<-NULL
  for(i in 1:nrow(cor_matrix)){
    if(cor_matrix$p[i]<0.01){
      tmp[i]<-"**"
    }
    if(cor_matrix$p[i]<0.05&cor_matrix$p[i]>=0.01){
      tmp[i]<-"*"
    }
    if(cor_matrix$p[i]>=0.05){
      tmp[i]<-" "
    }
  }
  if(star==TRUE){
    cor_matrix$sig<-tmp
  }else{
    cor_matrix$sig<-""
  }
  p<-ggplot(cor_matrix,aes(V1,V2))
  p<-p+geom_tile(aes(fill=estimate),color="white")+
    scale_fill_gradient2(low=low,mid=mid,high=high,midpoint=0,limit = c(-1,1), name=paste(.simpleCap(method),"Correlation",sep="\n"))+theme_minimal()
  if(label==TRUE){
    if(pvalue==TRUE){
      p<-p+geom_text(aes(label=paste(r,paste0("(",p,sig,")"),sep=label.sep)),size=label_size)
    }else{
      p<-p+geom_text(aes(label=r),size=label_size)
    }
  }
  p<-p+theme(axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank())+xlab("")+ylab("")
  if(upper==TRUE){
    p<-p+scale_x_discrete(position = "top")+theme(axis.text.x=element_text(angle=xangle,vjust=-0.1,hjust=0))
  }else{
    p<-p+scale_y_discrete(position = "right")+theme(axis.text.x=element_text(angle=xangle,vjust=1,hjust=1))
  }
  p
}
.cor_test<-function(x, y, method = "pearson",use = "complete.obs", label.sep = ", ", output.type = "expression"){
  require(tidyverse)
  .cor <- stats::cor.test(x, y, method = method, exact = FALSE, use = use)
  estimate <- p.value <- p <- r <- rr <-  NULL
  z <- data.frame(estimate = .cor$estimate, p.value = .cor$p.value, method = method) %>%
    mutate(
      r = signif(estimate, 2),
      rr = signif(estimate^2, 2),
      p = signif(p.value, 2)
    )
  z
}


.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

