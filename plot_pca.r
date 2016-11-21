pca_plot<-function(dat,group){
library(ggplot2);
pca<-prcomp(t(dat));
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
d<-data.frame(PC1=pca$x[,1],PC2=pca$x[,2],group=group,name=colnames(dat));
ggplot(data=d,aes(x=PC1,y=PC2,color=group))+geom_point(size=3)+labs(title="PCA")+
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  geom_text(aes(label=name),vjust=min(pca$x[,2])/20)
}
