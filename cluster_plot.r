#################################################################################
# Author:Kai Guo guokai8@gmail.com ##############################################
#
# Last modified: 2016-11-23 09:14######################################
#
# Filename: cluster_plot.r
#
# Description: This script is writed for ....####
##################################################################################
#!/usr/bin/R
############################check all below packages installed or not
packages <- c("ggplot2", "dplyr","devtools","factoextra")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())),repos="http://cran.r-project.org")  
}
###########################plot 2D-PCA
pca_plot<-function(dat,group,filename="./pca.pdf"){ ###########dat: the dataset you want do PCA analysis. group: experiment design.
library(ggplot2);
pca<-prcomp(t(dat));
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
d<-data.frame(PC1=pca$x[,1],PC2=pca$x[,2],group=group,name=colnames(dat));
p<-ggplot(data=d,aes(x=PC1,y=PC2,color=group))+geom_point(size=3)+labs(title="PCA")+
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  geom_text(aes(label=name),vjust=min(pca$x[,2])/20)
ggsave(p,file=filename);
##print(p) remove the comment if you want print it on the screen
}
#############################plot kmeans cluster
km_plot<-function(df,group){
library(factoextra)
########## df:the dataset you want to do kmeans analysis. group: experiment design
km<-eclust(t(df),"kmeans",k=length(unique(group)),nstart=25,graph=FALSE)
p<-fviz_cluster(km,frame.type = "norm", frame.level = 0.68,labelsize=4)
ggsave(p,file="./kmeans.pdf") ###remove this if you want print it on the screen
}
############################hclust plot(simple version)
plot_hc<-function(df,group) {
  #################df:the dataset you want to do kmeans analysis. group: experiment design
  cc<-cor(df)
  dd<-as.dist(1-cc)
  hc<-hclust(dd)
  pdf(file="./hclust.pdf") ####remove this if you want print it on the screen
  plot(hc,xlab="Hierarchical cluster",sub="",cex=0.5)
  dev.off()
}
###########################correlation heatmap 
library(ggplot2)
library(reshape2)
#Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA; 
  return(cormat);
} 
# Get upper triangle of the correlation matrix 
get_upper_tri <- function(cormat){ 
  cormat[lower.tri(cormat)]<- NA; 
  return(cormat); 
};
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2);
  hc <- hclust(dd);
  cormat <-cormat[hc$order, hc$order];
}
gg_heatmap<-function(df){ #dataset input
  cormat<-round(cor(df),2);
  cormat <- reorder_cormat(cormat);
  upper_tri <- get_upper_tri(cormat);
  melted_cormat <- melt(upper_tri);
  melted_cormat <- na.omit(melted_cormat);
  p <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+ 
    geom_tile(color = "white")+ 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), name="Pearson\nCorrelation")+ 
      theme_minimal()+theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+ coord_fixed()
  p<-p+geom_text(aes(Var2, Var1, label = value), color = "white", size = 2) + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          legend.justification = c(1, 0), 
          legend.position = c(0.6, 0.7), 
          legend.direction = "horizontal")+
          guides(fill = guide_colorbar(barwidth = 7,barheight = 1, title.position = "top", title.hjust = 0.5))
  ggsave(p,file="./ggheatmap.pdf") ###replace this if you want print it on the screen
}

#####################plot 3D pca !!!important: you may like to modify the figure by using mouse#######
if (!require(car,character.only=TRUE)){
    install.packages(pkgs="car",repos="http://cran.r-project.org")
}
if (!require(rgl,character.only=TRUE)){
    install.packages(pkgs="rgl",repos="http://cran.r-project.org")
}
pca3d<-function(df,group){ #df:the dataset you want to do kmeans analysis. group: experiment design
  library(car)
  library(rgl)
  pca<-prcomp(t(df))
  pc1<-as.vector(pca$x[,1])
  pc2<-as.vector(pca$x[,2])
  pc3<-as.vector(pca$x[,3])
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  scatter3d(pc1,pc2,pc3,group=factor(group),
            grid=F,surface=F,ellipsoid = TRUE,
            labels = colnames(df), id.n=ncol(df),
            xlab=paste0("PC1: ",round(percentVar[1] * 100),"% variance"),
            ylab=paste0("PC2: ",round(percentVar[2] * 100),"% variance"),
            zlab=paste0("PC3: ",round(percentVar[3] * 100),"% variance"),
            axis.scales = FALSE,surface.col = c("blue", "red"),ellipsoid.alpha=0.1,smooth=T)
  rgl.snapshot(filename = "./3dpca.png")
}

