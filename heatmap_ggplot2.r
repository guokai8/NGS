```
df <- mtcars[, c(1,3,4,5,6,7)]
cormat <- round(cor(df),2)
library(reshape2);
mc <- melt(cormat);
library(ggplot2)
ggplot(data = mc, aes(x=Var1, y=Var2, fill=value)) + geom_tile()
##
get_lower_tri<-function(cormat){ 
cormat[upper.tri(cormat)] <- NA
return(cormat)
} # Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
cormat[lower.tri(cormat)]<- NA 
return(cormat)
}
upper_tri <- get_upper_tri(cormat)
mc <- melt(upper_tri)
mc<-na.omit(mc)
ggplot(data = mc, aes(Var2, Var1, fill = value))+
geom_tile(color = "white")+
scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), name="Pearson\nCorrelation") +
theme_minimal()+
theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+ 
coord_fixed()
```
