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
