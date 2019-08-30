library(mvtnorm)
library(mclust)  #gmm


############################### fitting mixture model ###################################
#n is the number of components
gmm.fit = function(img,n=2){
  fit = fit.glm(img)
  cimg = as.cimg(fit$lr_response)   #convert to cimg 
  blur = isoblur(cimg,sigma=5)  #blur image
  eps=0.03
  img.blurred = blur>eps
  lab = bwlabel(img.blurred)
  x1 = which(lab=="1",arr.ind = TRUE)[,1:2]  ##all points inside the sandstorm
  
  gr = imgradient(img.blurred,"xy")   
  mag <- with(gr,sqrt(x^2+y^2))
  set = which(mag!=0,arr.ind = TRUE)
  set = as.matrix(set)[,1:2]  ##sandstorm boundary
  
  x1.plot = x1
  x1.plot[,1] = 600-x1.plot[,1]
  x1[,1] = x1.plot[,2]
  x1[,2] = x1.plot[,1]
  
  set.plot = set
  set.plot[,1] = 600-set[,1]
  set[,1] = set.plot[,2]
  set[,2] = set.plot[,1]
  
  #setting prior, n-component mixture
  defaultPrior(data= x1,G=2,"EEE")
  gmm.fit = densityMclust(data=x1,G=2,prior=priorControl(functionName = "defaultPrior"),modelNames="EEE")
  return(list(gmm.fit=gmm.fit,set=set))
}










