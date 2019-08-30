library(jpeg) #to load image


#generate binary labels
#input is image
#output it labels indicating magenta
bin_class = function(img){
  R = img[,,1]
  G = img[,,2]
  B = img[,,3]
  y = matrix(0,dim(img)[1],dim(img)[2])
  for (i in 1:dim(img)[1]) {
    for (j in 1:dim(img)[2]) {
      #difference to magenta within 15%, classify as 1
      if (abs(R[i,j]-1)<=0.15 &&  abs(G[i,j]-0)<=0.5 && abs(B[i,j]-1)<=0.15)
        y[i,j] = 1
    }
  }
  return(y)
}



##applying logistic regression
#output is a list with model and probabilities
fit.glm = function(img){
  R = img[,,1]
  G = img[,,2]
  B = img[,,3]
  R = as.vector(R)
  G = as.vector(G)
  B = as.vector(B)
  y = as.vector(bin_class(img))
  df = as.data.frame(cbind(R,G,B,y))
  fit=glm(y~R+G+B,data=df,family=binomial())
  lr_response = matrix(data=predict(fit,type="response"),nrow=dim(img)[1],ncol=dim(img)[2],byrow=FALSE)
  return(list(fit=fit,lr_response = lr_response))
}


library(imager)
### imaging blurring and edge detection
#fit is the fitted logistic regression model
#sigma is the s.d. of blur
#eps is the threshold to select after blurring
edge.detection = function(fit,sigma=5,eps=0.03){
  cimg = as.cimg(fit$lr_response)   #convert to cimg 
  blur = isoblur(cimg,sigma)  #blur image
  img.blurred = blur>eps #set threshold
  gr = imgradient(img.blurred,"xy")   #compute gradient after the blur
  mag <- with(gr,sqrt(x^2+y^2))
  
  #boundary is where mag nonzero
  set = which(mag!=0,arr.ind = TRUE)
  set = as.matrix(set)
  return(list(mag=mag,set=set))
}


library(EBImage)
library(rlist)  #need to use list.append
#image segmentation, if these are multiple sandstorms in one image
###installing EBImage to do image segmentation
##source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("EBImage")

#output l is a list of coordinates of n objects and total number of objects
image.seg = function(mag,set){
  x = as.matrix(bwlabel(mag))  #generate labels to objects in the binary image
  n = max(x) #number of vortices
  set = which(x==1,arr.ind = TRUE)
  l = list(set)
  if(n>1){
    for (i in 2:n) {
    set = which(x==i,arr.ind = TRUE)
    l = list.append(l,set)
    }
  }
  return(list(l=l,n=n,x=x))
}


