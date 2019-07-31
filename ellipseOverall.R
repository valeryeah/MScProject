library(imager)
#edge detection
#fit is the fitted logistic regression model
#sigma is the s.d. of blur
#eps is the threshold after the blur
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

#image segmentation, if these are multiple vortices
###installing EBImage to do image segmentation
##source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("EBImage")
library(EBImage)
library(rlist)  #to use list.append

#output l is a list of n boundaries 
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

#using this function to plot!!!
img.seg2 = function(mag,set){
  x = as.matrix(bwlabel(mag))  #generate labels to objects in the binary image
  n = max(x) #number of vortices
  set = which(x==1,arr.ind = TRUE)
  return(list(set=set,n=n,x=x))
}



##fit ellipse to each detected boundary
fitted.ellip = function(l,n){
  set = l[[1]]
  set[,1] = 600-set[,1]
  efit = fit.ellipse(x= set[,2],y=set[,1])
  #efit$angle = pi/2+efit$angle
  ep.para = list(efit)
  if(n>1){
    for (i in 2:n) {
      set = l[[i]]  #the ith element in list
      #set[,1] = 600-set[,1]
      efit = fit.ellipse(x=set[,2],y=set[,1])
      #efit$angle = pi/2+efit$angle
      ep.para = list.append(ep.para,efit)
    }
  }
  return(ep.para = ep.para)
}



get.ellipse <- function(fit, n=360 ) {
  tt <- seq(0, 2*pi, length=n) 
  sa <- sin(fit$angle) 
  ca <- cos(fit$angle) 
  ct <- cos(tt) 
  st <- sin(tt) 
  
  x <- fit$center[1] + fit$major * ct * ca - fit$minor * st * sa 
  y <- fit$center[2] + fit$major * ct * sa + fit$minor * st * ca 
  #y <- fit$center[1] + fit$major * ct * ca - fit$minor * st * sa 
  #x <- fit$center[2] + fit$major * ct * sa + fit$minor * st * ca 
  cbind(x=x, y=y) 
}


img = readJPEG("201706051530_dst.jpg")
###########using this one !
img = readJPEG("201706051100_dst.jpg")

v = fit.glm(img)
w = edge.detection(v)
v5 = image.seg(w$mag,w$set)

l = v5$l


ep.para= fitted.ellip(v5$l,v5$n)

##if want to plot
par(mfrow = c(1,1))
x = v5$x  ##generated labels 
set = which(x==1,arr.ind = TRUE)
set[,1] = 600-set[,1]

#fitted ellipse
e <- get.ellipse(ep.para[[1]])

plot(x = set[,2],y=set[,1],cex=.2,xlim=c(0,800),ylim=c(0,600),xlab="",ylab="",main="Fit Ellipse") 
#image = load.image("201706051100_dst.jpg") #using the imager package
#plot(image)

#lines(x = set[,2],y=set[,1],cex=.2,xlim=c(0,800),ylim=c(0,600),lty=2)
lines(x=e[,1],y=e[,2],col="red")
e2 <- get.ellipse(ep.para[[2]])
#add second vortex
set2 = which(x==2,arr.ind = TRUE)
set2[,1] = 600-set2[,1]
lines(x = set2[,2],y=set2[,1],cex=.2,xlim=c(0,800),ylim=c(0,600),lwd="0.5",lty=2) 
lines(x=e2[,1],y=600-e2[,2],col="blue",lwd="0.5")
legend("bottomleft", legend=c("Detected boundary","Ellipse 1","Ellipse 2"), pch=c(NA,NA), col=c("black","red","blue"), lty=c(1,1,2),cex=0.75,lwd=c(2,1,1))



cscale = function(x) ifelse(cimg[x,x,1,1]==0,rgb(1,1,1),rgb(x,0,x))

ycimg = as.cimg(y)
ycimg= imrotate(ycimg,90)
ycimg=mirror(ycimg,"x")
cscale = function(x) rgb(1,1-x,1)
plot(ycimg,colourscale=cscale,rescale=FALSE,main="Binary Sandstorm")

#y.lr = fit.glm(img)$lr_response
ycimg2 = as.cimg(y.lr)
ycimg2= imrotate(ycimg2,90)
ycimg2=mirror(ycimg2,"x")
plot(ycimg2,colourscale=cscale,rescale=FALSE,main="LR Smoothing")

blur = isoblur(ycimg2,4)  #blur image
blur= imrotate(blur,90)
blur=mirror(blur,"x")
plot(blur,colourscale=cscale,rescale=FALSE,main="Blurring")

#selecting threshold
img.blurred = blur>0.2
plot(img.blurred,colourscale=cscale,rescale=FALSE,main="Final Sandstorm")



gau = function(x,sigma){
  exp(-x^2/(2*sigma^2))/sqrt(2*pi*sigma^2)
}






