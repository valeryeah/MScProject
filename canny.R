library(imager)


#img = readJPEG("201605201415_dst.jpg") two sandstorms
img = readJPEG("201607180130_dst.jpg")
#v = fit.lr(img)
v = fit.glm(img)
#z = ifelse(v$lr_response>0.5,1,0)
cimg = as.cimg(v$lr_response)
gr = imgradient(cimg,"xy")
plot(gr,layout="row")

blur = isoblur(cimg,5)
plot(blur<0.2)



##STEP 1 DENOISING
img = readJPEG("201706052015_dst.jpg")
v = fit.lr(img)
cimg = as.cimg(v$lr_response)

##STEP2 computing image gradient
gr = imgradient(img.blurred,"xy")
plot(gr,layout="row")

mag <- with(gr,sqrt(x^2+y^2))
plot(mag)

threshold(mag) %>% plot

#Going along the (normalised) gradient
#Xc(im) is an image containing the x coordinates of the image
nX <- Xc(im) + gr$x/mag 
nY <- Yc(im) + gr$y/mag
#nX and nY are not integer values, so we can't use them directly as indices.
#We can use interpolation, though:
val.fwd <- interp(mag,data.frame(x=as.vector(nX),y=as.vector(nY)))

nX <- Xc(im) - gr$x/mag 
nY <- Yc(im) - gr$y/mag
val.bwd <- interp(mag,data.frame(x=as.vector(nX),y=as.vector(nY)))

throw <- (mag < val.bwd) | (mag < val.fwd)
mag[throw] <- 0
plot(mag)

##STEP4 strong and weak edges
t2 <- quantile(mag,.999)

strong <- mag>t2
plot(strong,main="Initial set of strong edges")
plot(mag,main="Initial boundary")



###installing EBImage to do image segmentation
##source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("EBImage")
library(EBImage)
x = bwlabel(mag)
x = as.matrix(x)
max(x)
set1 = which(x==1,arr.ind = TRUE)
set2 = which(x==2,arr.ind = TRUE)
set3 = which(x==3,arr.ind = TRUE)

set1[,1] = 600-set1[,1]
set2[,1] = 600-set2[,1]
set3[,1] = 600-set3[,1]

par(mfrow = c(1,1))
plot(x = set3[,2],y=set3[,1],cex=.2,xlim=c(0,800),ylim=c(0,600)) 

#fitted ellipse
efit = fit.ellipse(x=set1[,2],y=set1[,1])
efit = fit.ellipse(x=set2[,2],y=set2[,1])
efit = fit.ellipse(x=set3[,2],y=set3[,1])
#efit = fit.ellipse(x=set)
efit2 = efit

efit2$angle = pi/2+efit$angle
e <- get.ellipse(efit2)
lines(x=e[,1],y=e[,2],col="red") 









## Naive ellipse fitting
ellip_para2= function(lr_response,p,q,set){
  #center = apply(lr_response[set]*set,2,sum)/sum(lr_response[set])
  center = apply(set,2,mean)
  r=NULL
  for (i in 1:dim(set)[1]) {
    r[i] = dist(rbind(set[i,],center))
  }
  min.r = quantile(r,p)
  maj.r = quantile(r,q)
  a = 0.1
  set2 = which(abs(r-maj.r)<=a)   #to locate the furthest point from center
  if(length(set2)==0) {set2 = which(abs(r-maj.r)<=2*a)}
  #randomly pick a point from this set
  p = set[sample(set2,1),]
  p[1] = 600-p[1] 
  #compute angle of rotation
  theta = atan((p[1]-center[1])/(p[2]-center[2]))
  return(list(r=r,p=p,theta=theta,min.r=min.r,maj.r=maj.r,center = center))
}

v4 = ellip_para2(v$lr_response,0.1,0.9,set2)

##angle significantly affected by the chosen point p

