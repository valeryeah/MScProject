library(magick)
library(wvtool)
library(raster)

#img = readJPEG("201805280730_dst.jpg")
#img.raster = raster("201605201415_dst.jpg")
#plot(img.raster)
#v = fit.lr(img)
#v = fit.glm(img)
#v2=find_bound(v$lr_response)
#raster = as.raster(v2$bd)
#plot(raster)

v = fit.glm(img)

cimg = as.cimg(v$lr_response)
#gr = imgradient(cimg,"xy")
#plot(gr,layout="row")

blur = isoblur(cimg,5)  #prob still greater than .5 after blurred

#library(spatstat)  #to convert image into matrix
img.blurred = blur>0.05
m = blur[,,1,1]
image(1:800,1:600,t(apply(m,2, rev)),col=c("white","magenta"),main="Blurred image")

gr = imgradient(img.blurred,"xy")

mag <- with(gr,sqrt(x^2+y^2))
m2 = mag[,,1,1]
image(1:800,1:600,t(apply(m2,2, rev)),col=c("white","black"),main="Detected Boundary")

plot(mag)
set = which(mag!=0,arr.ind = TRUE)
set = as.matrix(set)


#center2 = ellip_center2(v$lr_response)



#set = which(v2$bd==1,arr.ind = TRUE)
set[,1] = 600-set[,1]

par(mfrow = c(1,1))
#plot(x = set[,2],y=set[,1],cex=.2,xlim=c(0,800),ylim=c(600,0)) 
plot(x = set[,2],y=set[,1],cex=.2,xlim=c(0,800),ylim=c(0,600)) 

#fitted ellipse
efit = fit.ellipse(x=set[,2],y=set[,1])
#efit = fit.ellipse(x=set)
#efit2 = efit
#because flipped y axis, mirror image of the original ellipse
#efit2$angle = -efit$angle
e <- get.ellipse(efit)

lines(x=e[,1],y=e[,2],col="red") 

#print(efit2)


#ellip_center1(v$lr_response)

library(imager)
#visualize detected boundary

image = load.image("201605201415_dst.jpg")
plot(image)
mag.raster = as.raster(mag)
img.raster = as.raster(img)
overlay(mag.raster,img.rast)

