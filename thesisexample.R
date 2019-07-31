library(ggplot2)

sandstorm1 = readJPEG("201706052115_dst.jpg") #one sandstorm
thesis = readJPEG("201706051300_dst.jpg")  ##two sandstorm

fit.thesis = fit.glm(thesis)
e.thesis = edge.detection(fit.thesis,sigma=5,eps = 0.03)
l = image.seg(e.thesis$mag,e.thesis$set)$l
set1 = l[[1]]
set2 = l[[2]]
set1[,1] = 600-set1[,1]
set2[,1] = 600-set2[,1]

d = data.frame(x=c(set1[,2],set2[,2]),y = c(set1[,1],set2[,1]),names=c(rep("1",1027),rep("2",240)))
gg = ggplot(d,aes(x,y))+geom_point(aes(color = names))+geom_text(aes(label=names),size = 1)+xlim(0,800)+ylim(0,600)+ggtitle("Segmented Image")
gg+theme(
  plot.title = element_text(size=15,hjust=0.5),
  axis.title.x = element_blank(),
  axis.title.y = element_blank())


plot(x= set[,2],y=set[,1],xlim = c(0,800),ylim=c(0,600),cex=0.05,col="magenta",main="Detected Boundary",xlab="",ylab="")

########################  blurring effect  ##########################

cscale = function(x) ifelse(cimg[x,x,1,1]==0,rgb(1,1,1),rgb(x,0,x))
img=readJPEG("201706060045_dst.jpg")
v = fit.glm(img)
y.lr = v$lr_response
ycimg2 = as.cimg(y.lr)
plot(ycimg2,colourscale=cscale,rescale=FALSE,main="LR Smoothing")

ycimg2= imrotate(ycimg2,90)
ycimg2=mirror(ycimg2,"x")
plot(ycimg2,colourscale=cscale,rescale=FALSE,main="LR Smoothing")

#sd = 1,4,7
blur = isoblur(ycimg2,4)  #blur image
#blur= imrotate(blur,90)
#blur=mirror(blur,"x")
plot(blur,colourscale=cscale,rescale=FALSE,main="Blurring with sd=4")

#eps=0.1
img.blurred = blur>0.007
plot(img.blurred,colourscale=cscale,rescale=FALSE,main="eps=0.01")

img.blurred = blur>0.1
plot(img.blurred,colourscale=cscale,rescale=FALSE,main="eps = 0.1")


#selecting threshold
img.blurred = blur>0.5
plot(img.blurred,colourscale=cscale,rescale=FALSE,main="eps = 0.5")




########################  fitting circles   ##########################


img = readJPEG("201706052115_dst.jpg")
img = readJPEG("201805281015_dst.jpg")
v = fit.glm(img)
w = edge.detection(v)
v5 = image.seg(w$mag,w$set)
l = v5$l
set = l[[1]]
set2 = set
set[,2] = 600-set2[,1]
set[,1] = set2[,2]

cent = circle_cent(v$lr_response)
cent2 = cent
cent[2] = 600-cent2[1]
cent[1] = cent2[2]
circle_rad = function(cent,set){
  d = rep(0,nrow(set))
  for (i in 1:nrow(set)) {
    d[i] = dist(rbind(cent,set[i,]))
  }
  r = max(d)
  return(r)
}
radius = circle_rad(cent,set)



plot(x = set[,1],y=set[,2],cex=.1,xlim=c(0,800),ylim=c(0,600),xlab="",ylab="",main="Fitting Circle") 
theta = seq(0, 2*pi, length = 200)
lines(x = cent[1]+radius * cos(theta), y = cent[2]+ radius * sin(theta),col="red",lwd=1.5,lty=2)
points(x = cent[1],y=cent[2],pch=20,cex=0.5,col="red")
legend("topleft", legend=c("fitted circle", "centre","boundary"), col=c("red", "red","black"), lty=c(2,NA,1), pch=c(NA,20,NA),lwd=c(1,NA,2),cex=0.8)



########################  fitting ellipses ##########################
img = readJPEG("201607181145_dst.jpg")
img = readJPEG("201503201415_dst.jpg")
img = readJPEG("201607171245_dst.jpg")
img = readJPEG("201605201415_dst.jpg")


img = readJPEG("201606291515_dst.jpg")
v = fit.glm(img)
w = edge.detection(v)
v5 = image.seg(w$mag,w$set)
v5$n
x = v5$x  ##generated labels 
set = which(x==1,arr.ind = TRUE)
set[,1] = 600-set[,1]

efit = fit.ellipse(x= set[,2],y=set[,1])
efit$angle = efit$angle + pi/2
e <- get.ellipse(efit)

plot(x = set[,2],y=set[,1],cex=.05,xlim=c(350,650),ylim=c(150,550),xlab="",ylab="",main="Multiple Sandstorms") 
lines(x=e[,1],y=e[,2],col="red",lwd=2,lty=2)


set2 = which(x==2,arr.ind = TRUE)
set2[,1] = 600-set2[,1]
points(x = set2[,2],y=set2[,1],cex=.05)

efit = fit.ellipse(x= set2[,2],y=set2[,1])
efit$angle = efit$angle + pi/2
e <- get.ellipse(efit)
lines(x=e[,1],y=e[,2],col="blue",lwd=2,lty=2)


set3 = which(x==3,arr.ind = TRUE)
set3[,1] = 600-set3[,1]
points(x = set3[,2],y=set3[,1],cex=.05)

efit = fit.ellipse(x= set3[,2],y=set3[,1])
efit$angle = efit$angle + pi/2
e <- get.ellipse(efit)
lines(x=e[,1],y=e[,2],col="green",lwd=2,lty=2)


legend("topright", legend=c("sandstorm1", "sandstorm2","sandstorm3","boundary"), col=c("red", "blue","green","black"), lty=c(2,2,2,1),lwd=c(2,2,2,1),cex=0.75)










########################  fitting GMM ##########################

sandstorm1 = readJPEG("201706052115_dst.jpg") #one sandstorm
fit = fit.glm(sandstorm1)
cimg = as.cimg(fit$lr_response)   #convert to cimg 
blur = isoblur(cimg,sigma=5)  #blur image
eps=0.03
img.blurred = blur>eps
lab = bwlabel(img.blurred)
x1 = which(lab=="1",arr.ind = TRUE)[,1:2]  ##all points inside

gr = imgradient(img.blurred,"xy")   
mag <- with(gr,sqrt(x^2+y^2))
set = which(mag!=0,arr.ind = TRUE)
set = as.matrix(set)[,1:2]  


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
gmm.fit4e = densityMclust(data=x1,G=4,prior=priorControl(functionName = "defaultPrior"),modelNames="EEE")


##plot the GMM with boundary 

plot(gmm.fit4e,what="density",xlim=c(400,670),ylim=c(170,400),main=TRUE)
points(set,cex=0.1,col="blue")

#plotting with ellipse
efit = fit.ellipse(x= set[,1],y=set[,2])
#efit$angle = efit$angle + pi/2
e <- get.ellipse(efit)

lines(x=e[,1],y=e[,2],col="red",lwd=2,lty=2)




#plot the specific contour
densty = dens(modelName = gmm.fit4e$modelName, data = set, parameters = gmm.fit4e$parameters)
m4e=mean(densty)
plotDensityMclust2(gmm.fit4e,data=set,nlevels=1,levels=m4e,points.cex = 0.1,col="red",lwd=2)

##model type 
gmm.fit$modelName  #VVV

##when model name ="EEE", equal covariance
gmm.fit$parameters$variance


############################### GMM MODEL COMPARISON #################################
##2 COMPONENT, VVV
plot(gmm.fit2,what="density",xlim=c(460,610),ylim=c(230,345),main=TRUE)
points(set,cex=0.1,col="blue")
plotDensityMclust2(gmm.fit2,data=set,nlevels=1,levels=m2,points.cex = 0.1,col="red",lwd=2)

##2 COMPONENT, EEE
plot(gmm.fit2e,what="density",xlim=c(450,610),ylim=c(230,345),main=TRUE)
points(set,cex=0.1,col="blue")
plotDensityMclust2(gmm.fit2e,data=set,nlevels=1,levels=m2e,points.cex = 0.1,col="red",lwd=2)

##4 COMPONENT, VVV
plot(gmm.fit4,what="density",xlim=c(460,610),ylim=c(230,345),main=TRUE)
points(set,cex=0.1,col="blue")
plotDensityMclust2(gmm.fit4,data=set,nlevels=1,levels=m4,points.cex = 0.1,col="red",lwd=2)

##4 COMPONENT, EEE
plot(gmm.fit4e,what="density",xlim=c(450,615),ylim=c(230,345),main=TRUE)
points(set,cex=0.1,col="blue")
plotDensityMclust2(gmm.fit4e,data=set,nlevels=1,levels=m4e,points.cex = 0.1,col="red",lwd=2,main=TRUE)



###############################  MH sampling for bounadry ################################
## predicting 10*15min ahead
set.seed(1234)

# h=2, WITH CONSTRAINT
plot(x.all2t,cex=0.2,main="h=2,CONSTRAINT=T",ylim=c(460,600))
points(true2,cex=0.2,col="blue")

# h=2 WITHOUT CONSTRAINT
plot(x.all2f,cex=0.2,main="h=2,CONSTRAINT=F")
points(true2,cex=0.2,col="blue")

# h=10, WITH CONSTRAINT
plot(x.all10t,cex=0.2,main="h=10,CONSTRAINT=T",ylim=c(460,600))
points(true10,cex=0.2,col="blue")

# h=10, WITHOUT CONSTRAINT
plot(x.all10f,cex=0.2,main="h=10,CONSTRAINT=F")
points(true10,cex=0.2,col="blue")


# compute acceptance rate
n=1.5e4
#with constraint
x = x.all2t[,1] #acceptance rate
mean(x[-1]!=x[-(n+1)]) #0.4823333
#without constraint
x = x.all2f[,1]
mean(x[-1]!=x[-(n+1)]) #0.6629  


############################### UNIVARIATE TS forecasting ################################
## UNIVARIATE forecasting major radius
majr.ts = arima(ep.train[,3],order= c(1,1,1)) 
minr.ts = arima(ep.train[,4],order= c(0,1,1))
auto.arima(ep.train[,3])

forecast(majr.ts,3)$mean
#true majr 
ep.test[1:3,3]

## UNI forecasting mean and sd of bd points in x,y directions 
m.x = arima(bd.train[,1],order= c(2,2,2))
forecast(m.x,5)$mean
bd.test[1:5,1]
auto.arima(bd.train[,1])  #(0,1,2)

m.y = arima(bd.train[,2],order= c(2,2,2))
forecast(m.y,5)$mean
bd.test[1:5,2]
auto.arima(bd.train[,2])   #(0,1,0)

s.d.x = arima(bd.train[,3],order= c(1,1,1))
forecast(s.d.x,5)$mean
bd.test[1:5,3]
auto.arima(bd.train[,3]) #(1,1,0)

s.d.y = arima(bd.train[,4],order= c(1,1,1))
forecast(s.d.y,5)$mean
bd.test[1:5,4]
auto.arima(bd.train[,4]) #(0,1,0)


############################### MULTIVARIATE TS forecasting ################################






