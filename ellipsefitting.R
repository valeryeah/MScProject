#load image
library(jpeg)
img1=readJPEG("201706241900_dst.jpg")
img=readJPEG("201706241615_dst.jpg")
img = readJPEG("201706052115_dst.jpg")
#magenta = (1,0,1)
#binary classification


#generate labels
bin_class = function(img){
  R = img[,,1]
  G = img[,,2]
  B = img[,,3]
  y = matrix(0,dim(img)[1],dim(img)[2])
  for (i in 1:dim(img)[1]) {
    for (j in 1:dim(img)[2]) {
      #if difference to (1,0,1) i.e. magenta within 10%, classify as 1
      if (abs(R[i,j]-1)<=0.15 &&  abs(G[i,j]-0)<=0.5 && abs(B[i,j]-1)<=0.15)
        y[i,j] = 1
    }
  }
  return(y)
}
#y = bin_class(img)
#image(1:800,1:600,t(apply(y,2, rev)),col=c("white","magenta"),main="Detect dust storm")
yr = as.raster(y)
yr[which(y=="1",arr.ind = TRUE)] <-"magenta"
yr[which(y=="0",arr.ind = TRUE)] <-"white"


plot(yr,main="binary image")



#convert matrices into vectors col by col
y = as.vector(bin_class(img))
R = as.vector(R)
G = as.vector(G)
B = as.vector(B)
df = as.data.frame(cbind(R,G,B,y))
#build logistic regression model
fit=glm(y~R+G+B,data=df,family=binomial())
#fitted probs 0 or 1



library(logistf)
# "separation" problem of fitting glm
# fit fith's logistic regression 
#return fitted probabilities in a matrix form, with same dimension as image
flr = function(img){
  y = as.vector(bin_class(img))
  R = as.vector(img[,,1])
  G = as.vector(img[,,2])
  B = as.vector(img[,,3])
  df = as.data.frame(cbind(R,G,B,y))
  fit = logistf(y~R+G+B,data=df)
  probs = matrix(data=fit$predict,nrow=dim(img)[1],ncol=dim(img)[2],byrow = FALSE)
  return(probs)
}


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

#fit firth's logistic regression
library(logistf)
fit.lr = function(img){
  R = img[,,1]
  G = img[,,2]
  B = img[,,3]
  R = as.vector(R)
  G = as.vector(G)
  B = as.vector(B)
  y = as.vector(bin_class(img))
  df = as.data.frame(cbind(R,G,B,y))
  fit = logistf(y~R+G+B,data=df)
  lr_response = matrix(data=fit$predict,nrow=dim(img)[1],ncol=dim(img)[2],byrow=FALSE)
  return(list(fit=fit,lr_response = lr_response))
}

img = readJPEG("201706051100_dst.jpg")
fit = fit.glm(img)
summary(fit$fit)
fit2 = fit.lr(img)
summary(fit2$fit)

fit3 = logistf(y~R+G+B,data=df)
#reponse vector:
probs = fit3$predict
#convert into matrix
lr_response = matrix(data=probs,nrow=dim(img2)[1],ncol=dim(img2)[2],byrow=FALSE)

#for pixel[i,j], p(y=1) = prob[600*(j-1)+i]


#############################    fitting circle ########################################
circle_center1 = function(lr_response){
  center=c(0,0)
  for (i in 1:dim(lr_response)[1]) {
    for (j in 1:dim(lr_response)[2]) {
      center= center+lr_response[i,j]*c(i,j)
    }
  }
  center=center/sum(lr_response)
  return(center)
}

row = 1:600
col = 1:800

circle_center1(fit$lr_response)

img = readJPEG("201706052115_dst.jpg")
v = fit.glm(img)
w = edge.detection(v)
v5 = image.seg(w$mag,w$set)
l = v5$l
set = l[[1]]
set2 = set
set[,2] = 600-set2[,1]
set[,1] = set2[,2]

cent = circle_center1(fit$lr_response)
cent2 = cent
cent[2] = 600-cent2[1]
cent[1] = cent2[2]
circle_rad = function(cent,set){
  d = rep(0,600)
  for (i in 1:600) {
    d[i] = dist(rbind(cent,set[i,]))
  }
  r = max(d)
  return(r)
}
radius = circle_rad(cent,set)



plot(x = set[,1],y=set[,2],cex=.2,xlim=c(0,800),ylim=c(0,600),xlab="",ylab="",main="Fit Ellipse") 
pp = 2*pi
theta = seq(0, pp, length = 200)
lines(x = radius * cos(theta), y = radius * sin(theta))







#image function turns it upside down
#z_thre = bin_class(img2)
#image(1:800,1:600,t(apply(z_thre,2, rev)),col=c("white","magenta"))

#
plot_fitted = function(lr_response){
  z_lr = ifelse(lr_response>0.5,1,0)
  image(1:800,1:600,t(apply(z_lr,2, rev)),col=c("white","magenta"))
}


##knowing center c find minor and major axes of ellipse
#z_lr = ifelse(lr_response>0.5,1,0)

##finding boundaries by converting interior points into 0
find_bound = function(lr_response){
  z_lr = ifelse(lr_response>0.5,1,0)
  bd = z_lr
  for (i in 2:(dim(bd)[1]-1)) {
    for (j in 2:(dim(bd)[2]-1)) {
      #point lying inside the ellipse
      if(z_lr[i-1,j]==1 && z_lr[i+1,j]==1 && z_lr[i,j-1]==1 && z_lr[i,j+1]==1){
        bd[i,j]=0
      }
      #point lying outside the ellipse
      if(z_lr[i-1,j]==0 && z_lr[i+1,j]==0 && z_lr[i,j-1]==0 && z_lr[i,j+1]==0 &&z_lr[i-1,j-1]==0 && z_lr[i-1,j+1]==0 && z_lr[i+1,j-1]==0 && z_lr[i+1,j+1]==0){
        bd[i,j]=0
      }
    }
  }
  return(list(bd=bd,z_lr=z_lr))
}

#plotting boundary
#bd = find_bound(lr_response)$bd #600*800 matrix 
#image(1:800,1:600,t(apply(bd,2, rev)),col=c("white","magenta"))

#y = function(x) x+ 310-529
#t = seq(528,530,by=1)
#lines(t,y(t))


#######wrong 
ellip_center2 = function(lr_response){
  bd = find_bound(lr_response)$bd
  s = which(bd==1,arr.ind=TRUE)
  center=c(0,0)
  d = 0
  for (i in 1:dim(s)[1]) {
    center = center + s[i,] * lr_response[s[i,][1],s[i,][2]]
    d = d + lr_response[s[i,][1],s[i,][2]]
  }
  center = center/d
  return(center)
}
#center = ellip_center(lr_response)



#z_lr is the binary classification of lr, c is the center
#z_lr = find_bound(lr_response)$z_lr

##minor and major radius are the pth and qth quantile of distances to center
ellip_para= function(lr_response,center,p,q){
  z_lr = find_bound(lr_response)$z_lr
  bd = find_bound(lr_response)$bd
  s = which(bd==1,arr.ind = TRUE)
  r=NULL
  for (i in 1:dim(s)[1]) {
    r[i] = dist(rbind(s[i,],center))
  }
  min.r = quantile(r,p)
  maj.r = quantile(r,q)
  a = 0.1
  set = which(abs(r-maj.r)<=a)   #to locate the furthest point from center
  if(length(set)==0) {set = which(abs(r-maj.r)<=2*a)}
  #randomly pick a point from this set
  p = s[sample(set,1),] 
  #compute angle of rotation
  theta = -atan((p[1]-center[1])/(p[2]-center[2]))
  return(list(r=r,p=p,theta=theta,min.r=min.r,maj.r=maj.r))
}

r = ellip_axes(z_lr,c(300,529))$r
p = ellip_axes(z_lr,c(300,529))$p
summary(r)


##try 
theta = atan((p[1]-center[1])/(p[2]-center[2]))

#plot ellipse
t = seq(0,2*pi,0.01)
major = 50
minor = 10
phi = atan((330-295)/(600-535))
#x axis is the column number
ellp_x = 535 + major*cos(t)*cos(phi) - minor*sin(t)*sin(phi)
ellp_y = 295 + major*cos(t)*cos(phi) + minor*sin(t)*cos(phi)
image(1:800,1:600,t(apply(bd,2, rev)),col=c("white","magenta"))
lines(ellp_x,ellp_y,col="black",lty=2,pch=2)


