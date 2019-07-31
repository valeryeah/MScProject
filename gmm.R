library(mvtnorm)
library(MASS) #to use kde2d
library(mclust)  #gmm


#kde
library(ks)


############################### GMM parameters ###################################


gmm.para = function(images){
  para=matrix(0,nrow=length(images),ncol=14)
  for (i in 1:length(images)) {
    img = images[[i]]
    gm = gmm.fit(img)$gmm.fit$parameters
    pi = gm$pro
    m = gm$mean #mean of each component in column
    v = gm$variance$sigma
    para[i,] = c(pi,m,v)  #column by column
  }
  return(para)
}

p = gmm.para(list_of_images)
#training and test data
gmm.train = p[1:69,]
gmm.test = p[70:85,]





##setting prior
defaultPrior(x1,G=2,"VVV")
dens1 = densityMclust(x1,G=2,prior=priorControl(functionName = "defaultPrior"))




############################### UNIVARIATE TS ###################################
#auto.arima
m.x = arima(bd.train[,1],order= c(2,2,2))
m.y = arima(bd.train[,2],order= c(2,2,2))
s.d.x = arima(bd.train[,3],order= c(1,1,1))
s.d.y = arima(bd.train[,4],order= c(1,1,1))







cimg = as.cimg(fit$lr_response)   #convert to cimg 
blur = isoblur(cimg,sigma=5)  #blur image
eps=0.05
img.blurred = blur>eps
lab = bwlabel(img.blurred)
x1 = which(lab=="1",arr.ind = TRUE)[,1:2]
gmm.fit2 = Mclust(x1,G=2)
gmm.fit3 = Mclust(x1,G=3)

gmm.fit = function(img,n=2){
  fit = fit.glm(img)
  cimg = as.cimg(fit$lr_response)   #convert to cimg 
  blur = isoblur(cimg,sigma=5)  #blur image
  eps=0.05
  img.blurred = blur>eps
  lab = bwlabel(img.blurred)
  x1 = which(lab=="1",arr.ind = TRUE)[,1:2]
  gmm.fit = densityMclust(x1,G=n)
  w = edge.detection(fit)
  v5 = image.seg(w$mag,w$set)
  x = v5$x  ##generated labels 
  set = which(x==1,arr.ind = TRUE)
  return(list(gmm.fit=gmm.fit,set=set))
}
g = gmm.fit(list_of_images[[20]],n=2)
plot(g$gmm.fit,what="density",xlim=c(200,400),ylim=c(300,600))
points(g$set,cex=0.1)

dens = g$gmm.fit$density
mean(dens)   ##level of contour which is closest to data 


###likelihood
img = readJPEG("201706052030_dst.jpg")
fit = fit.glm(img)
cimg = as.cimg(fit$lr_response)   #convert to cimg 
blur = isoblur(cimg,sigma=5)  #blur image
eps=0.03
img.blurred = blur>eps
lab = bwlabel(img.blurred)
x1 = which(lab=="1",arr.ind = TRUE)[,1:2]

gr = imgradient(img.blurred,"xy")   #compute gradient after the blur
mag <- with(gr,sqrt(x^2+y^2))
x = as.matrix(bwlabel(mag)) 
bd1 = which(x==1,arr.ind = TRUE)

##fit two ellipses using boundary
bd1[,1] = 600-bd1[,1]
med = median(bd[,1])
bd1.1 = bd1[which(bd1[,1]>med),]
bd1.2 = bd1[which(bd1[,1]<med),]
e1 = fit.ellipse(x=bd1.1[,2],y=bd1.1[,1])
e2 = fit.ellipse(x=bd1.2[,2],y=bd1.2[,1])
#likelihood function, theta=[pi_1,mu_1,mu_2,sigma]
mu_1.init = e1$center
mu_2.init = e2$center
sigma.init = c(e1$major,e2$minor,e2$major)
sig.init = matrix(data=c(e1$major,e2$minor,e2$minor,e2$major),2,2)
gmm.l = function(theta){
  pi_1 = theta[1]
  pi_2 = 1-pi_1
  mu_1 = theta[2:3]
  mu_2 = theta[4:5]
  sigma = matrix(data=c(theta[6],theta[7],theta[7],theta[8]),2,2,byrow=TRUE)
  gmm = function(x){pi_1*dmvnorm(x,mean=mu_1,sigma=sigma)+pi_2*dmvnorm(x,mean=mu_2,sigma=sigma)}
  l = sum(apply(x1, 1, gmm)) + dnorm(pi_1,mean=0.5,sd=0.5) + dmvnorm(mu_1,mean=mu_1.init,sig.init)+dmvnorm(mu_2,mean=mu_2.init,sig.init)
  return(-l)
}

optim(c(0.5,mu_1.init,mu_2.init,sigma.init),gmm.l)

















