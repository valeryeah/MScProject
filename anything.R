#boundary

bd = which(x==1,arr.ind = TRUE)
#sort according to x and y respective
#ind1 = sort(bd[,1],index.return=TRUE)$ix
#ind2 = sort(bd[,2])
ellip.2 <- function(bd){
  center = fit.ellipse(x=bd[,2],y=bd[,1])$center
  ind = which(bd[,1]<center[2],arr.ind = TRUE)
  index1 = sample(ind,round(length(ind)*0.3))
  bd1 = bd[index1,]
  ind2 = (1:nrow(bd))[-ind]
  index2 = sample(ind2,round(length(ind2)*0.3))
  bd2 = bd[index2,]
  #fit an ellipse to each bd (center,major axis,minor axis, angle)
  efit1 = fit.ellipse2(x= bd1[,2],y=bd1[,1])
  efit2 = fit.ellipse2(x=bd2[,2],y=bd2[,1])
  return(list(efit1=efit1,efit2=efit2))
}

ellip.mc <-function(N){
  para1 = matrix(0,nrow = N,ncol=5)
  para2 = matrix(0,nrow = N,ncol=5)
  for (i in 1:N) {
    v = ellip.2(bd)
    para1[i,]=v$efit1
    para2[i,]=v$efit2
  }
  return(list(para1=para1,para2=para2))
}

vv = ellip.mc(100)
para1 = apply(vv$para1,2,mean)
para2 = apply(vv$para2,2,mean)


##plotting
get.ellipse2 <- function(para, n=360 ) {
  center = para[1:2]
  major = para[3]
  minor = para[4]
  angle = para[5]
  tt <- seq(0, 2*pi, length=n) 
  sa <- sin(angle) 
  ca <- cos(angle) 
  ct <- cos(tt) 
  st <- sin(tt) 
  x <- center[1] + major * ct * ca - minor * st * sa 
  y <- center[2] + major * ct * sa + minor * st * ca 
  cbind(x=x, y=y) 
}

e1 <- get.ellipse2(para1)
e2 <- get.ellipse2(para2)

#bd
plot(x = set[,2],y=set[,1],cex=.2,xlim=c(0,800),ylim=c(0,600),xlab="",ylab="",main="Fitting Ellipse") 
#fit one ellipse
lines(x=e[,1],y=e[,2],col="red")

lines(x = e1[,1],y=600-e1[,2],col="green") 
lines(x=e2[,1],y=600-e2[,2],col="blue")

###contour plot
library(mvtnorm)
library(MASS)


bd2 = bd
bd2[,1] = bd[,2]
bd2[,2] = 600-bd[,1]

ggm = function(w){
  n = nrow(bd2)
  z1 = matrix(0,nrow=n,ncol=1)
  z2 = matrix(0,nrow=n,ncol=1)
  mu1 = para1[1:2]
  sigma1 = matrix(c(1,0,0,1),2,2)
  mu2 = para2[1:2]
  sigma2 = matrix(c(1,0,0,1),2,2)
  for (i in 1:n) {
      z1[i] = dmvnorm(bd2[i,],mean=mu1,sigma=sigma1)
      z2[i] = dmvnorm(bd2[i,],mean=mu2,sigma=sigma2)
  }
  z = w*z1+(1-w)*z2
  return(z)
}
  
z = ggm(0.5)
z1 = ggm(1)
sd(z1)


#x = seq(420,600,length.out = 200)
#y = seq(220,380,length.out = 200)
#contour(x,y,z)

#library(mclust)


#plot.Mclust(m,what="density",xlim=c(0,800),ylim=c(0,600))
#lines(x = set[,2],y=set[,1],cex=.2,col="blue") 
#plot(x = set[,2],y=set[,1],cex=.2,xlim=c(0,800),ylim=c(0,600),xlab="",ylab="",main="Fitting Ellipse") 
#lines(m)
bd2 = bd
bd2[,1] = bd[,2]
bd2[,2] = 600-bd[,1]
m=Mclust(bd2,G=2)



mu1 = para1[1:2]
sigma1 = matrix(c(para1[3],2,2,para1[4]),2,2)
mu2 = para2[1:2]
sigma2 = matrix(c(para2[3],2,2,para2[4]),2,2)
for (i in 1:200) {
  for (j in 1:200) {
    z1[i,j] = dmvnorm(c(x[i],y[j]),mean=mu1,sigma=sigma1)
    z2[i,j] = dmvnorm(c(x[i],y[j]),mean=mu2,sigma=sigma2)
  }
}











