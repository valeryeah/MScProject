###Time Series
#for now just consider the first vortex in the list
files <- list.files(path="/Users/mac/Desktop/MSc Project/ts", pattern=".jpg",all.files=T, full.names=T) 
list_of_images = lapply(files, readJPEG)


#input is a list of images
#output is parameters of fitted ellipses
extract.para = function(images){
  img = images[[1]]
  f = fit.glm(img)
  e = edge.detection(f)
  i = image.seg(e$mag,e$set)
  efit= fitted.ellip(i$l,i$n)
  ep.para = list(efit[[1]])
  for (i in 2:length(images)) {
    img = images[[i]]
    f = fit.glm(img)
    e = edge.detection(f)
    i = image.seg(e$mag,e$set)
    efit= fitted.ellip(i$l,i$n)
    ep.para = list.append(ep.para,efit[[1]])
  }
  return(ep.para)
}


ts = extract.para(list_of_images)



######################### EXTRACT ELLIPSE PARAMETER #################################


##parameters of fitted ellipse 1("x==1")
ts.para = function(mag,set){
  x = as.matrix(bwlabel(mag))  #generate labels to objects in the binary image
  set = which(x==1,arr.ind = TRUE)
  set[,1] = 600-set[,1]
  efit = fit.ellipse2(x= set[,2],y=set[,1])
  return(efit)
}

#input is list of images
#output is parameters of fitted ellipse
ep.para = function(images){
  para=matrix(0,nrow=length(images),ncol=5)
  for (i in 1:length(images)) {
    img = images[[i]]
    f = fit.glm(img)
    e = edge.detection(f)
    x = as.matrix(bwlabel(e$mag))  #generate labels to objects in the binary image
    set = which(x==1,arr.ind = TRUE)
    set[,1] = 600-set[,1]
    para[i,] = fit.ellipse2(x= set[,2],y=set[,1])
  }
  return(para)
}


ep=ep.para(list_of_images)

# outlier
img56 = list_of_images[[77]]
f56 = fit.glm(img56)
e56 = edge.detection(f56)
ep[77,] = ts.para(e56$mag,e56$set)



library(forecast)
ep.train = ep[c(1:69),]
ep.test = ep[70:85,]



################################ TRY TRY TRY #####################################
ar(ep.train[,3],order.max = 1,method="yule-walker")
ar.burg(ep.train,order.max = 1,method="yw")


#univariate times series for ellipse parameter
centx.ts = arima(ep.train[,1],order= c(0,1,0))
centy.ts = arima(ep.train[,2],order= c(0,1,0))
majr.ts = arima(ep.train[,3],order= c(2,2,2))
minr.ts = arima(ep.train[,4],order= c(0,1,0))
ang.ts = arima(ep.train[,5],order= c(2,2,2))

#forecasting on h=1
centx.fc = forecast(centx.ts,h=1)$mean[1]
centy.fc = forecast(centy.ts,h=1)$mean[1]
majr.fc = forecast(majr.ts,h=1)$mean[1]
minr.fc = forecast(minr.ts,h=1)$mean[1]
ang.fc = forecast(ang.ts,h=1)$mean[1]






###using the whole sandstorm, rather than the boundary!
ss.para = function(images){
  para=matrix(0,nrow=length(images),ncol=4)
  for (i in 1:length(images)) {
    img = images[[i]]
    fit = fit.glm(img)
    cimg = as.cimg(fit$lr_response)   #convert to cimg 
    blur = isoblur(cimg,sigma=5)  #blur image
    eps=0.03
    img.blurred = blur>eps
    lab = bwlabel(img.blurred)
    set = which(lab=="1",arr.ind = TRUE)[,1:2]
    sam.set = sample(nrow(set),round(0.7*nrow(set)))
    #compute mean and sd, both in x y directions
    para[i,1:2] = apply(set[sam.set,],2,mean)
    para[i,3:4] = apply(set[sam.set,],2,sd)
  }
  return(para)
}

ss.para = ss.para(list_of_images)
bd.train = ss.para[1:69,]
bd.test = ss.para[70:85,]



set.seed(1284)
lt = function(bd.train){
  n = dim(bd.train)[1]
  img = list_of_images[[n]]
  fit = fit.glm(img)
  cimg = as.cimg(fit$lr_response)   #convert to cimg 
  blur = isoblur(cimg,sigma=5)  #blur image
  eps=0.03
  img.blurred = blur>eps
  lab = bwlabel(img.blurred)
  set = which(lab=="1",arr.ind = TRUE)[,1:2]
  ind = sample(nrow(set),round(0.7*nrow(set)))
  sam.set = set[ind,]
  m = apply(set,2,mean)
  return(list(m=m,sam.set=sam.set))
}

m.x = arima(bd.train[,1],order= c(2,2,2))
m.y = arima(bd.train[,2],order= c(2,2,2))
s.d.x = arima(bd.train[,3],order= c(1,1,1))
s.d.y = arima(bd.train[,4],order= c(1,1,1))





############################### boundary mean & sd ###################################


#compute mean and sd of RANDOMLY SELECTED boundary points for each image
bd.para = function(images){
  para=matrix(0,nrow=length(images),ncol=4)
  for (i in 1:length(images)) {
    img = images[[i]]
    f = fit.glm(img)
    e = edge.detection(f)
    x = as.matrix(bwlabel(e$mag)) 
    set = which(x==1,arr.ind = TRUE)
    #randomly select 70% boundary points
    sam.set = sample(nrow(set),round(0.7*nrow(set)))
    #compute mean and sd, both in x y directions
    para[i,1:2] = apply(set[sam.set,],2,mean)
    para[i,3:4] = apply(set[sam.set,],2,sd)
  }
  return(para)
}

bd.para = bd.para(list_of_images)
x56 = as.matrix(bwlabel(e56$mag))
set56 = which(x56==2,arr.ind = TRUE)
sam.set56 = sample(nrow(set56),round(0.7*nrow(set56)))
bd.para[56,1:2] = apply(set[sam.set56,],2,mean)
bd.para[56,3:4] = apply(set[sam.set56,],2,sd)


#training and test data
#col 1&2 mean of bd
#col 3&4 sd of bd
bd.train = bd.para[1:69,]  
bd.test = bd.para[70:85,]


#last image in training set
set.seed(1284)
last.train = function(bd.train){
  n = dim(bd.train)[1]
  img = list_of_images[[n]]
  fit = fit.glm(img)
  e = edge.detection(fit)
  x = as.matrix(bwlabel(e$mag)) 
  set = which(x==1,arr.ind = TRUE)
  ind = sample(nrow(set),round(0.7*nrow(set)))
  sam.set = set[ind,]
  m = apply(set,2,mean)
  return(list(m=m,sam.set=sam.set))
}

lt69 = last.train(bd.train)
m69 = lt69$m
sam.set69 = lt69$sam.set

lt50 = list(m=m,sam.set=sam.set)

############################### UNIVARIATE TS ###################################

#auto.arima
m.x = arima(bd.train[,1],order= c(2,2,2))
m.y = arima(bd.train[,2],order= c(2,2,2))
s.d.x = arima(bd.train[,3],order= c(1,1,1))
s.d.y = arima(bd.train[,4],order= c(1,1,1))


#forecast h step ahead
#returns forecasts of mean and sd of boundary center 
mv.fc = function(h){
  m.x.fc = forecast(m.x,h)$mean[h]
  m.y.fc = forecast(m.y,h)$mean[h]
  sd.x.fc = forecast(s.d.x,h)$mean[h]
  sd.y.fc = forecast(s.d.y,h)$mean[h]
  m.fc = c(m.x.fc,m.y.fc)
  sd.fc = c(sd.x.fc,sd.y.fc)
  return(list(m.fc=m.fc,sd.fc=sd.fc))
}


library(fGarch)

g = garchFit(~garch(1,2),data=bd.train[,1])
predict(g,n.ahead=10,trace=FALSE)




############################### MULTIVARIATE TS ###################################
library(marima)
library(vars)


############################### Metropolis Hastings ###################################
##adaptive estimation of bd center 
set.seed(9637)
w = rnorm(100)
bdx = adpt.est(bd.train[1:50,1],0.05,w)$y.pred  ##predict t=51 
bdy = adpt.est(bd.train[1:50,2],0.05,w)$y.pred
bdxsd = adpt.est(bd.train[1:50,3],0.05,w)$y.pred  ##predict t=51 
bdysd = adpt.est(bd.train[1:50,4],0.05,w)$y.pred

m.fc = c(bdx[length(bdx)],bdy[length(bdy)])
sd.fc = c(bdxsd[length(bdx)],bdysd[length(bdy)])
fc = list(m.fc=m.fc,sd.fc=sd.fc)

#preperation for metropolis hastings
#returns functions f&q for MCMC
mh.prep = function(h,lt,fc=fc){
  #fc = mv.fc(h)
  #add adjustment for sampled points
  set.seed(9637)
  r = rnorm(2,sd=0.5)
  d.fc = fc$m.fc-lt$m + r
  #new.set = lt$sam.set + matrix(data=d.fc,nrow=dim(lt$sam.set)[1],ncol=2) 
  new.set = lt$sam.set + d.fc
  #new kde 
  dens = kde(new.set)
  #target function
  f <- function(y) predict(dens,x=y) 
  #proposal distribution, multivairate normal with sd predicted by TS
  prop <- function(x) rmvnorm(1,mean=x,sigma=diag(fc$sd.fc))   ###UNIvariate ts
  #sig = matrix(data=c(fc$sd.fc[1],5,5,fc$sd.fc[2]),2,2)
  #q <- function(x) rmvnorm(1,mean=x,sigma=sig)
  return(list(f=f,prop=prop,m.fc=fc$m.fc))
}



cent = mv.fc(h)$m.fc
##use predicted major and minor axes of fitted ellipse as constraint
#auto.arima
majr.ts = arima(ep.train[,3],order= c(1,1,1))
minr.ts = arima(ep.train[,4],order= c(0,1,1))


majr.ts = adpt.est(ep.train[1:50,3],0.05,w)$y.pred  ##predict t=51 
minr.ts = adpt.est(ep.train[1:50,4],0.05,w)$y.pred
majr.pd = majr.ts[length(majr.ts)]
minr.pd = minr.ts[length(minr.ts)]

mh = function(n=1e4,h,constraint=T){
  x.all <- matrix(0,nrow=n,ncol=2)
  prep = mh.prep(h,lt50,fc)
  x <-prep$m.fc
  prop = prep$prop
  f = prep$f
  #cent = mv.fc(h)$m.fc
  cent = m.fc
  #majr = forecast(majr.ts,h)$mean[h]
  #minr = forecast(minr.ts,h)$mean[h]
  for (i in 1:n) {
    y = prop(x)
    #To compare with or without constraint
    if(constraint){
      #allow some relaxation, 10%
      if (runif(1) <= min(f(y)/f(x),1) && dist(rbind(cent,y))<=(majr*1.1) && dist(rbind(cent,y))>=(minr*0.9))
        #if (runif(1) <= min(f(y)/f(x),1))
        x<-y
        x.all[i,] <- x
    }
    #No constraint
    else{
      if (runif(1) <= min(f(y)/f(x),1))
      x<-y
      x.all[i,] <- x
    }
  }
  x.all
}

#compare samples with true boundary
true.fc = function(h){
  img = list_of_images[[69+h]]
  fit = fit.glm(img)
  e = edge.detection(fit)
  x = as.matrix(bwlabel(e$mag)) 
  set = which(x==1,arr.ind = TRUE)
  return(set)
}


x.all =mh(n=1e4,3,constraint = T)
dim(x.all)

x.all2 =mh(n=1e4,3,constraint = F)
dim(x.all2)

true3 = true.fc(3)


set.seed(1234)
# with constraint
x.all = mh(n=1e4,10,constraint=T)
#acceptance rate
x = x.all[,1]
n=1e4
mean(x[-1]!=x[-(n+1)])
#0.4043

plot(x.all,cex=0.2,col="red")
points(true.fc(10))


# without constraint
set.seed(1234)
x.all = mh(n=1e4,10,constraint=F)
x = x.all[,1]
mean(x[-1]!=x[-(n+1)])
###0.6584
  

plot(x.all,cex=0.2,col="red")
true.fc(10)

#compare with last train
true.fc(0)




