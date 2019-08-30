################## extracting the parameters of each model to form times series


##time series of circle parameter
cir.para = function(images){
  para=matrix(0,nrow=length(images),ncol=3)
  for (i in 1:length(images)) {
    img = images[[i]]
    f = fit.glm(img)
    cent = circle_cent(f$lr_response)
    w = edge.detection(f)
    v5 = image.seg(w$mag,w$set)
    l = v5$l
    set = l[[1]]
    radius = circle_rad(cent,set)
    para[i,] = c(cent,radius)
  }
  return(para)
}

##ts of ellipse parameter
ep.para = function(images){
  para=matrix(0,nrow=length(images),ncol=5)
  for (i in 1:length(images)) {
    img = images[[i]]
    f = fit.glm(img)
    e = edge.detection(f)
    x = as.matrix(bwlabel(e$mag))  #generate labels to objects in the binary image
    set = which(x==1,arr.ind = TRUE)
    set[,1] = 600-set[,1]
    para[i,] = fit.ellipse(x= set[,2],y=set[,1])
  }
  return(para)
}

##exponential transformation to ellipse angle
ep = ep,para(list_of_images)
ep.ang = ep[,5]
ep.ang = 100*exp(ep.ang)


##gmm parameter
## using the whole region rather than the boundary
gmm.para = function(images){
  para=matrix(0,nrow=length(images),ncol=6)
  covar = matrix(0,nrow=length(images),ncol=4)
  for (i in 1:length(images)) {
    img = images[[i]]
    f = fit.glm(img)
    cimg = as.cimg(f$lr_response)   #convert to cimg 
    blur = isoblur(cimg,sigma=5)  #blur image
    img.blurred = blur>0.03
    lab = bwlabel(img.blurred)
    if(i==16|i==50){x1 = which(lab=="2",arr.ind = TRUE)[,1:2]}
    else{
      x1 = which(lab=="1",arr.ind = TRUE)[,1:2]}  ##all points inside
    gr = imgradient(img.blurred,"xy")   
    mag <- with(gr,sqrt(x^2+y^2))
    x = as.matrix(bwlabel(mag))
    if(i==16|i==50){  
      set = which(x==2,arr.ind = TRUE)}
    else{
      set = which(x==1,arr.ind = TRUE)   ##boundary 
    }
    x1.plot = x1
    x1.plot[,1] = 600-x1.plot[,1]
    x1[,1] = x1.plot[,2]
    x1[,2] = x1.plot[,1]
    set.plot = set
    set.plot[,1] = 600-set[,1]
    set[,1] = set.plot[,2]
    set[,2] = set.plot[,1]
    
    defaultPrior(data= x1,G=2,"EEE")
    gmm = densityMclust(data=x1,G=2,prior=priorControl(functionName = "defaultPrior"),modelNames="EEE")
    para[i,1] = gmm$parameters$pro[1]
    para[i,2:3] = gmm$parameters$mean[,1]
    para[i,4:5] = gmm$parameters$mean[,2]
    densty = dens(modelName = gmm$modelName, data = set, parameters = gmm$parameters)
    para[i,6] = mean(densty) 
    covar[i,] = gmm$parameters$variance$Sigma
  }
  return(list(para=para,covar=covar))
}

#exponential transformation to mixing coefficient
gmm = gmm.para(list_of_images)
gmm.coef = gmm$para[,1]
gmm.coef1 = 100*log(gmm.coef/(1-gmm.coef))





