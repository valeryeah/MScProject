##find the area of true sandstorm by counting pixels
sd.area = function(images){
  sd.area=NULL
  for (i in 1:length(images)) {
    img = images[[i]]
    f = fit.glm(img)
    cimg = as.cimg(f$lr_response)   #convert to cimg 
    blur = isoblur(cimg,sigma=5)  #blur image
    img.blurred = blur>0.03
    lab = bwlabel(img.blurred)
    if(i==16|i==50){x1 = which(lab=="2",arr.ind = TRUE)[,1:2]}
    else{x1 = which(lab=="1",arr.ind = TRUE)[,1:2]}##all points inside
    sd.area[i] = nrow(x1)
  }
  return(sd.area)
}

##find the area of mixture model by counting pixels
#input is the mixture model and the contour level representing the boundary
#model = densityMclust(data=x1,G=2,prior=priorControl(functionName = "defaultPrior"),modelNames="EEE")
#densty = dens(modelName = gmm$modelName, data = set, parameters = gmm$parameters) 
gmm.area.fn = function(model,densty){
  x1 = rep(1:600,1,each=800)
  x2 = rep(1:800,600)
  dat = cbind(x1,x2)
  d = dens(modelName=model$modelName,data=dat,parameters=model$parameters)
  return(sum(d>=densty))
}

