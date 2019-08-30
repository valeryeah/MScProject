
############################   fixed parameter  estimation ####################################
##ARIMA(1,1,1) forecasts when parameters are fixed
#n is the number of predictions, train is the training data, w is the noise vector
fixed.para = function(n,train,w){
  pred = rep(0,n)
  mod = arima(train,order=c(1,1,1),method="ML")
  theta = mod$coef
  pred[1] = train[length(train)] + theta[1]*(train[length(train)]-train[length(train)-1]) +w[3]+ theta[2]*w[2]
  pred[2] = pred[1] + theta[1]*(pred[1] - train[length(train)])+w[4] - theta[2]*w[3]
  for (i in 3:n) {
    pred[i] = pred[i-1] + theta[1]*(pred[i-1]-pred[i-2]) +w[i+2]+ theta[2]*w[i+1]
  }
  return(pred)
}


############################   sliding window    ####################################
##ARIMA(1,1,1) forecasts when parameters are estimated by sliding window
#l is the length of window
#n is the number of prediction
#dat is the data
sliding.window = function(l,n,dat){
  pred = rep(0,n)
  for (i in 1:n) {
    x = dat[i:(l+i-1)]
    mod = arima(x,order=c(1,1,1),method="ML")
    theta = mod$coef
    #using last two obs in window to make predictions 
    pred[i] = x[l] + theta[1] * (x[l]- x[l-1]) +w[i+2]+ theta[2]*w[i+1] #+ theta[3]*w[i]
  }
  return(pred)
}

############################   adaptive estimation   ################################
##ARIMA(1,1,1) forecasts when parameters are estimated by adaptive estimation
#delta is the step size
#y is the observation set, also including the last two obs in training data
#w is the noise
adpt.est = function(delta,train,y,w){
  y.pred = rep(0,length(y)-2)  ##prediction of y_t
  dif = rep(0,length(y)-2)
  theta.all = matrix(0,nrow=length(y),ncol=2)
  mod = arima(train,order=c(1,1,1),method="ML")
  theta = mod$coef
  theta.all[1,] = theta
  for (i in 3:length(y)) {
    #predict y_t, i.e. y_{i+2}
    y.pred[i-2] = y[i-1] + theta[1] * (y[i-1]- y[i-2]) +w[i] + theta[2]*w[i-1] #+ theta[3]*w[i-2]
    dif = y.pred[i-2] - y[i]
    d1 = abs(delta*(y[i-1] - y[i-2]))
    theta[1] <- theta[1] - d1*sign(y[i-1]-y[i-2])

    theta.all[i,] = theta
  }
  return(list(theta.all=theta.all,y.pred=y.pred))
}


