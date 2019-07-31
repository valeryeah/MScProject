
############################   fixed parameter   ####################################
fix.arima = arima(ep.cx,order=c(1,1,2))
fix.arima$coef


## major and minor axis
ep.maj = ep[,3]
plot(ep.maj,type="l",lty=2)
ep.min = ep[,4]
plot(ep.min,type="l",lty=2)

auto.arima(ep.maj)
arima(ep.maj,order=c(1,1,3))

auto.arima(ep.min)


############################   adaptive estimation   ################################
##parameter of ellipses in ep
#first plotting x y coordinate of center
ep.cx = ep[,1]
plot(ep.cx,type="l",lty=2)
ep.cy = ep[,2]
plot(ep.cy,type="l",lty=2)


library("tseries")
#first order difference of x 
ep.cxd = diff(ep.ex,differences=1) 
cx20 = ep.cx[21:40]
adf.test(cx20, alternative = "stationary")
cx20d = diff(cx20,differences=1)
adf.test(cx20d, alternative = "stationary")

#y_t - y_{t-1} = alpha*(y_{t-1} - y_{t-2}) + beta* w_t + gamma* w_{t-1}
#noise w



arima112 = function(theta,y){
  n = length(y)
  y = y[n] + theta[1] * (y[n]- y[n-1]) - theta[2]*w[n+1] - theta[3]*w[n]
  return(y)
}



set.seed(1023)
w = rnorm(100,sd=0.5)
train = ep.cy[40:60]
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
    d1 = ifelse(dif>0,-delta*(y[i-1] - y[i-2]),delta*(y[i-1] - y[i-2]))
    d2 = ifelse(dif>0,-delta*w[i-1],delta*w[i-1])
    #d3 = ifelse(dif>0,-delta*w[i-2],delta*w[i-2])
    theta[1] <- theta[1] + d1*abs(dif)
    theta[2] <- theta[2] + d2*abs(dif)
    #theta[3] <- theta[3] + d3*abs(dif)
    theta.all[i,] = theta
  }
  return(list(theta.all=theta.all,y.pred=y.pred))
  #s = sum((dif)^2)
  #return(s)
}
#optim(par=0.1,fn=adpt.est,y,w)

ae = adpt.est(0.08,train,ep.cy[59:65],w)$y.pred 
#theta.all =  adpt.est(0.08,ep.cy[49:85],w)$theta.all
#d = y.pred - ep.cy[51:85]
#sum(abs(d))



## major and minor axis
ep.maj = ep[,3]
plot(ep.maj,type="l",lty=2)
ep.min = ep[,4]
plot(ep.min,type="l",lty=2)




############################   sliding window    ####################################
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

sw= sliding.window(20,65-60,ep.cy[40:65])

plot(ep.cy[61:65],type="l",lty=2,lwd=2,ylim=c(274,288),main="Prediction of c_y at t=61:65")
#sliding window
lines(sw,type="l",col="red")
#adaptive 
lines(ae,col="blue")




############################   fixed parameter  ####################################

#n is the number of predictions
fixed.para = function(n,train){
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

fixed= fixed.para(65-60,ep.cy[40:60])

lines(fixed,col="green")


train = ep.cx[35:50]
ae = adpt.est(0.08,ep.cx[49:85],w)$y.pred
sw = sliding.window(15,85-50,ep.cx[35:85])
fixed = fixed.para(85-50,ep.cx[35:50])

plot(ep.cx[51:80],type="l",lwd=2)
lines(sw,type="l",col="red")
lines(ae,col="blue")
lines(fixed,col="green")
lines(ep.cx[35:50],col="mistyrose")

