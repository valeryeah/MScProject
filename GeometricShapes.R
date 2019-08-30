#############################    fitting circle ########################################
##fitting circle to detected sandstorm after image preprocessing
#centre as the bary-centre
circle_cent = function(lr_response){
  center=c(0,0)
  for (i in 1:dim(lr_response)[1]) {
    for (j in 1:dim(lr_response)[2]) {
      center= center+lr_response[i,j]*c(i,j)
    }
  }
  center=center/sum(lr_response)
  return(center)
}

#to find radius
circle_rad = function(cent,set){
  d = rep(0,nrow(set))
  for (i in 1:nrow(set)) {
    d[i] = dist(rbind(cent,set[i,]))
  }
  r = max(d)
  return(r)
}



#############################    fitting ellipse  ###################################

fit.ellipse <- function (x, y = NULL) {
  # from:http://r.789695.n4.nabble.com/Fitting-a-half-ellipse-curve-tp2719037p2720560.html
  EPS <- 1.0e-8 
  dat <- xy.coords(x, y) 
  
  D1 <- cbind(dat$x * dat$x, dat$x * dat$y, dat$y * dat$y) 
  D2 <- cbind(dat$x, dat$y, 1) 
  S1 <- t(D1) %*% D1 
  S2 <- t(D1) %*% D2 
  S3 <- t(D2) %*% D2 
  T <- -solve(S3) %*% t(S2) 
  M <- S1 + S2 %*% T 
  M <- rbind(M[3,] / 2, -M[2,], M[1,] / 2) 
  evec <- eigen(M)$vec 
  cond <- 4 * evec[1,] * evec[3,] - evec[2,]^2 
  a1 <- evec[, which(cond > 0)] 
  f <- c(a1, T %*% a1) 
  names(f) <- letters[1:6] 
  
  # calculate the center and lengths of the semi-axes 
  # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2288654/
  A <- matrix(c(2*f[1], f[2], f[2], 2*f[3]), nrow=2, ncol=2, byrow=T )
  b <- matrix(c(-f[4], -f[5]), nrow=2, ncol=1, byrow=T)
  soln <- solve(A) %*% b
  b2 <- f[2]^2 / 4
  center <- c(soln[1], soln[2]) 
  names(center) <- c("x", "y") 
  
  num  <- 2 * (f[1] * f[5]^2 / 4 + f[3] * f[4]^2 / 4 + f[6] * b2 - f[2]*f[4]*f[5]/4 - f[1]*f[3]*f[6]) 
  den1 <- (b2 - f[1]*f[3]) 
  den2 <- sqrt((f[1] - f[3])^2 + 4*b2) 
  den3 <- f[1] + f[3] 
  
  semi.axes <- sqrt(c( num / (den1 * (den2 - den3)),  num / (den1 * (-den2 - den3)) )) 
  
  # calculate the angle of rotation 
  term <- (f[1] - f[3]) / f[2] 
  angle <- atan(1 / term) / 2 
  
  list(coef=f, center = center, major = max(semi.axes), minor = min(semi.axes), angle = unname(angle)) 
}

##function used for plotting 
get.ellipse <- function(fit, n=360 ) {
  tt <- seq(0, 2*pi, length=n) 
  sa <- sin(fit$angle) 
  ca <- cos(fit$angle) 
  ct <- cos(tt) 
  st <- sin(tt) 
  
  x <- fit$center[1] + fit$major * ct * ca - fit$minor * st * sa 
  y <- fit$center[2] + fit$major * ct * sa + fit$minor * st * ca 
  #y <- fit$center[1] + fit$major * ct * ca - fit$minor * st * sa 
  #x <- fit$center[2] + fit$major * ct * sa + fit$minor * st * ca 
  cbind(x=x, y=y) 
}

