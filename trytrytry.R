img3 = readJPEG("201706241615_dst.jpg")
v = fit.lr(img3)
v = fit.glm(img3)
v2 = find_bound(v$lr_response)
center1 = ellip_center1(v$lr_response)
#center2 = ellip_center2(v$lr_response)

z_thre = bin_class(img3)
image(1:800,1:600,t(apply(z_thre,2, rev)),col=c("white","magenta"))

plot_fitted(v$lr_response)


image(1:800,1:600,t(apply(v2$bd,2, rev)),col=c("white","magenta"))
#locate center
center2[1] = 600-center2[1]
y = function(x) x+ center2[1]-center2[2]
t = seq(center2[2]-1,center2[2]+1,by=1)
lines(t,y(t),col="blue",lwd=3)

z_lr = v2$z_lr
r =  dist.to.center(lr_response,center)
summary(r)

##with points lying inside and outside the ellipse, use quantile p and q for minor and major axes

set.seed(5234)
center1= ellip_center1(v$lr_response)
v3 = ellip_para(v$lr_response,center1,0.1,0.9)
p = v3$p
p[1] = 600-p[1]
center1[1] = 600-center1[1]

#plotting fitted ellipse
t = seq(0,2*pi,0.01)
major = v3$maj.r
minor = v3$min.r
phi = v3$theta
#x axis is the column number
ellp_x = center1[2] + major*cos(t)*cos(phi) - minor*sin(t)*sin(phi)
ellp_y = center1[1] + major*cos(t)*sin(phi) + minor*sin(t)*cos(phi)
image(1:800,1:600,t(apply(v2$bd,2, rev)),col=c("white","magenta"))
lines(ellp_x,ellp_y,col="black",lty=2,pch=2)

plot(ellp_x,ellp_y,col="magenta",lty=2,lwd=0.2,xlim=c(0,800),ylim=c(0,600))
rasterImage(raster,0,800,0,600)



