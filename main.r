library(rgl)
library(colorspace)
library(plotrix)
source('get_geometry.r')
source('diffusion_functions.r')


# get geometry for diffusion surface
geo = ellipsoid(1,1)
colf = geo[,4:5]
geo = geo[,1:3]
write.csv(geo,"sim/geometry.csv")
u = matrix(runif(length(unique(colf[,1]))*length(unique(colf[,2])))*1+10,length(unique(colf[,2])),length(unique(colf[,1])),dimnames=list(unique(colf[,1]),unique(colf[,2])))
v_theta = matrix(rep(0,length(colf[,2])),length(unique(colf[,1])),dimnames=list(unique(colf[,1]),unique(colf[,2])))
v_gamma = matrix(rep(0,length(colf[,2])),length(unique(colf[,1])),dimnames=list(unique(colf[,1]),unique(colf[,2])))

# set parameters and apply diffusion equations on the set geometry
D = 1
eta = 5
zeta = 100
COF = 1 
t = 0
while(t<5000)
{
	if (t %% 100 == 0)
	{
		write.csv(u,paste0("sim/dif",t,".csv"))
		v = sqrt(v_gamma^2 + v_theta^2)
		write.csv(v,paste0("sim/vel",t,".csv"))
	}

	div_flux = D * lap_t_u(u)  - div(v_theta, v_gamma)*u
	c = div_flux
	print(head(c))
	sig_aa = eta * del_th(v_theta) + c /(1 + c)
	sig_ab = eta * del_gam(v_theta) + c /(1 + c)
	sig_ba = eta * del_th(v_gamma) + c /(1 + c)
	sig_bb = eta * del_gam(v_gamma) + c /(1 + c)
	div_sig_theta = div(sig_aa, sig_ba) 
	div_sig_gamma = div(sig_ab, sig_bb)
	v_theta =  div_sig_theta / zeta
	v_gamma =  div_sig_gamma / zeta
	u = c	 
	
	#print(head(c))
	t = t + 1
}

# get diffusion (differential equations) images 3D
show_sim = function()
{
	t = 0
	geo = read.csv("sim/geometry.csv",row.names=1)
	while (T)
	{
	c = read.csv(paste0("sim/dif",t,".csv"),row.names=1)
	d=c()
	for (i in 1:ncol(c))
	d = c(d,c[,i])
	leg = round(seq(min(d),max(d),(max(d)-min(d))/6),digits=4)
	b = (d - mean(d))/(max(d) - mean(d))*255/2+255/2
	b[b<0]=0
	b[b>255]=255
	palette <- colorRampPalette(c("yellow", "green",  "red","blue", "gray", "black"))
	col.table <- palette(256)
	col.index <- cut(b, 256)
	par3d("windowRect" = c(0,0,700,700))		
	layout3d(matrix(1:2, 1,2), c(0.8, 0.2), 1)	
	plot3d(geo,col = col.table[col.index],aspect = c(1,1,1))		
	next3d()
	bgplot3d({
  	plot.new()
  	color.legend(0.1, 0.1, 0.9, 0.9, rect.col=col.table, legend=leg, gradient="y", cex = 1.5)})	
	rgl.snapshot(paste0("makegif/im",t,".png"),"png")
	t=t+100
	Sys.sleep(.5)
	}
}

# get diffusion (differential equations) images 2D
show_sim_2D = function()
{
	t = 0
	while (T)
	{
	c = read.csv(paste0("sim/dif",t,".csv"),row.names=1)
	d=c()
	for (i in 1:ncol(c))
	d = c(d,c[,i])
	leg = round(seq(min(d),max(d),(max(d)-min(d))/4),digits=4)-10
	b = (d - mean(d))/(max(d) - mean(d))*255/2+255/2
	b[b<0]=0
	b[b>255]=255
	palette <- colorRampPalette(c( "red", "orange", "yellow","white"))
	col.table <- palette(256)
	col.index <- cut(b, 256)
	x = rep(1:ncol(c),ncol(c))
	y=c()
	for(i in 1: ncol(c))
	y = c(y, rep(i,ncol(c)))
	par3d("windowRect" = c(0,0,700,700))		
	layout3d(matrix(1:2, 1,2), c(0.8, 0.2), 1)
	plot3d(y,0,x,col = col.table[col.index],aspect = c(1,1,1))
	next3d()
	bgplot3d({
  	plot.new()
  	color.legend(0.1, 0.1, 0.9, 0.9, rect.col=col.table, legend=leg, gradient="y", cex = 1.5)})	
	rgl.snapshot(paste0("makegif/im",t,".png"),"png")
	t=t+100
	Sys.sleep(.5)
	}
}

# get diffusion (differential equations) contour images
show_sim_contour = function()
{
	t = 0
	while (T)
	{
	d = read.csv(paste0("sim/dif",t,".csv"),row.names=1)
	vel = read.csv(paste0("sim/vel",t,".csv"),row.names=1)
	c = v = c()
	for (i in 1:ncol(d))
	{
		c = c(c,d[,i])
		v = c(v,vel[,i])
	}
	c = matrix(c, nrow(d), byrow=T)
	v = matrix(v, nrow(vel), byrow=T)
	jpeg(paste0("fim",t,".jpg"),width=720,height=720,pointsize=17)
	layout(matrix(1:2, 1,2), c(0.8, 0.2), 1)
	image(c,  ylab = "y / l", main = "Advect", xlab = "x / l")
	par(new = TRUE)
	contour(v,col='blue',lwd =1,nlevels = 6, lty = 2,add = T, xaxt = "n", yaxt = "n",
      ylab = "", xlab = "")
	axis(side = 4)
	#mtext("velocity", side = 4, line = 3)
	legend("topleft", c("c", "vel."), col = c("red", "blue"), lty = c(1, 2), lwd = c(2,2))	
	dev.off()
	t=t+100
	Sys.sleep(.5)
	}
}

#generate gifs from generated images
gg =function()
{
	t = 0
	m = NULL
	geo = read.csv("sim/geometry.csv",row.names=1)
	while (t<102e3)
	{
	c = read.csv(paste0("sim/dif",t,".csv"),row.names=1)
	d=c()
	for (i in 1:ncol(c))
	d = c(d,c[,i])
	cols = rgb(20,(d - abs(min(d)))/(max(d) - abs(min(d)))*255,0,maxColorValue=255)
	par3d("windowRect" = c(0,0,700,700))
	plot3d(geo,col = cols,aspect = c(1,1,1))
	rgl.snapshot(paste0("makegif/im",t,".png"),"png")
	frame <- magick::image_read(paste0("makegif/im",t,".png"))
	if (is.null(m)) 
	m <- frame
	else m <- c(m, frame)		
	t=t+1000
	#Sys.sleep(.5)
	}
	m <- magick::image_animate(m, fps = 5, dispose = "previous")
	magick::image_write(m, "makegif/movie1.gif")
}
