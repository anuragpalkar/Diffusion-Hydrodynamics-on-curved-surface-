sqr = function(s)
{
	a = seq(-s/2,s/2,.1)
	x = as.vector(rep(a,length(a)))
	y = as.vector(sapply(a,function(i)rep(i,length(a))))
	#z = rgb(255,(50 - x^2- y^2)/50*255,255,maxColorValue=255)
	coord = cbind(x,y)
	return (coord)
}

rec = function(a,b)
{
	a = seq(-a/2,a/2,.1)
	b = seq(-b/2,b/2,.1)
	x = as.vector(rep(a,length(b)))
	y = as.vector(sapply(b,function(i)rep(i,length(a))))
	#z = rgb(255,(50 - x^2- y^2)/50*255,255,maxColorValue=255)
	coord = cbind(x,y)
	return (coord)
}

sphere = function(r)
{
	gam = seq(0, pi, .05)
	thet = seq(0,2*pi,.05)
	gama = as.vector(rep(gam,length(thet)))
	theta = as.vector(sapply(thet,function(i)rep(i,length(gam))))
	x = r*sin(gama)*cos(theta)
	y = r*sin(gama)*sin(theta)
	z = r*cos(gama)
	coord = cbind(x, y, z, gama, theta)
	return (coord)
}

ellipsoid = function(a,b)
{
	gam = seq(0, pi, .02 * min(c( a, b))/max(c( a, b)))
	thet= seq(0,2*pi,.04)
	gama = as.vector(rep(gam,length(thet)))
	theta = as.vector(sapply(thet,function(i)rep(i,length(gam))))
	x = a*sin(gama)*cos(theta)
	y = a*sin(gama)*sin(theta)
	z = b*cos(gama)
	coord = cbind(x, y, z, gama, theta)
	return (coord)
}
