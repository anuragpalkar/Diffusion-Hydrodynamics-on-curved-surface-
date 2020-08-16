lap_t_u = function(u)
{
	h = 0.01
	k = 0.01
	dt = 0.1
	a = dt/2/h^2
	a = .1
	u_next = u
	m = nrow(u)
	n = ncol(u)
	Trm = matrix(rep(0,n*n),n,byrow=T)
	for (i in 2:(n-1)) {Trm[i,i-1] = Trm[i,i+1] = a; Trm[i,i] = 1 - 4*a}
	Trm[1,1] = Trm[n,n] = 1 - 4*a
	Trm[n,1] = Trm[1,n] = Trm[1,2] = Trm[n,n-1] = a
	u_next1 = u %*% Trm
	u1 = u * 0
	u1[1,] = u[2,]
	u1[m,] = u[m-1,]
	Trm1 = diag(a, n)
	for (i in 1:n) 
	{
		if (n %% 2 == 0)
		{
			if (i == n/2)
			inext = n
			else
			inext = (i + n/2) %% n
			Trm1[inext,i] = a
		}
		else
		{
			if (i <= (n - 1) /2)
			Trm1[(i + (n - 1)/2),i] = a
			else
			Trm1[(i + (n + 1)/2) %% n,i] = a
		}
	}
	u_next1 = u_next1 + u1 %*% Trm1
	
	u2 = t(u)
	Trm2 = matrix(rep(0,m*m),m)
	for (i in 2:(m-1)) {Trm2[i+1,i] = Trm2[i-1,i] = a}
	u_next1 = u_next1 + t(u2 %*% Trm2)
	u_next1[which(u_next1 < 1e-17)] = 0	
	return (u_next1)
}

del_th = function(v)
{
	Trm = matrix(0,ncol(v),ncol(v))
	for ( i in 1:(ncol(v) - 1))
	Trm[i + 1, i] = 1
	for ( i in 1:(ncol(v) - 1))
	Trm[i, i + 1] = -1
	Trm[ncol(v), 1] = -1
	Trm[1, ncol(v)] = 1
	dv = v %*% Trm
	return (dv)
}

del_gam = function(v)
{
	v1 = v * 0
	m = nrow(v)
	n = ncol(v)
	v1[1,] = v[2,]
	v1[m,] = v[m-1,]
	Trm1 = diag(1, n)
	for (i in 1:n) 
	{
		if (n %% 2 == 0)
		{
		if (i == n/2)
		inext = n
		else
		inext = (i + n/2) %% n
		Trm1[inext,i] = -1
		}
		else
		Trm1[(i + (n + 1)/2) %% n,i] = -1
	}
	
	v2 = t(v)
	Trm2 = matrix(rep(0,m*m),m)
	for (i in 2:(m-1)) 
	{
		Trm2[i+1,i] = 1 
		Trm2[i-1,i] = -1
	}
	
	dv = v1 %*% Trm1 + t(v2 %*% Trm2)
	return (dv) 
}


div = function(v_theta, v_gamma)
{
	div = del_th(v_theta) + del_gam(v_gamma)
	return(div)
}

