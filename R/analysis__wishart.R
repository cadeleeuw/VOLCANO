### All functions used for p-value computation ###

# the integral.func argument should just be one of bivariate.integral, multivariate.integral or pcov.integral
# for pcov.integral, omega and sigma should be formatted such that the first two pheno are X and Y, and the remainder is Z
integral.p = function(integral.func, K, omega, sigma, min.iter=10000, adap.thresh=c(1e-4, 1e-6)) {
	tot.iter = min.iter * 10^(0:length(adap.thresh))
	adap.thresh = c(adap.thresh,0)  # adding dummy 0 at the end to simplify loop code

	p = 1; curr.iter = 0
	for (i in 1:length(tot.iter)) {
		add.iter = tot.iter[i] - curr.iter
		add.p = integral.func(K, omega, sigma, n.iter=add.iter)
		p = (curr.iter*p + add.iter*add.p) / tot.iter[i]
		curr.iter = tot.iter[i]
		if (all(is.na(p)) || all(p[!is.na(p)] >= adap.thresh[i])) break
	}
	return(p)
}


# tests gamma_j = 0 for each element j of gamma
# input omega and sigma should be the whole matrix, including Y
# computing observed and null gammas internally, to save on number of required input arguments
# not doing checks on eg. invertability here, is assumed to have been checked and addressed before
multivariate.integral = function(K, omega, sigma, n.iter=1000) {
	P = dim(sigma)[1]; Px = P - 1
	omega.x = omega[-P,-P]; omega.xy = omega[-P,P]
	sig.xys = solve(sigma[1:Px,1:Px]) %*% sigma[1:Px,P] / K
	var.y = as.numeric(sigma[P,P] - sigma[P,1:Px] %*% solve(sigma[1:Px,1:Px]) %*% sigma[1:Px,P]) / K^2

	sigma.use = matrix(0,P+Px,P+Px); sigma.use[1:Px,1:Px] = sigma[1:Px,1:Px]
	theta = matrix(0,P+Px,P+Px); theta[1:Px+Px,1:Px+Px] = K*omega.x; theta[P+Px,P+Px] = K

	draws = matrixsampling::rwishart(n.iter, K, Sigma=sigma.use, Theta=theta)
	param = apply(draws, 3, multi.cond.stats, K, sigma, sig.xys, var.y)
	C1 = param[1:(Px^2),]
	C2 = param[(1:Px)+(Px^2),]
	C3 = param[(1:Px)+(Px^2+Px),]
	sds = param[(1:Px)+(Px^2+2*Px),]

	gamma.ss.obs = diag(sqrt(diag(omega.x))) %*% solve(omega.x) %*% omega.xy
	p.out = rep(NA, Px)
	for (index in 1:Px) {
		gamma.null = rep(0,Px); gamma.null[-index] = solve(omega.x[-index,-index]) %*% omega.xy[-index]
		tau.null = as.numeric(omega[P,P] - t(omega.xy[-index]) %*% solve(omega.x[-index,-index]) %*% omega.xy[-index])
		tau.null.sqrt = ifelse(tau.null > 0, sqrt(tau.null), 0) # setting negative tau to 0 here

		M = gamma.null %*% C1[1:Px+(index-1)*Px,] + tau.null.sqrt * C2[index,] + C3[index,]
		p.out[index] = conditional.norm(gamma.ss.obs[index], M, sds[index,])
	}
	return(p.out)
}


multi.cond.stats = function(draw, K, sigma, sig.xys, var.y) {
	Px = dim(sigma)[1]-1; i.eps = 1:Px; i.delta = i.eps + Px; i.y = 2*Px+1
	dtd.x = matrix(draw[i.eps,i.eps] + draw[i.eps,i.delta] + draw[i.delta,i.eps] + draw[i.delta,i.delta], ncol=Px)

	omega.x = matrix(dtd.x/K - sigma[1:Px,1:Px], ncol=Px)
	omega.x.inv = tryCatch(solve(omega.x),error=function(x){omega.x*NA}) # silently put to NA if not invertible
	O.x = suppressWarnings(diag(sqrt(diag(omega.x)), ncol=Px) %*% omega.x.inv)

	sds = suppressWarnings(sqrt(diag(var.y * O.x %*% dtd.x %*% t(O.x))))

	C1 = (draw[i.delta,i.delta] + draw[i.delta,i.eps]) %*% t(O.x) / K
	C2 = O.x %*% (draw[i.delta,i.y] + draw[i.eps,i.y])/K
	C3 = O.x %*% ((draw[i.delta,i.eps] + draw[i.eps,i.eps]) %*% sig.xys - sigma[1:Px,Px+1])

	return(c(C1,C2,C3,sds))
}


bivariate.integral = function(K, omega, sigma, n.iter=1000, add.reverse=T) {
	if (!add.reverse) {
		omega.null = diag(diag(omega))
		sig.use = matrix(0,3,3); sig.use[1,1] = sigma[1,1]
		theta = matrix(0,3,3); theta[-1,-1] = omega.null*K

		sig.xy = sigma[1,2]
		sig.xys = sig.xy/sigma[1,1]
		var.y = sigma[2,2] - (sigma[1,2]^2)/sigma[1,1]

		params = apply(matrixsampling::rwishart(n.iter, K, Sigma=sig.use, Theta=theta), 3, bivar.cond.stats, K=K, sig.xy, sig.xys, var.y)   # first row is means, second is SDs
		return(conditional.norm(omega[1,2], params[1,], params[2,]))
	} else {
		p1 = bivariate.integral(K, omega, sigma, n.iter/2, add.reverse=F)
		p2 = bivariate.integral(K, omega[2:1,2:1], sigma[2:1,2:1], n.iter/2, add.reverse=F)
		return((p1+p2)/2)
	}
}

# this is an internal function for the apply in integral.p(), defined here for clarity
# draw will be the 3x3 matrix drawn from the wishart
bivar.cond.stats = function(draw, K, sig.xy, sig.xys, var.y) {
	m = draw[2,3] + draw[1,3] + sig.xys*(draw[1,2] + draw[1,1])
	m = m/K - sig.xy

	v = var.y * (draw[2,2] + 2*draw[1,2] + draw[1,1])
	v = v / K^2
	v = ifelse(v <= 0, NA, sqrt(v))
	return(c(m,v))
}

conditional.norm = function(obs, means, sds) {
	obs = abs(obs)
	prob = suppressWarnings(pnorm(obs, mean=means, sd=sds, lower.tail=F))
	prob = prob + suppressWarnings(pnorm(-obs, mean=means, sd=sds, lower.tail=T))
	return(mean(prob, na.rm=T))
}

# x and y should be at the end of omega/sigma
pcov.integral = function(K, omega, sigma, n.iter=1000, add.reverse=T) {
	if (nrow(omega) <= 2) stop("nothing to condition on")

	if (!add.reverse) {
		P = nrow(omega); Pw = P - 1; Pz = Pw - 1;

		fit.x = omega[1:Pz,Pw] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,Pw]
		fit.y = omega[1:Pz,P] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P]
		if (omega[Pw,Pw] <= fit.x) omega[Pw,Pw] = fit.x / 0.99999               # scale var up to have R2 slightly below 1
		fit.xy = omega[1:Pw,P] %*% solve(omega[1:Pw,1:Pw]) %*% omega[1:Pw,P]
		if (omega[P,P] <= fit.xy) omega[P,P] = fit.xy / 0.99999                 # scale var up to have R2 slightly below 1

		var.y = as.numeric(sigma[P,P] - sigma[P,1:Pw] %*% solve(sigma[1:Pw,1:Pw]) %*% sigma[1:Pw,P]) / K^2
		sig.xys = solve(sigma[1:Pw,1:Pw]) %*% sigma[1:Pw,P] / K
		gamma.parts = list(
			x = solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,Pw],
			y = solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P],
			fit.x = fit.x*K,
			dw.dz.gamma = c(omega[1:Pz,P], omega[1:Pz,Pw] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P])*K
		)

		sigma.use = matrix(0,P+Pw,P+Pw); sigma.use[1:Pw,1:Pw] = sigma[1:Pw,1:Pw]
		theta = matrix(0,P+Pw,P+Pw); theta[1:Pz+Pw,1:Pz+Pw] = K*omega[1:Pz,1:Pz];
		diag(theta)[Pw+Pw:P] = K*c(omega[Pw,Pw] - fit.x, omega[P,P] - fit.y)

		draws = matrixsampling::rwishart(n.iter, K, Sigma=sigma.use, Theta=theta)
		params = apply(draws, 3, pcov.cond.stats, K, sigma, gamma.parts, sig.xys, var.y)

		pcov.obs = omega[Pw,P] - t(omega[1:Pz,Pw]) %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P]
		p = conditional.norm(pcov.obs, params[1,], params[2,])
		return(p)
	} else {
		p1 = pcov.integral(K, omega, sigma, n.iter/2, add.reverse=F)
		Pz = nrow(omega)-2; index = c(1:Pz, Pz+2, Pz+1)
		p2 = pcov.integral(K, omega[index,index], sigma[index,index], n.iter/2, add.reverse=F)
		return((p1+p2)/2)
	}
}

pcov.cond.stats = function(draw, K, sigma, gamma, sig.xys, var.y) {
	Pw = dim(sigma)[1]-1; Pz = Pw - 1
	i.eps = 1:Pw; i.eps.z = 1:Pz; i.eps.x = Pw
	i.delta = 1:Pw + Pw; i.delta.z = 1:Pz + Pw; i.x = 2*Pw; i.y = i.x+1

	dtd.w = matrix(draw[i.eps,i.eps] + draw[i.eps,i.delta] + draw[i.delta,i.eps] + draw[i.delta,i.delta], ncol=Pw)
	dtd.w[Pw,Pw] = dtd.w[Pw,Pw] + 2*t(gamma$x) %*% draw[i.delta.z,i.eps.x] + gamma$fit.x
	dtd.w[-Pw,Pw] = dtd.w[-Pw,Pw] + (draw[i.delta.z,i.delta.z] + draw[i.eps.z,i.delta.z]) %*% gamma$x
	dtd.w[Pw,-Pw] = dtd.w[-Pw,Pw]

	omega.w = matrix(dtd.w/K - sigma[1:Pw,1:Pw], ncol=Pw)
	omega.z.inv = tryCatch(solve(omega.w[-Pw,-Pw]),error=function(x){omega.w[-Pw,-Pw]*NA})  # silently put to NA if not invertible
	b = matrix(c(-(omega.z.inv %*% omega.w[1:Pz,Pw]), 1), ncol=1)

	dhw.dy = gamma$dw.dz.gamma + draw[i.eps,i.delta.z] %*% gamma$y + draw[i.eps,i.y]
	dw.ew = rbind(draw[i.delta.z,i.eps], t(gamma$x) %*% draw[i.delta.z,i.eps] + draw[i.x,i.eps])

	M = dhw.dy/K + (dw.ew + draw[i.eps,i.eps]) %*% sig.xys - sigma[1:Pw,Pw+1]
	V = t(b) %*% dtd.w %*% b * var.y

	return(c(t(b) %*% M, ifelse(V >= 0, sqrt(V), NA)))
}



