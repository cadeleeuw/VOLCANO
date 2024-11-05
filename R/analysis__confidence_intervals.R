# The functions do some bounds checking on correlations truncating them to just below the bound to stop the matrixsampling::rwishart function from aborting
# In general, its assumed that the input has been validated in advance (ie. omega.x is invertible, etc.)

### BIVARIATE ###
ci.bivariate = function(K, omega, sigma, n.iter=10000) {
	S = diag(sqrt(diag(omega)))
	corrs = solve(S) %*% omega %*% solve(S)
	corrs[corrs >= 1] = 0.99999; corrs[corrs <= -1] = -0.99999; diag(corrs) = 1
	omega = S %*% corrs %*% S

	P = dim(omega)[1]; tri = lower.tri(corrs)
	out = list()
	draws = matrixsampling::rwishart(n.iter, K, Sigma=sigma/K, Theta=omega)
	if (P == 2) {
		func = function(draw, sigma) {o = draw-sigma; o[1,2]/sqrt(o[1,1]*o[2,2])}
		r = matrix(suppressWarnings(apply(draws, 3, func, sigma)), nrow=1)
	} else {
		func = function(draw, sigma) {cov2cor(draw-sigma)[lower.tri(sigma)]}
		r = suppressWarnings(apply(draws, 3, func, sigma))
	}

	out$rho = apply(r, 1, quantile, c(0.025, 0.975), na.rm=T)
	out$r2 = apply(r^2, 1, quantile, c(0.025, 0.975), na.rm=T)
	if (sign(out$rho[1]) != sign(out$rho[2])) out$r2 = c(0, max(0, out$r2[2]))  # set r2 lower to 0 if rho CI spans 0

	return(out)
}


### MULTIPLE REG ###
# expects omega.x to be invertible
ci.multivariate = function(K, omega, sigma, n.iter=10000) {
	P = dim(omega)[1]
	S = diag(sqrt(diag(omega)))
	corrs = solve(S) %*% omega %*% solve(S)
	corrs[corrs >= 1] = 0.99999; corrs[corrs <= -1] = -0.99999; diag(corrs) = 1
	omega = S %*% corrs %*% S
	fit = omega[-P,P] %*% solve(omega[-P,-P]) %*% omega[-P,P]
	if (fit >= omega[P,P]) omega[P,P] = fit/0.99999
	# increasing omega_Y to fit if r2 > 1; setting r2 slightly below 1 in that case, since otherwise the matrixsampling::rwishart function will fail

	gamma.ss = solve(corrs[-P,-P]) %*% corrs[-P,P]
	r2 = max(0, fit/omega[P,P])

	draws = matrixsampling::rwishart(n.iter, K, Sigma=sigma/K, Theta=omega)
	est = apply(draws, 3, estimate.std, sigma)
	failed = apply(is.na(est), 2, any)
	if (mean(failed) > 0.5) stop("too many missing values in permutation output")
	qq = apply(est, 1, quantile, c(0.025, 0.975), na.rm=T)

	return(list(gamma.lower=qq[1,-P], gamma.upper=qq[2,-P], r2.interval=qq[,P]))
}


estimate.std = function(draw, sigma) {
	P = dim(sigma)[1]
	o = cov2cor(draw-sigma)
	g = try(solve(o[-P,-P]) %*% o[-P,P], silent=T)
	if (class(g)[1] == "try-error") return(rep(NA, P))
	r2 = o[-P,P] %*% g
	return(c(g,r2))
}



### PARTIAL COR ###
# omega and sigma should have outcomes x and y at the end
ci.pcor = function(K, omega, sigma, n.iter=10000) {
	P = dim(omega)[1]; Pw = P - 1; Pz = Pw - 1;
	S = diag(sqrt(diag(omega)))
	corrs = solve(S) %*% omega %*% solve(S)
	corrs[corrs >= 1] = 0.99999; corrs[corrs <= -1] = -0.99999; diag(corrs) = 1
	omega = S %*% corrs %*% S

	fit.x = omega[1:Pz,Pw] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,Pw]
	fit.y = omega[1:Pz,P] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P]
	if (omega[Pw,Pw] <= fit.x) omega[Pw,Pw] = fit.x / 0.99999 #scale var up to have R2 slightly below 1
	fit.xy = omega[1:Pw,P] %*% solve(omega[1:Pw,1:Pw]) %*% omega[1:Pw,P]
	if (omega[P,P] <= fit.xy) omega[P,P] = fit.xy / 0.99999 #scale var up to have R2 slightly below 1

	p.cov = omega[Pw,P] - t(omega[1:Pz,Pw]) %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P]
	p.vars = c(omega[Pw,Pw] - fit.x, omega[P,P] - fit.y)

	draws = matrixsampling::rwishart(n.iter, K, Sigma=sigma/K, Theta=omega)
	est = apply(draws, 3, estimate.pcor, sigma)
	if (mean(is.na(est)) > 0.5) stop("too many missing values in permutation output")
	out = quantile(est, c(0.025, 0.975), na.rm=T)
	return(out)
}


estimate.pcor = function(draw, sigma) {
	i.y = dim(sigma)[1]; i.x = i.y - 1; i.z = 1:(i.x-1)
	o = draw - sigma

	inv = try(solve(o[i.z,i.z]), silent=T)
	if (class(inv)[1] == "try-error") return(NA)

	cov = o[i.x,i.y] - t(o[i.z,i.x]) %*% inv %*% o[i.z,i.y]
	vars = c(
		o[i.x,i.x] - t(o[i.z,i.x]) %*% inv %*% o[i.z,i.x],
		o[i.y,i.y] - t(o[i.z,i.y]) %*% inv %*% o[i.z,i.y]
	)
	if (any(vars <= 0)) return(NA)
	return(cov / sqrt(prod(vars)))
}

