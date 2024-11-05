# This function reconstructs the b0 and b1 for the logistic regression of a single SNP from its test statistic
find.beta = function(x, N1, N.orig, stat, tolerance=1e-5, reduction=0.25) {
	tbl.x = table(x); tbl.x = tbl.x / sum(tbl.x) * N.orig;
	val.x = as.numeric(names(tbl.x))
	n.x = as.numeric(tbl.x)
	x = cbind(1, val.x)

	make.b1 = function(curr, value, step, N1, stat, x, n.x, reduction, tolerance) {
		curr$b1 = value
		curr$b0 = find.b0(curr$b0, curr$b1, N1, n.x, x[,2], step=abs(step), reduction=reduction, tolerance=max(abs(step), tolerance))
		curr$err = stat.error(curr$b0, curr$b1, stat, x, n.x)
		return(curr)
	}

	update.b1 = function(curr, step, N1, stat, x, n.x, reduction, tolerance) {
		return(make.b1(curr, curr$b1 + step, step, N1, stat, x, n.x, reduction, tolerance))
	}

	step = 0.01 * sign(stat); try.reverse = F
	curr = make.b1(list(b0 = log(N1/(N.orig-N1))), 0, step, N1, stat, x, n.x, reduction, tolerance)
	while (abs(step) > tolerance/10 && curr$err > tolerance) {
		prop = update.b1(curr, step, N1, stat, x, n.x, reduction, tolerance)
		if (prop$err < curr$err) {
			curr = prop
			try.reverse = F
		} else {
			if (try.reverse) {
				step = -step
				try.reverse = F
			} else {
				step = step * reduction
				try.reverse = (curr$b1 != 0)
			}
		}
	}

	return(c(curr$b0, curr$b1))
}

# finds the b0 corresponding to the given b0, minimizing the difference between sum of model-fitted probabilities and sum of sample probabilities (ie. sum.error())
find.b0 = function(b0, b1, N1, n.x, val.x, step=0.01, reduction=0.25, tolerance=0.01) {
	err = sum.error(b0, b1, N1, n.x, val.x)
	while (step != 0 && abs(err) > tolerance) {
		prop.b0 = b0 + -step*sign(err)
		prop.err = sum.error(prop.b0, b1, N1, n.x, val.x)
		if (abs(prop.err) < abs(err)) {
			b0 = prop.b0
			err = prop.err
		} else {
			step = step * reduction
		}
	}
	return(b0)
}

stat.error = function(b0, b1, stat, x, n.x) {
	mu = 1/(1+exp(-(b0 + b1*x[,2])))
	xsx = t(x) %*% diag(mu*(1-mu)*n.x) %*% x
	var = xsx[1] / (xsx[1]*xsx[4] - xsx[2]^2)
	return(ifelse(var > 0, abs(b1/sqrt(var)-stat)/abs(stat), Inf))
}

sum.error = function(b0, b1, N1, n.x, val.x) {
	mu = 1/(1+exp(-(b0 + b1*val.x)))
	return((sum(mu*n.x) - N1)/N1)
}
