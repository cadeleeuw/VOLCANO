## SumstatBlock objects contain summary statistics and LD reference, computes relevant summary statistics
## UnivariateEstimates wraps around SumstatBlock and fits univariate models


##############################

## interface for block summary statistics data
SumstatBlock = R6::R6Class("SumstatBlock",
	private = list(
		pheno.info = NULL, ## Phenotype object
		ld.block = NULL, ## LDblock object
		snp.data = NULL,	## data.frame with processed summary statistics, harmonized with ld.block
		computed = list(), ## list containing already computed internal statistics (to avoid recomputation on successive calls)

		fatal = function(...) {fatal.error(..., .label=list(phenotype=private$pheno.info$get.name()))},

		resolve.level = function(name, output.type=c("automatic", "global", "snp")) {
			output.type = match.arg(output.type)
			if (name %in% names(private$snp.data)) {
				if (output.type == "global") {
					if (!(name %in% private$computed)) private$computed[[name]] = mean(private$snp.data[[name]], na.rm=T)
					return(private$computed$sample.size)
				} else return(private$snp.data[[name]])
			} else {
				if (output.type == "snp") return(rep(private$pheno.info$get.global(name), nrow(private$snp.data)))
				else return(private$pheno.info$get.global(name))
			}
		},

		compute.xty = function() {undefined.error("compute.xty", class(self)[1])},
		compute.covar = function() {undefined.error("compute.covar", class(self)[1])}
	),
	public = list(
		initialize = function(ld.block, snp.data, pheno.info) {
			private$ld.block = ld.block; private$pheno.info = pheno.info
			private$snp.data = tryCatch(ld.block$align.data(snp.data, mode="sort"), error=function(err) {private$fatal("invalid input to SumstatBlock object, ", err)})
		},

		print = function(...) {
			printer = ObjectPrinter$new(class(self)[1])
			printer$parameter.list(phenotype=self$phenotype())$add.parameter("input type", self$data.type(as.string=T))$add.parameter("no. snps", self$no.snps())
			cat(printer$to.string())
		},

		phenotype = function() {return(private$pheno.info$get.name())},
		data.type = function(as.string=F, sep=" ") {return(private$pheno.info$get.traits(as.string=as.string, sep=sep))},
		is.binary = function() {return(private$pheno.info$is("binary"))},
		phenotype.info = function() {return(private$pheno.info)},

		no.snps = function() {return(nrow(private$snp.data))},
		no.components = function() {return(private$ld.block$no.components())},
		pheno.variance = function() {undefined.error("pheno.variance", class(self)[1])},

		snp.id = function() {return(private$ld.block$get.snp.info()$snp.id)},
		internal.id = function() {return(private$snp.data$internal.id)},
		alleles = function() {return(private$ld.block$get.snp.info()[,c("allele1", "allele2")])},
		sample.size = function(output.type=c("automatic", "global", "snp")) {return(private$resolve.level("sample.size", output.type))},
		statistic = function() {return(private$snp.data$statistic)},
		p.value = function() {return(2*pnorm(abs(private$snp.data$statistic), lower.tail=F))},
		xty = function() {
			if (!("xty" %in% names(private$computed))) private$compute.xty()
			return(private$computed$xty)
		},
		covariance = function(standardize=T) {
			if (!("covar.raw" %in% names(private$computed))) private$compute.covar()
			covar = private$computed$covar.raw
			if (standardize && self$pheno.variance() != 1) covar = covar / sqrt(self$pheno.variance())
			return(covar)
		},

		get.ld = function() {return(private$ld.block)},
		get.snp.data = function() {return(cbind(private$ld.block$get.snp.info()[,c("snp.id", "allele1", "allele2")], subset(private$snp.data, select=-internal.id)))},
		get.pheno.info = function() {return(private$pheno.info)}
	)
)

ContinuousSumstatBlock = R6::R6Class("ContinuousSumstatBlock",
	inherit = SumstatBlock,
	private = list(
		compute.xty = function() {private$computed$xty = self$covariance(standardize=F) * (self$sample.size() - 1)},
		compute.covar = function() {
			beta = private$snp.data$statistic / sqrt(private$snp.data$statistic^2 + self$sample.size() - 2)
			if (self$pheno.variance() != 1) beta = beta * sqrt(self$pheno.variance())
			private$computed$covar.raw = beta
		}
	),
	public = list(
		initialize = function(ld.block, snp.data, pheno.info) {
			if (pheno.info$is("binary")) private$fatal("binary phenotype input into ContinuousSumstatBlock object")
			super$initialize(ld.block, snp.data, pheno.info)
		},

		pheno.variance = function() {return(1)}
	)
)


## currently requires raw genotype data from reference to run
BinarySumstatBlock = R6::R6Class("BinarySumstatBlock",
	inherit = SumstatBlock,
	private = list(
		compute.xty = function() {
			genotypes = private$ld.block$get.data()
			sample.size = self$sample.size(output.type="snp"); case.proportion = self$case.proportion(output.type="snp")

			private$computed$xty = rep(0, nrow(private$snp.data))
			for (i in 1:nrow(private$snp.data)) {
				beta = find.beta(genotypes[,i], case.proportion[i]*sample.size[i], sample.size[i], private$snp.data$statistic[i])
				mu = 1 / (1+exp(-(beta[1] + beta[2]*genotypes[,i])))
				private$computed$xty[i] = sum(genotypes[,i] * mu)
			}
		},
		compute.covar = function() {private$computed$covar.raw = self$get.xty() / (self$sample.size() - 1)}
	),
	public = list(
		initialize = function(ld.block, snp.data, pheno.info) {
			if (pheno.info$is("continuous")) private$fatal("non-binary phenotype input into BinarySumstatBlock object")
			super$initialize(ld.block, snp.data, pheno.info)
		},

		pheno.variance = function() {
			sample.size = self$sample.size(output.type="global"); proportion = self$case.proportion(output.type="global")
			return(proportion * (1 - proportion) * sample.size / (sample.size - 1))
		},
		prevalence = function() {return(private$pheno.info$get.global("prevalence"))},
		case.proportion = function(output.type=c("automatic", "global", "snp")) {return(private$resolve.level("case.proportion", output.type))}
	)
)



##############################


## estimates are scaled to a phenotypic variance of one for continuous and binary phenotypes
## composite phenotypes are on an undefined scale, and currently do not define h2 values
MarginalEstimates = R6::R6Class("MarginalEstimates",
	private = list(
		info = list(phenotype = NULL, no.snps = NULL, no.components = NULL, sample.size = NULL),
		estimates = list(delta = NULL, eta = NULL, sigma = NULL, omega = NULL, h2 = NULL, h2.latent = NULL, p.value = NULL),

		failure = function(...) {data.error(..., .label=list(phenotype=private$info$phenotype))},

		compute.chisq = function(delta, sigma) {return(pchisq(ifelse(sigma > 0, sum(delta^2) / sigma, NA), df=nrow(delta), lower.tail=F))},
		compute.F = function(delta, sigma, sample.size, no.components) {
			if (sample.size - no.components - 1 <= 0) input.error("degrees of freedom for F-test are zero or negative")
			return(pf(ifelse(sigma > 0, sum(delta^2) / sigma / no.components, NA), no.components, sample.size - no.components - 1, lower.tail=F))
		}
	),
	public = list(
		phenotype = function() {return(private$info$phenotype)},
		no.snps = function() {return(private$info$no.snps)},
		no.components = function() {return(private$info$no.components)},
		ld.block = function() {undefined.error("ld.block", class(self)[1])},

		## for parameters: delta, sigma, eta, omega, h2, h2.latent, p.value; returns NA if not present
		get = function(parameter) {return(if (!is.null(private$estimates[[parameter]])) private$estimates[[parameter]] else NA)},
		get.estimates = function() {return(private$estimates)},

		is.binary = function() {return(F)},
		is.composite = function() {return(F)}
	)
)


MarginalDataEstimates = R6::R6Class("MarginalDataEstimates",
	inherit = MarginalEstimates,
	private = list(
		sum.stats = NULL, ## SumstatBlock object

		estimate = function() {undefined.error("estimate", class(self)[1])}
	),
	public = list(
		initialize = function(sum.stats) {
			check.types(sum.stats="SumstatBlock", .class=class(self))
			private$sum.stats = sum.stats
			private$info = list(
				phenotype = sum.stats$phenotype(),
				no.snps = sum.stats$no.snps(),
				no.components = sum.stats$no.components(),
				sample.size = sum.stats$sample.size(output.type="global")
			)
			private$estimate()

			for (param in c("h2", "h2.latent")) {
				if (param %in% names(private$estimates) && !is.na(private$estimates[[param]])) {
					if (private$info$no.components / private$info$sample.size > globals$h2.component.ratio) private$estimates[[param]] = NA
					else if (private$estimates[[param]] < 0) private$estimates[[param]] = 0
				}
			}
		},

		ld.block = function() {return(private$sum.stats$get.ld())}
	)
)


ContinuousMarginal = R6::R6Class("ContinuousMarginal",
	inherit = MarginalDataEstimates,
	private = list(
		estimate = function() {
			N = private$info$sample.size; K = private$info$no.components
			delta = t(private$sum.stats$get.ld()$get.inverse.root()) %*% private$sum.stats$covariance(standardize=T) / sqrt(private$sum.stats$pheno.variance())
			eta = (1 - sum(delta^2)) * (N - 1) / (N - K - 1)
			sigma = eta / (N - 1)
			omega = sum(delta^2) - K*sigma
			h2 = 1 - eta  ## will be the same as omega
			p.value = private$compute.F(delta, sigma, N, K)

			private$estimates = list(delta=delta, eta=eta, sigma=sigma, omega=omega, h2=h2, p.value=p.value)
		}
	)
)


BinaryMarginal = R6::R6Class("BinaryMarginal",
	inherit = MarginalDataEstimates,
	public = list(is.binary = function() {return(T)})
)


LogisticBinaryMarginal = R6::R6Class("LogisticBinaryMarginal",
	inherit = BinaryMarginal,
	private = list(
		estimate = function() {no.iter=25; svar.thresh=1.75 #leaving in hard-coded parameters for now, remove later
			ld.block = private$sum.stats$get.ld(); genotypes = private$ld.block$get.data(); components = private$ld.block$get.components()
			N = private$info$sample.size; N.ref = nrow(genotypes); N.case = N * sum.stats$case.proportion(output.type="global")
			K = private$info$no.components; K.snp = private$info.no.snps

			wty0 = N.case * N.ref / N
			wty1 = t(ld.block$get.inverse.root()) %*% private$sum.stats$xty()

			comp.use = 1:K
			while (length(comp.use) > 1) {
				W = cbind(1, components[,comp.use])
				wty = c(wty0, wty1[comp.use])

				beta = c(log(no.cases/(N - no.cases)), rep(0, length(comp.use)))
				for (i in 1:no.iter) {
					mu = 1 / (1+exp(-W%*%beta))
					s = as.numeric(mu * (1-mu))
					wsw.inv = try(solve(t(W) %*% diag(s) %*% W), silent=T)
					if (class(wsw.inv)[1] == "try-error") private$failure("inversion error when fitting multiple logistic regression model")
					beta = beta + wsw.inv %*% (wty - t(W) %*% mu)
				}
				V = diag(wsw.inv)[-1]
				svar.ratio = max(V) / quantile(V, 0.5)	# using ratio of max to median; will normally be very close to 1
				if (svar.ratio < svar.thresh) break 	# if above threshold, drop the PC with highest sampling variance and rerun
				comp.use = comp.use[-which.max(V)]
			}
			if (length(comp.use) < K) private$failure("unstable components in multiple logistic regression model", .sub.type="dropped.components", .data=list(dropped=which(!(1:K %in% comp.use))))

			linear.variance = sum(wty1^2) / ((N - 1) * N.ref/N)^2   ## ie. var(G) under linear regression of Y on W, G = W*delta.linear
			pheno.variance = private$sum.stats$pheno.variance()


			delta = beta[-1,] / sqrt(pheno.variance)
			sigma = mean(V) * N.ref / N / pheno.variance
			omega = sum(delta^2) - K*sigma
			h2 = 1 - (1 - linear.variance / pheno.variance) * (N - 1) / (N - K - 1)   ## observed h2
			p.value = private$compute.chisq(delta, sigma)
			private$estimates = list(delta=delta, sigma=sigma, omega=omega, h2=h2, p.value=p.value)

			prevalence = private$sum.stats$prevalence()
			if (!is.na(prevalence)) private$estimates$h2.latent = h2 * prevalence/(1-prevalence) / dnorm(qnorm(prevalence))^2
		}
	)
)


## NB: currently does not set h2
# CompositeEstimator = R6::R6Class("CompositeEstimator",
# 	inherit = MarginalEstimates,
# 	private = list(
# 		combine = function(processed.locus, weights) {
# 			if (is.null(weights)) weights = rep(1, processed.locus$no.pheno())
# 			weights = weights / sqrt(sum(weights^2))
#
# 			private$estimates$delta = processed.locus$get.delta() %*% weights
# 			private$estimates$sigma = t(weights) %*% processed.locus$get.sigma() %*% weights
# 		}
# 	),
# 	public = list(
# 		initialize = function(processed.locus, weights=NULL) {
# 			check.types(processed.locus="ProcessedLocus", .class=class(self))
# 			if (processed.locus$no.pheno("failed") > 0) throw.failure("cannot create composite, input contains failed phenotypes")
# 			if (processed.locus$no.pheno() < 2) throw.error("fewer than two phenotypes provided, cannot merge")
# 			if (length(weights) > 0 && length(weights) != processed.locus$no.pheno()) throw.error("invalid weights vector")
#
# 			private$combine(processed.locus, weights)
# 		},
#
# 		is.composite = function() {return(T)},
# 		get.pval = function() {return(private$compute.chisq(private$estimates$delta, private$estimates$sigma))}
# 	)
# )

