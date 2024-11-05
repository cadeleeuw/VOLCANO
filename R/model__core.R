

# temporary wrapper function (input ProcessedLocus object)
compute.rg = function(data) {
	out = NULL
	if (data$no.pheno("positive") > 1) {
		model = LocusModel$new(data)
		phenotypes = model$phenotypes(); no.pheno = length(phenotypes)
		omega = model$get.omega(); rg = model$get.rg()

		for (i in 1:(no.pheno-1)) {
			for (j in (i+1):no.pheno) {
				if (!grepl("unavailable", model$get.type(phenotypes[c(i,j)]))) {
					V = model$get.variance(phenotypes[c(i,j)])
					p.value = pnorm(-abs(omega[i,j] / sqrt(V)))*2
					out = rbind(out, data.frame(
						locus=data$get.locus()$get.name(),
						pheno1=phenotypes[i], pheno2=phenotypes[j],
						omega.var1=omega[i,i], omega.var2=omega[j,j],
						omega.cov=omega[i,j], samp.variance=V, rg = rg[i,j],
						p.value = p.value, stringsAsFactors=F
					))
				}
			}
		}
	}
	return(out)
}


## given a ProcessedLocus object, exposes cov(G) estimates and related moments, as well as distribution generators
## input ProcessedLocus object is stripped of all unusable phenotypes upon initialization
LocusModel = R6::R6Class("LocusModel",
	private = list(
		ld.block = NULL, ## global LD block for locus
		marginal = list(), ## named list of MarginalEstimates objects
		pair.type = matrix(), ## named symmetric matrix specifying type for each phenotype pair, currently: "harmonized", "truncated", "unavailable: [reason]"
		estimates = list(sigma = NULL, omega = NULL),
		cache = list(component.products=list()),

		pair.id = function(phenotypes) {return(paste0(sort(phenotypes), collapse=" - "))},

		## product = M = cov(W1,W2); trace = tr(t(M) %*% M)
		component.product = function(phenotypes, compute.trace=F) {
			marginal = private$marginal[phenotypes]
			id = private$pair.id(phenotypes)

			if (!is.null(private$cache$component.product[[id]])) {
				cross = private$cache$component.product[[id]]
				if (nrow(cross$product) != marginal[[1]]$no.components()) cross$product = t(cross$product)
			} else {
				cross = marginal[[1]]$ld.block()$component.correlations(marginal[[2]]$ld.block())
				private$cache$component.product[[id]] = list(product = cross)
			}

			if (compute.trace && is.null(cross$trace)) {
				cross$trace = sum(diag(cross$product %*% t(cross$product)))
				private$cache$component.product[[id]]$trace = cross$trace
			}

			return(cross)
		},

		## match partial names, check and throw error for ambiguous and unknown, remove duplicates/expand wildcards
		check.phenotypes = function(phenotypes) {
			if (length(phenotypes) > 0) return(check.phenotypes(phenotypes, names(private$marginal), label="locus model"))
			else if (is.null(phenotypes)) return(NULL)
			else return(character(0))
		},

		check.pair = function(ph1, ph2) {
			phenotypes = private$check.phenotypes(c(ph1, ph2))
			if (length(phenotypes) != 2) input.error(ifelse(length(phenotypes) > 2, "more", "fewer"), " than two phenotypes specified")
			return(phenotypes)
		},

		process = function(data) {
			private$ld.block = data$get.ld()
			private$marginal = data$get.estimates()

			sigma.scale = sqrt(sapply(private$marginal, function(est) {est$get("sigma")}))
			sigma.corr = data$get.correlations()

			private$estimates$sigma = sweep(sweep(sigma.corr, 1, sigma.scale, FUN="*"), 1, sigma.scale, FUN="*")
			private$estimates$omega = diag.matrix(sapply(private$marginal, function(est) {est$get("omega")}), add.names=names(private$marginal))

			no.pheno = length(private$marginal)
			private$pair.type = diag.matrix(rep(NA, no.pheno), add.names=names(private$marginal))
			if (no.pheno > 1) {
				for (i in 1:(no.pheno-1)) {
					for (j in (i+1):no.pheno) {
						res = private$process.pair(names(private$marginal)[c(i,j)])
						private$estimates$omega[i,j] = private$estimates$omega[j,i] = res$omega
						private$pair.type[i,j] = private$pair.type[j,i] = res$type
					}
				}
			}
		},

		## determines pair type and computes genetic covariance for pair of input phenotypes
		process.pair = function(phenotypes) {
			marginal = private$marginal[phenotypes]
			sigma = private$estimates$sigma[phenotypes, phenotypes]
			delta = lapply(marginal, function(est) {est$get("delta")})
			K = sapply(marginal, function(est) {est$no.components()})

			if (is.na(sigma[1,2])) return(list(omega=NA, type="unavailable: sample overlap is unknown"))

			global.ld = sapply(marginal, function(est) {est$ld.block()$equals(private$ld.block)})
			if (all(global.ld) || marginal[[1]]$ld.block()$equals(marginal[[2]]$ld.block())) { ## same SNPs and PCs
				omega = sum(delta[[1]] * delta[[2]]) - K[1] * sigma[1,2]
				return(list(omega=omega, type="harmonized"))
			} else if (any(global.ld)) {
				if (sigma[1,2] == 0) {
					cross = private$component.product(phenotypes)
					omega = t(delta[[1]]) %*% cross %*% delta[[2]]
					return(list(omega=omega, type="truncated"))
				} else return(list(omega=NA, type="unavailable: truncated phenotype with non-zero sample overlap"))
			} else return(list(omega=NA, type="unavailable: both phenotypes are truncated"))
		},

		compute.variance = function(phenotypes, omega.null=0) {
			marginal = private$marginal[phenotypes]
			sigma = private$estimates$sigma[phenotypes, phenotypes]
			omega = private$estimates$omega[phenotypes, phenotypes]
			delta = lapply(marginal, function(est) {est$get("delta")})
			K = sapply(marginal, function(est) {est$no.components()})

			type = self$get.type(phenotypes)
			if (type == "harmonized") {
				var.base = K[1]*sigma[1,1]*sigma[2,2] + sigma[1,1]*omega[2,2] + sigma[2,2]*omega[1,1]
				var.overlap = 2*omega.null*sigma[1,2] + K[1]*sigma[1,2]^2
				return(as.numeric(var.base + var.overlap))
			} else if (type == "truncated") {
				cross = private$component.product(phenotypes, compute.trace=T)
				var.proj = c(
					sum((t(cross$product) %*% delta[[1]])^2) - sigma[1,1] * cross$trace,
					sum((cross$product %*% delta[[2]])^2) - sigma[2,2] * cross$trace
				)
				var.proj[var.proj < 0] = 0
				var.base = cross$trace*sigma[1,1]*sigma[2,2] + sigma[1,1]*var.proj[2] + sigma[2,2]*var.proj[1]
				return(var.base)
			} else return(NA)
		}
	),
	public = list(
		initialize = function(data) {
			check.types(data="ProcessedLocus", .class=class(self))
			private$process(data$filter("available"))
		},

		phenotypes = function() {return(names(private$marginal))},

		get.type = function(ph1, ph2=NULL) {
			phenotypes = private$check.pair(ph1, ph2)
			return(private$pair.type[phenotypes[1],phenotypes[2]])
		},

		get.omega = function(phenotypes=NULL) {
			if (!is.null(phenotypes)) {
				phenotypes = private$check.phenotypes(phenotypes)
				return(private$estimates$omega[phenotypes,phenotypes])
			} else return(private$estimates$omega)
		},

		get.rg = function(phenotypes=NULL) {return(cov2cor(self$get.omega(phenotypes)))},

		get.variance = function(ph1, ph2=NULL, omega.null=0) {return(private$compute.variance(private$check.pair(ph1, ph2), omega.null))}
	)
)







