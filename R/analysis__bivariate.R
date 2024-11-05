BivariateAnalysis = R6::R6Class("BivariateAnalysis",
	inherit = AnalysisProcessor,
	private = list(
		fit.model = function(targets, phenotypes) {
			rg = self$get.rg(pheno.row=targets, pheno.col=phenotypes, rg.limit=private$settings$rg.limit)
			results = data.frame(phenotype1=rownames(rg)[row(rg)], phenotype2=colnames(rg)[col(rg)], stringsAsFactors=F)
			results$rho = as.numeric(rg)
			results$r2 = results$rho^2

			results = results[order(match(results$phenotype1, targets), match(results$phenotype2, phenotypes)),]
			drop = results$phenotype1 == results$phenotype2 | duplicated(apply(results[,1:2], 1, function(x) {paste(sort(x), collapse=" ")}))
			results = results[!drop,]; rownames(results) = NULL

			if (private$settings$add.ci) results[c("rho.lower", "rho.upper", "r2.lower", "r2.upper")] = NA
			if (private$settings$add.pval) results$p.value = NA

			omega = private$locus.data$get.omega(scale="delta"); sigma = private$locus.data$get.sigma(); K = private$locus.data$no.components()
			for (i in which(!is.na(results$rho))) {
				curr.pheno = as.character(results[i,c("phenotype1", "phenotype2")])

				if (private$settings$add.ci) {
					res = tryCatch(ci.bivariate(K=K, omega=omega[curr.pheno,curr.pheno], sigma=sigma[curr.pheno,curr.pheno]),
												 error = function(err) {private$warning("unable to compute confidence intervals for ", quote.vector(curr.pheno, sep=" ~ "), "; ", err)})
					if (!is.null(res)) {
						results[i,c("rho.lower", "rho.upper")] = res$rho
						results[i,c("r2.lower", "r2.upper")] = res$r2
					}
				}

				if (private$settings$add.pval) {
					res = tryCatch(integral.p(bivariate.integral, K=K, omega=omega[curr.pheno,curr.pheno], sigma=sigma[curr.pheno,curr.pheno], adap.thresh=private$settings$adaptive.threshold),
												 error = function(err) {private$warning("unable to compute p-value for ", quote.vector(curr.pheno, sep=" ~ "), "; ", err)})
					if (!is.null(res)) results$p.value[i] = res
				}
			}

			private$results = private$format.results(results, rho=1, r2=c(0,1), p.value=Inf)
		}
	),
	public = list(
		initialize = function(locus.data, ...) {
			super$initialize(locus.data, ...,
				default.settings=list(
					rg.limit = globals$rg.invalid.bounds
				)
			)
		},

		compute = function(phenotypes=NULL, targets=NULL) {
			targets = private$check.phenotypes(targets)
			phenotypes = private$check.status(private$check.phenotypes(phenotypes), check.univ=T)$available
			if (length(phenotypes) == 0) private$failure("no valid phenotypes to analyse")

			if (length(targets) > 0) {
				status = private$check.status(targets, check.univ=T)
				if (length(status$available) == 0) private$failure("unable to analyse any of target phenotype(s)", lines.list=status$failed)
				else if (length(status$failed) > 0) private$warning("unable to analyse following target phenotype(s)", lines.list=status$failed)
				targets = status$available
			} else {
				if (is.null(targets)) targets = phenotypes
				else private$failure("no valid target phenotypes to analyse")
			}

			private$results = private$fit.model(targets, phenotypes)
			return(private$results)
		}
	)
)

