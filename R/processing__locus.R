## combined input data for specific locus
LocusData = R6::R6Class("LocusData",
	private = list(
		settings = list(minimum.snps = globals$locus.minimum.snps, minimum.components = globals$locus.minimum.components),
		locus = NULL, ## LocusDefinition object
		ld.block = NULL, ## LDblock object
		sum.stats = list(), ## named list of SumstatBlock object per phenotype
		sampling.correlation = NULL, ## row/col named matrix of sampling correlations
		harmonized = list(), ## named list indicating whether sumstats for a given phenotype are aligned to the global ld.block stored above

		subset.internal = function(phenotypes) {
			private$sum.stats = private$sum.stats[phenotypes]
			private$sampling.correlation = matrix.subset(private$sampling.correlation, cols=phenotypes, rows=phenotypes)
			private$harmonized = private$harmonized[phenotypes]
			return(self)
		}
	),
	public = list(
		initialize = function(locus, ld.block, sum.stats, sampling.correlation, settings.override=list()) {
			if (any(duplicated(names(sum.stats)))) input.error("input to LocusData object contains duplicate phenotypes")
			for (s in names(settings.override)) {if (s %in% names(private$settings)) private$settings[[s]] = settings.override[[s]]}
			private$locus = locus; private$ld.block = ld.block; private$sum.stats = sum.stats

			failed = c()
			for (ph in names(private$sum.stats)) {
				private$harmonized[[ph]] = private$ld.block$equals(private$sum.stats[[ph]]$get.ld())
				if (!private$harmonized[[ph]] && !private$ld.block$get.snp.index()$contains(private$sum.stats[[ph]]$internal.id())) failed = c(failed, ph)
			}
			if (length(failed) > 0) input.error("invalid input to LocusData object, summary statistics for %phenotype% ", items.and=failed, " are not fully contained in provided LD block")

			check.covariance(sampling.correlation, allow.NA=T, label="sampling correlation matrix input to LocusData object")
			if (is.null(colnames(sampling.correlation)) || is.null(rownames(sampling.correlation)) || ncol(sampling.correlation) != length(sum.stats) || !all(colnames(sampling.correlation) == names(sum.stats))) input.error("sampling correlation matrix input to LocusData object does not match input summary statistics")
			private$sampling.correlation = sampling.correlation
		},

		print = function(...) {
			printer = ObjectPrinter$new("LocusData")
			printer$add.parameter("locus", private$locus$get.name())$add.parameter("no. snps", self$no.snps())$add.parameter("no. components", self$no.components())
			printer$add.parameter("phenotypes", length(private$sum.stats))
			for (i in 1:length(private$sum.stats)) printer$add.line(names(private$sum.stats)[[i]], indent=2)
			cat(printer$to.string())
		},

		process = function() {return(ProcessedLocus$new(self, private$settings))},

		phenotypes = function() {return(names(private$sum.stats))},

		no.pheno = function() {return(length(private$sum.stats))},
		no.snps = function(total=T) {return(if (!total) unlist(lapply(private$sum.stats, function(s) {s$no.snps()})) else private$ld.block$no.snps())},
		no.components = function(total=T) {return(if (!total) unlist(lapply(private$sum.stats, function(s) {s$no.components()})) else private$ld.block$no.components())},

		get.locus = function() {return(private$locus)},
		get.ld = function() {return(private$ld.block)},
		get.sumstats = function() {return(private$sum.stats)},
		get.correlations = function() {return(private$sampling.correlation)},

		is.harmonized = function(total=T) {return(if (!total) unlist(private$harmonized) else all(unlist(private$harmonized)))},

		## keep only specified phenotypes if include=T, discard specified phenotypes otherwise
		subset = function(phenotypes, include=T) {
			if (length(phenotypes) > 0) {
				phenotypes = check.phenotypes(phenotypes, self$phenotypes(), label="locus data")
				if (!include) phenotypes = self$phenotypes()[!(self$phenotypes() %in% phenotypes)]
			}
			if (length(phenotypes) == 0) input.error("cannot create data subset, phenotype selection is empty")
			return(self$clone()$.__enclos_env__$private$subset.internal(phenotypes))
		}
	)
)




## processing is considered to have failed for a phenotype if an error occurred during processing or the sigma estimate is negative, but not if the omega variance estimate is negative

ProcessedLocus = R6::R6Class("ProcessedLocus",
	private = list(
		settings = list(minimum.snps = globals$locus.minimum.snps, minimum.components = globals$locus.minimum.components),
		phenotype.info = list(
			info = list(),	## named list of PhenotypeInfo objects
			status = list() ## named list with status (failed/negative/available) per phenotype
		),
		data = list(
			locus = NULL,  ## LocusDefinition object
			ld.block = NULL,  ## LDblock object
			correlations = NULL  ## sampling correlation matrix, subset to non-failed phenotypes
		),
		estimates = list(), ## named list of MarginalEstimates object per phenotype (excluding failed)

		warning = function(..., phenotype=NULL) {throw.warning(..., .label=list(phenotype=phenotype, "locus ID"=private$data$locus$get.name()))},
		failure = function(..., phenotype=NULL) {data.error(..., .label=list(phenotype=phenotype, "locus ID"=private$data$locus$get.name()))},

		## match partial names, check and throw error for ambiguous and unknown, remove duplicates/expand wildcards
		check.phenotypes = function(phenotypes) {
			if (length(phenotypes) > 0) return(check.phenotypes(phenotypes, names(private$phenotype.info$status), label="processed locus"))
			else if (is.null(phenotypes)) return(NULL)
			else return(character(0))
		},

		## NB: currently using trim.components on component drop, which changes shared ld.block object
		process = function(locus.data) {
			phenotypes = locus.data$phenotypes(); estimates = list()
			status = as.list(rep("available", length(phenotypes))); names(status) = phenotypes
			while (TRUE) { ## rerun estimators for all (non-failed) phenotypes until all completed or failed
				for (ph in phenotypes[unlist(status) != "failed"]) {
					curr.stats = locus.data$get.sumstats()[[ph]]
					EstimatorType = if (curr.stats$is.binary()) LogisticBinaryMarginal else ContinuousMarginal
					res = tryCatch(
						EstimatorType$new(curr.stats),
						dropped.components = function(err) {return(err)},
						data.error = function(err) {private$warning(err, phenotype=ph); return(err)},
						error = function(err) {private$failure(err, phenotype=ph)}
					)
					if (is.error(res)) {
						if (is.error(res, "dropped.components")) {   ## break out of for-loop, goes to next iteration of while loop and reruns all (non-failed) phenotypes
							ld.block = locus.data$get.ld()$trim.components(drop=res$get.data()$dropped) ## NB: currently, this modified the ld.block IN PLACE
							if (ld.block$no.components() < private$settings$minimum.components) private$failure("fewer than ", private$settings$minimum.components, " components in locus after pruning unstable components")
							break
						}
						if (is.error(res, "failure")) status[[ph]] = "failed"
					}	else estimates[[ph]] = res
				}
				break  ## all phenotypes processed or failed, break out of while loop
			}

			estimates = estimates[phenotypes[unlist(status) != "failed"]]
			neg.sigma = sapply(estimates, function(est) {est$get("sigma") <= 0})
			if (any(neg.sigma)) private$warning("%phenotype% ", items.and=names(estimates)[neg.sigma], " %has% a negative or zero sampling variance")
			status[names(status) %in% names(estimates)[neg.sigma]] = "failed"

			neg.omega = !neg.sigma & sapply(estimates, function(est) {est$get("omega") <= 0})
			if (any(neg.omega)) private$warning("%phenotype% ", items.and=names(estimates)[neg.omega], " %has% a negative or zero genetic variance estimate")
			status[names(status) %in% names(estimates)[neg.omega]] = "negative"

			estimates = estimates[phenotypes[unlist(status) != "failed"]]
			if (length(estimates) == 0) private$failure("processing has failed for all phenotypes")

			private$estimates = estimates
			private$phenotype.info$status = status
			private$data$correlations = locus.data$subset(names(estimates))$get.correlations()
		},

		subset.internal = function(phenotypes) {
			for (param in names(private$phenotype.info)) private$phenotype.info[[param]] = private$phenotype.info[[param]][phenotypes]

			pheno.processed = phenotypes[phenotypes %in% names(private$estimates)]
			private$estimates = private$estimates[pheno.processed]
			private$data$correlations = matrix.subset(private$data$correlations, cols=pheno.processed, rows=pheno.processed)
			return(self)
		}
	),
	public = list(
		initialize = function(locus.data, settings.override=list()) {
			check.types(locus.data="LocusData", .class=class(self))
			for (s in names(settings.override)) {if (s %in% names(private$settings)) private$settings[[s]] = settings.override[[s]]}

			private$phenotype.info$info = lapply(locus.data$get.sumstats(), function(ss) {ss$phenotype.info()})
			private$data$locus = locus.data$get.locus()
			private$data$ld.block = locus.data$get.ld()

			private$process(locus.data)
		},

		print = function(...) {
			summary = self$summarize()
			printer = ObjectPrinter$new("ProcessedLocus")$parameter.list(snps=private$data$ld.block$no.snps(), components=private$data$ld.block$no.components())
			printer$add.parameter("phenotypes", length(private$phenotype.info$status))
			if (!is.null(summary) && nrow(summary) > 0) {
				for (i in 1:nrow(summary)) {
					phenotype = summary$phenotype[i]
					type = summary$type[i]
					status = if (summary$status[i] != "available") toupper(summary$status[i])
					printer$add.line(phenotype, parentheses(type), parentheses(status), indent=2)
				}
			}
			cat(printer$to.string())
		},

		## phenotypes are subdivided into available/positive, negative (genetic variance) and failed; processed = available/positive + negative)
		## positive/negative refers to the genetic variance estimate; negative sigma estimates are considered 'failed'
		phenotypes = function(include=c("all", "processed", "available", "positive", "failed")) {
			include = match.arg(include); phenotypes = names(private$phenotype.info$status)
			if (length(phenotypes) > 0) {
				if (include == "all") return(phenotypes)
				if (include == "processed") return(phenotypes[unlist(private$phenotype.info$status) != "failed"])
				if (include == "positive" || include == "available") return(phenotypes[unlist(private$phenotype.info$status) == "available"])
				if (include == "failed") return(phenotypes[unlist(private$phenotype.info$status) == "failed"])
			} else return(character(0))
		},

		## return list of matched phenotype names (check validity, convert partial names to full, remove duplicates, expand wildcards)
		match.phenotypes = function(phenotypes) {return(if (length(phenotypes) > 0) private$check.phenotypes(phenotypes) else character(0))},

		no.pheno = function(include=c("all", "processed", "available", "positive", "failed")) {return(length(self$phenotypes(include)))},
		no.snps = function() {return(private$data$ld.block$no.snps())},
		no.components = function() {return(private$data$ld.block$no.components())},

		get.locus = function() {return(private$data$locus)},
		get.ld = function() {return(private$data$ld.block)},

		get.correlations = function(incl.negative=F) {
			if (!incl.negative) {
				keep = self$phenotypes("positive"); if (is.null(keep)) keep = character(0)
				return(matrix.subset(private$data$correlations, cols=keep, rows=keep))
			} else return(private$data$correlations)
		},
		get.estimates = function(incl.negative=F) {return(if (!incl.negative) private$estimates[self$phenotypes("positive")] else private$estimates)},

		summarize = function() {
			if (length(private$phenotype.info$status) > 0) {
				summary = data.frame(
					locus = private$data$locus$get.name(),
					phenotype = names(private$phenotype.info$status),
					type = unlist(lapply(private$phenotype.info$info, function(pi) {pi$get.traits(sep=";")})),
					status = unlist(private$phenotype.info$status),
					no.snps = NA, no.components = NA,
					h2.observed=NA, h2.latent=NA, p.univariate=NA,
					stringsAsFactors=F
				)

				index = match(names(private$estimates), summary$phenotype)
				summary$h2.observed[index] = unlist(lapply(private$estimates, function(est) {est$get("h2")}))
				summary$p.univariate[index] = unlist(lapply(private$estimates, function(est) {est$get("p.value")}))
				summary$no.snps[index] = unlist(lapply(private$estimates, function(est) {est$no.snps()}))
				summary$no.components[index] = unlist(lapply(private$estimates, function(est) {est$no.components()}))

				if (any(sapply(private$phenotype.info$info, function(pi) {pi$is("binary")}))) summary$h2.latent[index] = unlist(lapply(private$estimates, function(est) {est$get("h2.latent")}))
				else summary$h2.latent = NULL

				rownames(summary) = NULL
				return(summary)
			} else return(NULL)
		},

		## keep only specified phenotypes if include=T, discard specified phenotypes otherwise
		## if phenotypes=NULL, no subsetting is performed and a full copy is returned
		subset = function(phenotypes, include=T) {
			if (!is.null(phenotypes)) {
				phenotypes = private$check.phenotypes(phenotypes)
				if (!include) phenotypes = names(private$phenotype.info$status)[!(names(private$phenotype.info$status) %in% phenotypes)]
				if (length(phenotypes) == 0) private$failure("no phenotypes remaining after subsetting")

				return(self$clone()$.__enclos_env__$private$subset.internal(phenotypes))
			} else return(self$clone())
		},

		filter = function(include=c("available", "positive", "processed", "failed")) {return(self$subset(self$phenotypes(match.arg(include)), include=T))},
		filter.h2 = function(min.h2) {
			summary = self$summarize(); phenotypes = c()
			if (!is.null(summary)) phenotypes = summary$phenotype[!is.na(summary$h2) & summary$h2 >= min.h2]
			return(self$subset(phenotypes, include=T))
		},
		filter.pval = function(max.pval) {
			summary = self$summarize(); phenotypes = c()
			if (!is.null(summary)) phenotypes = summary$phenotype[!is.na(summary$p.univariate) & summary$p.univariate <= max.pval]
			return(self$subset(phenotypes, include=T))
		}
	)
)






		## weights argument can be "ivw" or vector of numeric weights; all weights are set to 1 if
		# make.composite = function(name, phenotypes, weights=NULL, drop.components=F) {
		# 	if (name %in% private$phenotype.info) private$error("phenotype name '", name, "' already in use")
		#
		# 	phenotypes = private$check.phenotypes(phenotypes)
		# 	if (length(phenotypes) < 2) private$error("fewer than two phenotypes specified")
		# 	failed = phenotypes %in% private$estimates$failed
		# 	if (any(failed)) private$failure("cannot create composite, locus processing failed for component %phenotype% ", items.and=phenotypes[failed])
		# 	if (any(sapply(private$pheno.info[phenotypes], function(pi) {pi$is("composite")}))) private$error("existing composite phenotypes cannot be specified as input for new composite phenotype")
		#
		# 	if (!is.null(weights)) {
		# 		if (length(weights) > 1) {
		# 			if (length(weights) != length(phenotypes)) private$error("invalid 'weights' argument, length does not match number of phenotypes specified")
		# 			if (!is.numeric(weights)) private$error("invalid 'weights' argument, weights are not numeric")
		# 		} else {
		# 			if (weights[1] == "ivw") weights = 1/diag(private$estimates$sigma[phenotypes,phenotypes])
		# 			else private$error("invalid 'weights' argument")
		# 		}
		# 	} else weights = rep(1, length(phenotypes))
		#
		# 	estimate = CompositeEstimator$new(self$subset(phenotypes), weights)
		# 	if (estimate$get.sigma() <= 0) private$failure("sampling variance of composite phenotype '", name, "' is negative")
		#
		# 	private$phenotype.info = c(private$phenotype.info, name)
		# 	private$pheno.info[[name]] = CompositePhenotype$new(name, private$pheno.info[phenotypes])
		#
		# 	corrs = private$sampling.correlation[,phenotypes] %*% weights / sqrt(sum(weights^2))
		# 	private$sampling.correlation = matrix.expand(private$sampling.correlation, name, c(corrs, 1))
		#
		# 	private$estimates$marginal[[name]] = estimate
		# 	private$process.marginal()
		#
		# 	if (drop.components) private$subset.internal(private$phenotype.info[!(private$phenotype.info %in% phenotypes)])
		# 	invisible(self)
		# }

		## data is scaled to a phenotypic variance of one
#		sum.phenotypes = function(name, phenotypes, weights=NULL, drop.components=F) {return(self$clone()$.__enclos_env__$private$make.composite(name, phenotypes, weights=weights, drop.components=drop.components))},
#		meta.analyse = function(name, phenotypes, drop.components=T) {return(self$clone()$.__enclos_env__$private$make.composite(name, phenotypes, weights="ivw", drop.components=drop.components))}






