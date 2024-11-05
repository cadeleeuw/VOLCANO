RegressionAnalysis = R6::R6Class("RegressionAnalysis",
	inherit = AnalysisProcessor,
	private = list(
		model.name = function(outcome, predictors) {return(paste0(outcome, " | ", paste(predictors, collapse=" + ")))},

		fit.model = function(outcome, predictors) {
			if (length(outcome) != 1 || length(predictors) == 0) private$error("invalid input specification for fit.model()")
			omega = private$locus.data$get.omega(scale="delta"); sigma = private$locus.data$get.sigma()
			K = private$locus.data$no.components(); no.predictors = length(predictors)
			model.name = private$model.name(outcome, predictors)

			coefficients = data.frame(predictor=predictors, gamma.raw=NA, gamma=NA)
			parameters = list(tau.raw=NA, tau=NA, r2=NA)
			if (private$settings$add.ci) {coefficients[c("gamma.lower", "gamma.upper")] = NA; parameters[c("r2.lower", "r2.upper")] = NA}

			if (no.predictors > 1) {
				omega.x = omega[predictors,predictors]; omega.xy = omega[predictors,outcome]; omega.y = omega[outcome,outcome]
				omega.x.inv = tryCatch(eigen.invert(omega.x), error=function(err) {private$failure("predictor genetic covariance matrix is not invertible for model ", model.name, sub.type="inversion")})

				coefficients$gamma.raw = omega.x.inv %*% omega.xy
				coefficients$gamma = coefficients$gamma.raw * sqrt(diag(omega.x)/omega.y)
				parameters$tau.raw = as.numeric(omega.y - t(omega.xy) %*% omega.x.inv %*% omega.xy)

				if (private$settings$add.pval) coefficients$p.value = NA
				invalid = abs(coefficients$gamma) > private$settings$coefficient.limit
				if (!any(invalid)) {
					phenotypes = c(predictors, outcome)

					if (private$settings$add.ci) {
						res = tryCatch(suppressWarnings(ci.multivariate(K=K, omega=omega[phenotypes,phenotypes], sigma=sigma[phenotypes,phenotypes])),
													 error = function(err) {private$warning("unable to compute confidence intervals for model ", model.name, "; ", err)})
						if (!is.null(res)) {
							coefficients$gamma.lower = res$gamma.lower
							coefficients$gamma.upper = res$gamma.upper
							parameters[c("r2.lower", "r2.upper")] = res$r2.interval
						}
					}

					if (private$settings$add.pval) {
						res = tryCatch(suppressWarnings(integral.p(multivariate.integral, K=K, omega=omega[phenotypes,phenotypes], sigma=sigma[phenotypes,phenotypes], adap.thresh=private$settings$adaptive.threshold)),
													 error = function(err)	{private$warning("unable to compute p-values for model ", model.name, "; ", err)})
						if (!is.null(res)) coefficients$p.value = res
					}
				} else {
					private$warning("model ", model.name, " failed, coefficient estimates out of bounds (+/- ", private$settings$coefficient.limit, ")")
					coefficients$gamma[invalid] = NA; parameters$tau.raw = NA
				}

				parameters$tau = parameters$tau.raw / omega.y
				parameters$r2 = 1 - parameters$tau
			} else {
				bivariate = BivariateAnalysis$new(private$locus.data, import.settings=private$settings, rg.limit=private$settings$coefficient.limit)
				results = bivariate$compute(predictors, outcome)

				coefficients$gamma.raw = omega[predictors,outcome] / omega[predictors,predictors]
				coefficients$gamma = results$rho
				parameters$tau.raw = omega[outcome,outcome] - omega[predictors,outcome]^2 / omega[predictors,predictors]
				parameters$tau = 1 - results$r2
				parameters$r2 = results$r2

				if (private$settings$add.ci) {
					coefficients[c("gamma.lower", "gamma.upper")] = results[c("rho.lower", "rho.upper")]
					parameters[c("r2.lower", "r2.upper")] = results[c("r2.lower", "r2.upper")]
				}

				if (private$settings$add.pval) coefficients$p.value = results$p.value
			}

			results = list(
				coefficients = private$format.results(coefficients, gamma=Inf, p.value=Inf),
				parameters = private$format.results(data.frame(parameters), tau=c(0,1), r2=c(0,1))
			)
			return(results)
		}
	),
	public = list(
		initialize = function(locus.data, ...) {
			super$initialize(locus.data, ...,
				default.settings=list(
					coefficient.limit = globals$regression.invalid.bounds,
					filter.predictors = F
				)
			)
		},

		compute = function(outcome, predictors=NULL) {
			outcome = private$check.phenotypes(outcome)
			if (length(outcome) == 1) {
				status = private$check.status(outcome, check.univ=T)
				if (length(status$failed) > 0) private$failure("unable to analyse outcome phenotype '", outcome, "': ", status$failed[[outcome]])
			} else {
				if (length(outcome) == 0) private$error("no outcome phenotype specified")
				else private$error("more than one outcome phenotype specified")
			}

			if (length(predictors) > 0) {
				predictors = private$check.phenotypes(predictors)
				if (outcome %in% predictors) private$error("outcome phenotype cannot be specified as a predictor")
			} else if (is.null(predictors))	predictors = setdiff(private$locus.data$phenotypes(), outcome)

			if (length(predictors) > 0) {
				status = private$check.status(predictors, check.univ=T)
				if (length(status$failed) > 0) {
					if (private$settings$filter.predictors) {
						if (length(status$available) > 0) {
							private$warning("unable to analyse following predictor phenotype(s), reducing model", lines.list=status$failed)
							predictors = status$available
						} else private$failure("no predictor phenotypes remaining after filtering")
					} else private$failure("not all predictor phenotype(s) are available", lines.list=status$failed)
				}
			} else private$error("no predictor phenotypes have been specified")

			model = private$fit.model(outcome, predictors)
			private$results = c(list(outcome=outcome, predictors=predictors), model)
			return(private$results)
		},

		compute.submodels = function(include.marginal=F) {
			if (is.null(private$results)) private$failure("cannot compute submodels, no valid main model has been computed")

			main.model = private$results; predictors = main.model$predictors; no.predictors = length(predictors)
			if (no.predictors == 2 && !include.marginal) private$failure("main model contains only two predictors, no submodels to compute", sub.type="no.submodels")
			if (no.predictors == 1) private$failure("main model contains only one predictor, no submodels to compute", sub.type="no.submodels")

			output = list(outcome=main.model$outcome, predictors=predictors)
			predictor.sets = coefficients = parameters = list(); failed = NULL
			for (no.sub in ifelse(include.marginal, 1, 2):(no.predictors-1)) {
				pattern.matrix = t(combn(1:no.predictors, no.sub))
				for (i in 1:nrow(pattern.matrix)) {
					curr.predictors = predictors[pattern.matrix[i,]]
					curr.name = private$model.name(main.model$outcome, curr.predictors)

					curr.model = tryCatch(private$fit.model(main.model$outcome, curr.predictors),
																error = function(err) {private$warning("could not fit submodel ", curr.name, "; ", err)})
					if (!is.null(curr.model)) {
						model.id = length(predictor.sets) + 1
						predictor.sets[[model.id]] = curr.predictors
						coefficients[[model.id]] = curr.model$coefficients
						parameters[[model.id]] = curr.model$parameters
					} else failed = c(failed, curr.name)
				}
			}

			output$no.models = length(predictor.sets)+1
			if (!is.null(failed)) output$failed = failed

			predictor.sets[[output$no.models]] = predictors
			coefficients[[output$no.models]] = main.model$coefficients
			parameters[[output$no.models]] = main.model$parameters

			output$coefficients = fill.rowmerge(lapply(seq_along(coefficients), function(i) {data.frame(model.id=i, coefficients[[i]])}))
			output$parameters = data.frame(
				model.id = seq_along(parameters),
				no.predictors = sapply(predictor.sets, length),
				predictors = sapply(predictor.sets, paste, collapse=", "),
				fill.rowmerge(parameters),
				stringsAsFactors=F
			)

			return(output)
		}
	)
)




PartialCorrelationAnalysis = R6::R6Class("PartialCorrelationAnalysis",
	inherit = AnalysisProcessor,
	private = list(
		model.name = function(outcomes, predictors) {return(paste0(outcomes[1], " ~ ", outcomes[2], " | ", paste(predictors, collapse=" + ")))},

		fit.model = function(outcomes, predictors) {
			if (length(outcomes) < 2 || length(predictors) == 0) private$error("invalid input specification for fit.model()")

			rg = self$get.rg(pheno.row=outcomes, pheno.col=outcomes, rg.limit=private$settings$rg.limit, verbose=F); tri = upper.tri(rg)
			results = data.frame(phenotype1=rownames(rg)[row(rg)[tri]], phenotype2=colnames(rg)[col(rg)[tri]], r2.marginal1=NA, r2.marginal2=NA, rho=rg[tri], rho.partial=NA, stringsAsFactors=F)
			if (private$settings$add.ci) results[c("rho.partial.lower", "rho.partial.upper")] = NA
			if (private$settings$add.pval) results$p.value = NA
			rownames(results) = NULL

			omega = private$locus.data$get.omega(scale="delta"); sigma = private$locus.data$get.sigma()
			K = private$locus.data$no.components(); no.predictors = length(predictors)
			omega.pred.inv = tryCatch(eigen.invert(omega[predictors,predictors]), error=function(err) {private$failure("genetic covariance matrix of conditioned-on phenotypes is not invertible", .sub.type="inversion")})

			reg.model = RegressionAnalysis$new(private$locus.data, coefficient.limit=Inf, add.pval=F, add.ci=F, signif.decimals=private$settings$signif.decimals)
			reg.results = univ.r2 = partial.var = failed = list()
			for (o in outcomes) {
				reg.results[[o]] = reg.model$compute(o, predictors)
				univ.r2[[o]] = ifelse(!is.null(reg.results[[o]]) && !is.null(reg.results[[o]]$parameters$r2), reg.results[[o]]$parameters$r2, NA)
				if (!is.na(univ.r2[[o]])) {
					if (univ.r2[[o]] <= private$settings$maximum.r2) {
						partial.var[[o]] = omega[o,o] - omega[o,predictors] %*% omega.pred.inv %*% omega[predictors,o]
						if (is.na(partial.var[[o]]) || partial.var[[o]] <= 0) failed[[o]] = "invalid partial variance"
					} else failed[[o]] = paste0("marginal r2 exceeds maximum of ", private$settings$maximum.r2)
				} else failed[[o]] = "unable to compute marginal r2"
			}

			if (length(failed) > 0) {
				if (length(outcomes) - length(failed) < 2) private$failure("marginal model failed for following	outcomes, no correlations to compute", lines.list=failed)
				else private$warning("marginal model failed for following	outcomes", lines.list=failed)
			}

			results$r2.marginal1 = unlist(univ.r2[results$phenotype1])
			results$r2.marginal2 = unlist(univ.r2[results$phenotype2])

			for (i in 1:nrow(results)) {
				curr.outcomes = unlist(results[i,c("phenotype1", "phenotype2")])
				if (!any(curr.outcomes %in% names(failed))) {
					model.name = private$model.name(curr.outcomes, predictors)
					vars = unlist(partial.var[curr.outcomes])
					partial.cor = as.numeric(omega[curr.outcomes[1],curr.outcomes[2]] - t(omega[outcomes[1],predictors]) %*% omega.pred.inv %*% omega[predictors,outcomes[2]]) / sqrt(prod(vars))
					if (abs(partial.cor) <= private$settings$rg.limit) {
						results$rho.partial[i] = partial.cor
		 				phenotypes = c(predictors, curr.outcomes)

						if (private$settings$add.ci) {
							res = tryCatch(suppressWarnings(ci.pcor(K=K, omega=omega[phenotypes,phenotypes], sigma=sigma[phenotypes,phenotypes])),
														 error = function(err) {private$warning("unable to compute confidence intervals for model ", model.name, "; ", err)})
							if (!is.null(res)) results[i,c("rho.partial.lower", "rho.partial.upper")] = res
						}

						if (private$settings$add.pval) {
							res = tryCatch(suppressWarnings(integral.p(pcov.integral, K=K, omega=omega[phenotypes,phenotypes], sigma=sigma[phenotypes,phenotypes], adap.thresh=private$settings$adaptive.threshold)),
														 error = function(err)	{private$warning("unable to compute p-values for model ", model.name, "; ", err)})
							if (!is.null(res)) results$p.value[i] = res
						}
					} else private$warning("partial correlation estimate too far out of bounds (+/-", private$settings$rg.limit, ") for model ", model.name)
				}
			}
			return(private$format.results(data.frame(results), rho.partial=1, p.value=Inf))
		}
	),
	public = list(
		initialize = function(locus.data, ...) {
			super$initialize(locus.data, ...,
				default.settings=list(
					rg.limit = globals$rg.invalid.bounds,
					maximum.r2 = globals$partial.cor.max.r2,
					filter.variables = F
				)
			)
		},

		compute = function(outcomes, predictors=NULL) {
			outcomes = private$check.phenotypes(outcomes)
			if (length(outcomes) >= 2) {
				status = private$check.status(outcomes, check.univ=T)
				if (length(status$failed) > 0) {
					if (length(status$available) >= 2) {
						private$warning("unable to analyse following outcome phenotypes", lines.list=status$failed)
						outcomes = status$available
					} else 	private$failure("unable to analyse following outcome phenotypes, no correlations to compute", lines.list=status$failed)
				}
			} else {
				if (length(outcomes) == 0) private$error("no outcome phenotypes specified")
				else private$error("only one outcome phenotype specified")
			}

			if (length(predictors) > 0) {
				predictors = private$check.phenotypes(predictors)
				if (any(outcomes %in% predictors)) private$failure("outcome phenotypes cannot be specified to be conditioned on")
			} else if (is.null(predictors)) predictors = setdiff(private$locus.data$phenotypes(), outcomes)

			if (length(predictors) > 0) {
				status = private$check.status(predictors, check.univ=T)
				if (length(status$failed) > 0) {
					if (private$settings$filter.variables) {
						if (length(status$available) > 0) {
							private$warning("unable to condition on following phenotype(s), reducing model", lines.list=status$failed)
							predictors = status$available
						} else private$failure("no phenotypes to condition on are remaining after filtering")
					} else private$failure("not all conditioned-on phenotype(s) are available", lines.list=status$failed)
				}
			} else private$error("no phenotypes to condition on have been specified")

			model = private$fit.model(outcomes,  predictors)
			private$results = list(outcomes=outcomes, conditioned=predictors, results=model)
			return(private$results)
		},

		compute.submodels = function() {
			if (is.null(private$results)) private$failure("cannot compute submodels, no valid main model has been computed")
			if (length(private$results$outcomes) != 2) private$error("cannot compute submodels, multiple main models are currently set")

			main.model = private$results; predictors = main.model$conditioned; no.predictors = length(predictors)
			if (no.predictors == 1) private$failure("main model conditions on only one phenotype, no submodels to compute", sub.type="no.submodels")

			output = list(outcomes=main.model$outcome, conditioned=predictors)
			results = failed = NULL
			for (no.sub in 1:(no.predictors-1)) {
				pattern.matrix = t(combn(1:no.predictors, no.sub))
				for (i in 1:nrow(pattern.matrix)) {
					curr.predictors = predictors[pattern.matrix[i,]]
					curr.name = private$model.name(main.model$outcomes, curr.predictors)

					curr.model = tryCatch(private$fit.model(main.model$outcomes, curr.predictors), error = function(err) {private$warning("could not fit submodel ", curr.name, "; ", err)})
					if (!is.null(curr.model)) results = rbind(results, data.frame(no.conditioned=no.sub, conditioned=paste(curr.predictors, collapse=", "), curr.model, stringsAsFactors=F))
					else failed = c(failed, curr.name)
				}
			}

			output$no.models = length(results) + 1
			if (!is.null(failed)) output$failed = failed
			output$results = fill.rowmerge(list(results, data.frame(no.conditioned=length(predictors), conditioned=paste(predictors, collapse=", "), main.model$results, stringsAsFactors=F)))
			output$results = subset(output$results, select=-c(phenotype1, phenotype2))
			return(output)
		}
	)
)

