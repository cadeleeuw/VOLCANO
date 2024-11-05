#' Bivariate local genetic correlation analysis
#'
#' Performs bivariate local genetic correlation analysis between two phenotypes.
#' By default, the bivariate test will be performed for all combinations of phenotypes in the locus,
#' but this can be modified using the 'phenos' and 'target' arguments (see below)
#'
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param phenos Subset of phenotypes to analyse. If NULL, all phenotypes in the locus object will be analysed
#' @param target Target phenotype of interest. If NULL, bivariate correlations between all pairs of phenotypes will be computed;
#' Otherwise, only the relations between the target phenotype and the other phenotypes will be tested.
#' @param adap.thresh The thresholds at which to increase the number of iterations for the p-value generation.
#' Default number of iterations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of iterations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param param.lim The +- threshold at which estimated parameters are considered to be too far out of bounds. If the estimated parameter exceeds this threshold, it is considered unreliable and will be set to NA.
#'
#' @return Data frame with the columns:
#' \itemize{
#'     \item phen1 / phen2 - analysed phenotypes
#'     \item rho - standardised coefficient for the local genetic correlation
#'     \item rho.lower / rho.upper - 95\% confidence intervals for rho
#'     \item r2 - proportion of variance in genetic signal for phen1 explained by phen2 (and vice versa)
#'     \item r2.lower / r2.upper - 95\% confidence intervals for the r2
#'     \item p - simulation p-values for the local genetic correlation
#' }
#' @export



#TODO: update help meta-data
# NB: internal failures are downgraded to printed error message if catch.failures=T (
run.bivar = function(locus.data, phenotypes=NULL, targets=NULL, add.pval=T, add.ci=T, output.univ=F, univ.thresh=1, adap.thresh=globals$adaptive.pval.threshold, param.lim=globals$rg.invalid.bounds, catch.failures=T) {
	analysis = BivariateAnalysis$new(locus.data, univ.threshold=univ.thresh, adaptive.threshold=ifelse(is.null(adap.thresh), 0, adap.thresh), add.pval=add.pval, add.ci=add.ci, rg.limit=param.lim)
	results = catchErrors(analysis$compute(phenotypes, targets), all.failures=catch.failures)

	if (!is.null(results)) {
		if (output.univ) {
			if (is.null(phenotypes)) phenotypes = analysis$phenotypes()
			univariate = analysis$get.univariate(c(targets, phenotypes))
			return(list(univariate=univariate, bivariate=results))
		} else return(results)
	} else return(NULL)
}




#' Local genetic multiple regression analysis
#'
#' Will perform a local genetic multiple regression analysis, which models the genetic signal for a single outcome phenotype of interest using two or more predictor phenotypes.
#' Here, the genetic correlations between all predictors will be accounted for, and their genetic relation with the outcome will be conditioned on one another.
#'
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param target Outcome phenotype of interest (all other phenotypes will be considered predictors)
#' @param phenos Subset of phenotypes to analyse. If NULL, all phenotypes in the locus object will be analysed.
#' @param adap.thresh The thresholds at which to increase the number of iterations for the p-value generation.
#' Default number of iterations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of iterations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param param.lim The +- threshold at which estimated parameters are considered to be too far out of bounds. If the estimated parameter exceeds this threshold, it is considered unreliable and will be set to NA.
#'
#' @return Data frame with the columns:
#' \itemize{
#'     \item predictors / outcome - analysed phenotypes
#'     \item gamma - standardised multiple regression coefficient
#'     \item gamma.lower / gamma.upper - 95\% confidence intervals for gamma
#'     \item r2 - proportion of variance in genetic signal for the outcome phenotype explained by all predictor phenotypes simultaneously
#'     \item r2.lower / r2.upper - 95\% confidence intervals for the r2
#'     \item p - simulation p-values for the gammas
#' }
#' @export
#'
#'
#'



run.multireg = function(locus.data, outcome, predictors=NULL, add.pval=T, add.ci=T, output.marginal=F, output.submodels=F, univ.thresh=1, filter.pred=(univ.thresh < 1), adap.thresh=globals$adaptive.pval.threshold, param.lim=globals$coefficient.invalid.bounds, catch.failures=T) {
	analysis = RegressionAnalysis$new(locus.data, univ.threshold=univ.thresh, adaptive.threshold=ifelse(is.null(adap.thresh), 0, adap.thresh), add.pval=add.pval, add.ci=add.ci, coefficient.limit=param.lim, filter.predictors=filter.pred)
	results = catchErrors(analysis$compute(outcome, predictors), all.failures=catch.failures)

	if (!is.null(results)) {
		results$coefficients = subset(results$coefficients, select=-gamma.raw)
		results$parameters = subset(results$parameters, select=-c(tau, tau.raw))

		if (output.marginal) {
			results$univariate = analysis$get.univariate(c(results$outcome, results$predictors))
			results$predictor.rg = analysis$get.rg(results$predictors, results$predictors, verbose=F)
		}

		if (output.submodels) {
			results$submodels = catchErrors(analysis$compute.submodels(include.marginal=F), as.warning="no.submodels", all.failures=catch.failures)
			if (!is.null(results$submodels)) {
				results$submodels$coefficients = subset(results$submodels$coefficients, select=-gamma.raw)
				results$submodels$parameters = subset(results$submodels$parameters, select=-c(tau, tau.raw))
			}
		}
		return(results)
	}
}



#' Local partial genetic correlation analysis
#'
#' Will perform a local partial genetic correlation between the first two phenotypes (phen1, phen2) conditioned on the rest (Z).
#' Phenotype order is based on that within the locus object by default, but can be changed by passing a phenotype vector with the desired order to the 'phenos' argument.
#'
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param target The two target phenotypes of interest for which the partial correlation will be computed. All other phenotypes will be conditioned on.
#' @param phenos Subset of phenotypes to analyse. If NULL, all phenotypes in the locus object will be analysed.
#' @param adapt.thresh The thresholds at which to increase the number of iterations for the p-value generation.
#' Default number of iterations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of iterations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param max.r2 Max r2 threshold for the regression of phen1~Z and phen2~Z. If any of these r2's are too high, the partial correlation becomes unstable, and analysis is therefore aborted.
#' @param param.lim The +- threshold at which estimated parameters are considered to be too far out of bounds. If the estimated parameter exceeds this threshold, it is considered unreliable and will be set to NA.
#'
#' @return Data frame with the columns:
#' \itemize{
#'     \item phen1 / phen2 - target phenotypes
#'     \item z - phenotype(s) which the genetic correlation between the target phenotypes were conditioned on
#'     \item r2.phen1_z / r2.phen2_z - the proportion of genetic signal in target phenotypes explained by z. Note: if either of these exceed the threshold specified by the max.r2 argument, the analysis will be aborted (to prevent unstable estimates)
#'     \item pcor - the partial correlation between phen1 and phen2 conditioned on z
#'     \item ci.lower / ci.upper - 95\% confidence intervals for the partial genetic correlation
#'     \item p - simulation p-values for the partial genetic correlation
#' }
#' @export

# NB: can have more than two outcomes now
run.pcor = function(locus.data, outcomes, condition.on=NULL, add.pval=T, add.ci=T, output.marginal=F, output.submodels=F, univ.thresh=1, filter.cond=(univ.thresh < 1), adap.thresh=globals$adaptive.pval.threshold, max.r2=globals$partial.cor.max.r2, param.lim=globals$rg.invalid.bounds, catch.failures=T) {
	if (output.submodels && length(locus.data$match.phenotypes(outcomes)) > 2) input.error("cannot compute submodels if more than two outcomes have been specified")
	analysis = PartialCorrelationAnalysis$new(locus.data, univ.threshold=univ.thresh, adaptive.threshold=ifelse(is.null(adap.thresh), 0, adap.thresh), add.pval=add.pval, add.ci=add.ci, rg.limit=param.lim, maximum.r2=max.r2, filter.variables=filter.cond)
	results = catchErrors(analysis$compute(outcomes, condition.on), all.failures=catch.failures)

	if (!is.null(results)) {
		if (output.marginal) {
			results$univariate = analysis$get.univariate(c(results$outcome, results$predictors))
			if (length(results$conditioned) > 1) results$conditioned.rg = analysis$get.rg(results$conditioned, results$conditioned, rg.limit=param.lim, verbose=F)
		}

		if (output.submodels) results$submodels = catchErrors(analysis$compute.submodels(), as.warning="no.submodels", all.failures=catch.failures)
		return(results)
	}
}







