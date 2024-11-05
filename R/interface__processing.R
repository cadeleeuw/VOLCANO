#TODO update
process.locus = function(locus.info, data.set, phenotypes=NULL, units=NULL, prune.threshold=0.99, catch.failures=T) {
	input = tryCatch({
			out = data.set$get.interface(phenotypes=phenotypes, prune.threshold=prune.threshold)$set.locus(locus.info, detect.units=T)
			if (!is.null(units)) out$set.units(units)
			out
		}, failure = function(err) {if (catch.failures) throw.warning(err) else throw.error(err)}
	)

	if (!is.null(input)) {
		output = tryCatch(input$process(), failure = function(err) {if (catch.failures) throw.warning(err) else throw.error(err)})
		return(output)
	}
}



#' Re-process locus to meta-analyse of selected phenotypes
#'
#' Will combine all elements of the requested phenotypes using standard inverse variance weighting, allowing them to be analysed as a single phenotype via the  multivariate analysis functions.
#' Note that the univariate test cannot currently be applied to meta-analysed phenotypes, so please do that beforehand on each phenotype individually.
#'
#' @param locus Locus object defined using the \code{\link{process.locus}} function.
#' @param meta.phenos Phenotypes you want to meta-analyse
#'
#' #'
#' @return This function returns an object just like that \code{\link{process.locus}} function, containing general locus info, the relevant processed sumstats, and info about the input phenotypes.
#'
#' \itemize{
#'     \item id - locus ID
#'     \item chr/start/stop - locus coordinates
#'     \item snps - list of locus SNPs
#'     \item N.snps - number of SNPs
#'     \item K - number of PCs
#'     \item delta - PC projected joint SNP effects for each phenotype
#'     \item sigma - sampling covariance matrix
#'     \item omega - genetic covariance matrix
#'     \item omega.cor - genetic correlation matrix
#'     \item N - vector of average N across locus SNPs for each phenotype
#'     \item phenos - phenotype IDs
#'     \item binary - boolean vector indicating whether phentoypes are binary
#' }
#'
#' @export


## phenotypes argument can be either a vector of names, or a named list
## if names is set, this overrides any names in the list
meta.analyse = function(locus.data, phenotypes, meta.names=NULL, drop.components=T) {
	if (!is.list(phenotypes)) phenotypes = list(phenotypes)
	if (!is.null(meta.names)) {
		if (length(meta.names) != length(phenotypes)) input.error("length of phenotypes argument does not match number of names provided")
		names(phenotypes) = meta.names
	}

	if (is.null(names(phenotypes))) names(phenotypes) = rep("", length(phenotypes))
	unnamed = names(phenotypes) == ""
	if (any(unnamed)) {
		used.index = as.numeric(gsub("meta", "", grep("^meta[0-9]*$", c(locus.data$phenotypes(), names(phenotypes)), value=T)))
		available = which(!(1:(sum(unnamed) + length(used.index)) %in% used.index))
		names(phenotypes)[unnamed] = paste0("meta", available[1:sum(unnamed)])
	}

	phenotypes = lapply(phenotypes, function(p) {locus.data$match.phenotypes(p)})
	for (meta in names(phenotypes)) locus.data = locus.data$meta.analyse(meta, phenotypes[[meta]], drop.components=F)
	if (drop.components) locus.data = locus.data$subset(unlist(phenotypes), include=F)

	return(locus.data)
}





