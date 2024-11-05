## input.info can be either:
## - a (vector of) filename(s)
##   - input.dir can be a vector of same length, to specify different directories for each file (use "" or "." for current directory)
## - a data.frame with the requisite columns
## - a list of PhenotypeInfo objects (output from read.info() / create.info())

#' @export
#' @useDynLib LAVA


#test.add = function(a, b) {.Call("add", a, b, PACKAGE="LAVA")}

process.input = function(input.info, sample.overlap.file, ref.prefix, phenotypes=NULL, chromosomes="all", input.dir=NULL, defer.loading=F) {
	log.message("Processing input", .indent=0)
	if (is.data.frame(input.info)) {
		if (length(input.dir) > 1) input.error("input.dir argument contains multiple elements, with data.frame value for input.info")
		log.message("loading phenotype info from input data.frame")
		pheno.info = validate.info(input.info, phenotypes=phenotypes, input.dir=input.dir)
	} else if (is.list(input.info)) { ## treat as output from create.info
		if (!all(sapply(input.info, has.type, "RawPhenotype"))) input.error("not all elements in input.info list are RawPhenotype information objects")
		pheno.info = input.info
	} else { ## treat as (vector of) filenames
		if (!is.null(input.dir)) {
			input.dir[input.dir == ""] = "."
			if (length(input.dir) == 1) input.dir = rep(input.dir, length(input.info))
			else if (length(input.dir) != length(input.info)) input.error("number of elements of input.dir argument is inconsistent with number of file names in input.info")
		} else input.dir = rep(".", length(input.info))

		check.files.exist(input.info)
		pheno.info = list()
		for (i in seq_along(input.info)) {
			log.message("reading phenotype info file ", input.info[i])
			pheno.info = c(pheno.info, validate.info(input.info[i], phenotypes=phenotypes, input.dir=input.dir[i], skip.unknown=(length(input.info) > 1)))
		}
	}

	names(pheno.info) = sapply(pheno.info, function(pi) {pi$get.name()})
	if (any(duplicated(names(pheno.info)))) input.error("duplicate phenotypes in input")
	if (!is.null(phenotypes)) pheno.info = pheno.info[check.phenotypes(phenotypes, names(pheno.info))]


	sampling.covariance = if (!is.null(sample.overlap.file)) read.sampling.covariance(sample.overlap.file, phenotypes=names(pheno.info))
	ld.reference = read.plink.data(ref.prefix)

	data = DataSet$new(ld.reference)
	for (i in 1:length(pheno.info)) data$add.phenotype(pheno.info[[i]])
	if (!is.null(sampling.covariance)) data$set.sampling.correlation(sampling.covariance)
	data$set.chromosomes(chromosomes)

	if (!defer.loading) {
		data$load()
		if (length(pheno.info) > 1) data$check.overlap()
	}

	return(data)
}


## phenotype can be a named list of phenotype=filename pairs, in which case the filename argument can be omitted
create.info = function(phenotype, filename, ..., input.dir=NULL) {
	if (is.list(phenotype) && !is.null(names(phenotype))) {
		filename = unlist(phenotype)
		phenotype = names(phenotype)
	}

	input.info = list(phenotype=phenotype, filename=filename, ...)
	lengths = sapply(input.info, length)
	if (!all(lengths == 1 | lengths == max(lengths))) input.error("lengths of input are inconsistent")
	return(validate.info(data.frame(input.info), input.dir=input.dir))
}


read.sampling.covariance = function(covar.file, phenotypes=NULL, trim=F) {
	log.message("reading sampling covariance file ", covar.file)
	check.files.exist(covar.file)
	covar = check.covariance(as.matrix(read.table(covar.file, header=T, check.names=F)), allow.NA=T, label="sampling covariance matrix")

	if (!is.null(phenotypes)) {
		exist = phenotypes %in% colnames(covar)
		if (!any(exist)) input.error("none of the input phenotypes are listed in sampling covariance file; please verify that the correct file has been provided")
		if (!all(exist)) throw.warning("%phenotype% ", items.and=phenotypes[!exist], " %is% not listed in sampling covariance file, assuming independent samples")
		if (trim) {
			select = colnames(covar) %in% phenotypes
			covar = covar[select,select]
		}
	}

	return(covar)
}


read.plink.data = function(ref.prefix) {
	log.message("reading PLINK reference data set ", ref.prefix)
	ld.ref = PlinkData$new(ref.prefix)
	log.message("individuals: ", ld.ref$no.indiv(), .indent=2)
	log.message("variants: ", ld.ref$no.snps(), .indent=2)
	return(ld.ref)
}



read.loci = function(loc.file, unit.type=NULL) {
	check.files.exist(loc.file)
	input = TextFileInput$new(loc.file, globals$annot.header.index, input.type="locus")

	if (input$has.parameters("snps")) return(ListAnnotation$new(input$get.subset("locus", "snps"), source=loc.file, type=unit.type))
	else return(RangeAnnotation$new(input$get.subset("locus", "chromosome", "start", "stop"), source=loc.file, type=unit.type))
}








