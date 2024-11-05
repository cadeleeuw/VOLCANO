## validates input data.frame or file and creates corresponding PhenotypeInfo objects
validate.info = function(input, phenotypes=NULL, input.dir=NULL, skip.unknown=F) {
	if (!is.null(input.dir)) check.dirs.exist(input.dir)

	## load in and subset data.frame with info
	if (is.data.frame(input)) input.data = DataframeInput$new(input, globals$info.header.index, input.type="info")
	else input.data = TextFileInput$new(input, globals$info.header.index, input.type="phenotype info")
	input.info = input.data$get.subset("phenotype", "filename", optional=names(globals$info.header.index))  ## select parameters and return as data.frame
	input.info[input.info == "." | input.info == ""] = NA   ## treat . and empty as missing value code,


	## process global parameters, add them as constant columns to
	if ("get.globals" %in% names(input.data)) {
		global.params = input.data$get.globals()
		if (length(global.params) > 0) {
			map = InputInterface$new(names(global.params), globals$info.header.index)$get.map(use="last")
			existing = names(map) %in% names(input.info)
			if (any(existing)) throw.warning("both global %value% and input %column% %has% been specified for %parameter% ", items.and=names(map)[existing], ", global %value% will be ignored")
			params = lapply(map[!existing], function(col) {global.params[[col]]})
			for (param in c("unit", "parameters", "properties")) {if (param %in% names(params)) params[[param]] = paste(params[[param]], collapse=";")}

			multiple = sapply(params, length) > 1
			if (any(multiple)) input.error("global %value% specified for %parameter% ", items.and=names(params)[multiple], " %[each]% %has% multiple elements")
			for (par in names(params)) input.info[[par]] = params[[par]]
		}
	}


	## validate phenotype names and input files
	if (any(duplicated(input.info$phenotype))) input.error("duplicate phenotypes in input")
	if (any(grepl("*", input.info$phenotype, fixed=T))) input.error("character '*' is not allowed in phenotype names")
	if (!is.null(phenotypes)) {
		phenotypes = check.phenotypes(phenotypes, input.info$phenotype, discard.unknown=skip.unknown)
		input.info = input.info[match(phenotypes, input.info$phenotype),]
	}
	if (!is.null(input.dir) && input.dir != ".") input.info$filename = paste0(input.dir, "/", input.info$filename)
	check.files.exist(input.info$filename, resolve.chr=T)


	## determine continuous vs binary, and validate associated parameters
	input.types = c("continuous", "binary")
	if ("input.type" %in% names(input.info)) {
		match = pmatch(input.info$input.type, input.types, duplicates.ok=T)
		invalid = unique(input.info$input.type[!is.na(input.info$input.type) & is.na(match)])
		if (length(invalid) > 0) input.error("unknown input.type %value% ", items.and=invalid)
		input.info$input.type = input.types[match]
	}

	binary.param = c("no.cases", "no.controls", "case.proportion", "prevalence")
	binary.param = binary.param[binary.param %in% names(input.info)]
	for (param in c("sample.size", binary.param)) {if (param %in% names(input.info)) input.info[[param]] = suppressWarnings(as.numeric(input.info[[param]]))}

	if (length(binary.param) > 0) {
		for (param in binary.param) input.info[[param]][!is.na(input.info[[param]]) & input.info[[param]] <= 0] = NA
		is.binary = apply(!is.na(input.info[binary.param]), 1, any)

		if ("input.type" %in% names(input.info)) input.info$input.type[is.na(input.info$input.type)] = input.types[1+is.binary[is.na(input.info$input.type)]]
		else input.info$input.type = input.types[1+is.binary]
	}

	if (!("input.type" %in% names(input.info))) input.error("unable to determine binary/continuous status of any phenotype, please specify")
	unknown = is.na(input.info$input.type)
	if (any(unknown)) input.error("unable to determine binary/continuous status of %phenotype% ", items.and=input.info$phenotype[unknown])


	## create PhenotypeInfo objects
	pheno.info = list()
	for (i in 1:nrow(input.info)) pheno.info[[input.info$phenotype[i]]] = RawPhenotype$new(input.info[i,])

	return(pheno.info)
}


#########################

## helper class to bundle and validate phenotype info; external input needs to contain at least 'phenotype' and 'filename' arguments
PhenotypeInfo = R6::R6Class("PhenotypeInfo",
	private = list(
		name = NULL, raw.info = list(),
		traits = list(continuous=F, binary=F, composite=F, localized=F, truncated=F),

		error = function(...) {input.error(..., .label=list(phenotype=self$get.name()))},

		get.value = function(parameter) {
			if (!(parameter %in% names(private$raw.info))) private$error("missing parameter '", parameter, "' for ", class(self)[1], " object")
			return(private$raw.info[[parameter]])
		},

		set.trait = function(name, bool.value=T) {
			if (!(name %in% names(private$traits))) fatal.error("trying to set unknown phenotype trait '", name, "'")
			private$traits[[name]] = bool.value
			if (name == "continuous") private$traits$binary = !bool.value
			if (name == "binary") private$traits$continuous = !bool.value
		},

		get.printer = function() {undefined.error("get.printer", class(self)[1])}
	),
	public = list(
		initialize = function(...) {
			private$raw.info = flatten.arglist(..., filter.NA=T)
			private$name = private$get.value("phenotype")

			unknown = !(names(private$raw.info) %in% names(globals$info.header.index))
			if (any(unknown)) private$error("unknown %parameter% ", items.and=names(private$raw.info)[unknown], " for ", class(self)[1], " object")
		},

		abbreviate = function() {return(paste0(self$get.name(), ", ", self$get.traits()))},
		print = function(...) {cat(private$get.printer()$to.string())},

		get.name = function() {return(private$name)},
		get.traits = function(as.string=T, sep=" ") {
			traits = names(private$traits)[unlist(private$traits)]
			if (as.string) return(ifelse(length(traits) > 0, paste(traits, collapse=sep), "generic"))
			else return(traits)
		},

		is = function(trait) {
			if (!(trait %in% names(private$traits))) fatal.error("trying to query unknown phenotype trait '", trait, "'")
			return(private$traits[[trait]])
		}
	)
)


## info for phenotype directly read from data
## validation of input parameters is performed *here* (SummaryStatistics objects assume them to be valid without further checks), with get.configuration() exposing a list with all available parameters
## underlying data is exposed via get.data()

## possible configuration components
## - core: snp.id, allele1, allele2
## - metrics: sample.size, no.cases, no.controls, case.proportion (max. 2)
## - statistics: statistic -or- p.value, beta/log.odds/odds.ratio
## - units: [localized phenotype units]

## valid binary column configurations (in order checked)
## - two columns of sample.size, no.cases, no.controls
## - sample.size and case.proportion columns
## - sample size column and either global case.proportion or two of global sample.size, no.cases, no.controls (will store global case.proportion)
## - two of global sample.size, no.cases, no.controls (will store sample.size and case.proportion)
## - global sample.size and case.proportion
RawPhenotype = R6::R6Class("RawPhenotype",
	inherit = PhenotypeInfo,
	private = list(
		filename = NULL,
		unit.info = list(names=NULL, columns=list(), values=list()),  ## units for localized phenotypes
		sample.metrics = list(sample.size=NULL, case.proportion=NULL, prevalence=NULL),  ## global values, if provided and not available through input column in data

		configuration = list(core=NULL, metrics=NULL, statistics=NULL, units=NULL), ## list with available parameters in input data, by category (see comment above for possible parameters)
		data.interface = list(main=NULL, chromosomes=list()),  ## TextFileInputInterface object for input file, separate objects per chromosome if split

		## define additional parameters for InputInterface objects
		add.parameters = function(update.params) {
			if (length(update.params) > 0) {
				private$data.interface$main$update.parameters(update.params, append.existing=T)
				for (i in seq_along(private$data.interface$chromosomes)) {
					if (!is.null(private$data.interface$chromosomes[[i]])) private$data.interface$chromosomes[[i]]$update.parameters(update.params, append.existing=T)
				}
			}
		},

		## validate input file/files (expanding [CHR] placeholder) and create TextFileInputInterface object(s)
		## TextFileInputInterface objects read in file header for validation of required columns, and provides interface to full data file object
		validate.files = function() {
			private$filename = private$get.value("filename")

			check.files.exist(private$filename, resolve.chr=T)
			if (self$by.chromosome()) {
				chr.files = chromosome.files(private$filename, prune.missing=T); columns = list()
				for (i in 1:nrow(chr.files)) {
					curr = TextFileInputInterface$new(chr.files$file[i], globals$sumstats.header.index)
					private$data.interface$chromosomes[[chr.files$chromosome[i]]] = curr
					columns[[i]] = sort(unique(tolower(curr$get.variables())))
				}
				if (length(columns) > 1 && !all(sapply(columns[-1], function(col) {length(col) == length(columns[[1]]) && all(col == columns[[1]])}))) private$error("headers of chromosome files for ", private$filename, " are not consistent")
				private$data.interface$main = private$data.interface$chromosomes[[chr.files$chromosome[1]]]  ## set to first chromosome for subsequent column validation
			} else private$data.interface$main = TextFileInputInterface$new(private$filename, globals$sumstats.header.index)
		},

		## load in and check general properties ('properties' parameter)
		validate.properties = function(recognized=c()) {
			properties = tryCatch(parse.value.string(private$raw.info$properties), error=function(err) {private$error(err)})
			unknown = !(names(properties) %in% recognized)
			if (any(unknown)) private$error("%property% ", items.and=names(properties)[unknown], " specified in 'properties' %is% unknown")
			if ("truncated" %in% names(properties)) private$set.trait("truncated")
		},

		## validate user-defined input columns
		validate.columns = function() {
			update.params = tryCatch(parse.value.string(private$raw.info$parameters, allow.unnamed=F, infer.names=F), error=function(err) {private$error(err)})
			unknown = !private$data.interface$main$has.parameters(names(update.params), aggregate=F)
			if (any(unknown)) private$error("%parameter% ", items.and=names(update.params)[unknown], " specified in 'parameters' %is% unknown")
			private$add.parameters(update.params)
		},

		## validate localized phenotype unit settings, and define corresponding parameters in InputInterface
		validate.localized = function() {
			private$set.trait("localized")
			units = tryCatch(parse.value.string(private$raw.info$unit, allow.unnamed=F, infer.names=T), error=function(err) {private$error(err)})
			units = lapply(units, function(u) {unique(tolower(u))})

			duplicate.names = duplicated(names(units))
			if (any(duplicate.names)) private$error("unit %name% ", items.and=unique(names(units)[duplicate.names]), " %is% used more than once")
			duplicate.columns = duplicated(unlist(units))
			if (any(duplicate.columns)) private$error("%column% ", items.and=unique(unlist(units)[duplicate.columns]), " %is% %[each]% specified for multiple units")

			exist = private$data.interface$main$update.parameters(units, append.existing=F)
			if (length(exist) > 0) private$error("parameter %name% ", items.and=exist, " %is% already in use")
			available = private$data.interface$main$has.available(names(units))
			if (!all(available)) private$error("%column% specified for %unit% ", items.and=names(units)[!available], " %does% not exist")

			private$add.parameters(units)
			private$unit.info = list(names=names(units), columns=units, values=list())
			private$configuration$units = names(units)
		},

		determine.columns = function() {
			## validate core parameters
			private$configuration$core = c("snp.id", "allele1", "allele2")
			available = private$data.interface$main$has.available(private$configuration$core)
			if (!all(available)) private$error("%column% for required %parameter% ", items.and=private$configuration$core[!available], " %is% missing")


			## validate statistics
			if (!private$data.interface$main$has.available("statistic")) {
				if (private$data.interface$main$has.available("p.value")) {
					effect.params = c(if (self$is("binary")) c("log.odds", "odds.ratio"), "beta"); available = private$data.interface$main$has.available(effect.params)
					if (any(available)) private$configuration$statistics = c("p.value", effect.params[available][1])
					else private$error("column for %[+one of]% %parameter% ", items.or=effect.params, " is required if no column for parameter 'statistic' is present")
				} else private$error("column for parameter 'p.value' is required if no column for parameter 'statistic' is present")
			} else private$configuration$statistics = "statistic"


			## validate metrics
			if (self$is("binary")) {
				binary.params = c("sample.size", "no.cases", "no.controls", "case.proportion")
				has.column = private$data.interface$main$has.available(binary.params)

				if (sum(has.column[-4]) >= 2) private$configuration$metrics = head(binary.params[has.column], 2)   ## use first two columns of sample.size, no.cases, no.controls
				else if (has.column[1] && has.column[4]) private$configuration$metrics = c("sample.size", "case.proportion")   ## use columns for sample.size and case proportion
				else {  ## don't have enough information in input file columns, add global values
					has.global = binary.params %in% names(private$raw.info)

					## get sample size if directly provided, or if have no.case and no.control
					sample.size = ifelse(has.global[1], private$raw.info$sample.size, ifelse(all(has.global[2:3]), private$raw.info$no.cases + private$raw.info$no.controls, NA))

					## compute
					case.proportion = ifelse(all(has.global[1:3] == c(F,T,T)), private$raw.info$no.cases/sample.size, ifelse(has.global[4], private$raw.info$case.proportion, NA))

					if (is.na(case.proportion) && !is.na(sample.size)) {
						if (has.global[2]) case.proportion = private$raw.info$no.cases / sample.size
						else if (has.global[3]) case.proportion = 1 - private$raw.info$no.cases / sample.size
					}

					if (is.na(case.proportion) || (!has.column[1] && is.na(sample.size))) private$error("please provide valid combination of input columns and/or global values for %parameter% ", items.and=binary.params)
					if (has.column[1]) private$configuration$metrics = "sample.size"
					else private$sample.metrics$sample.size = sample.size

					private$sample.metrics$case.proportion = case.proportion
				}

				if ("prevalence" %in% names(private$raw.info)) private$sample.metrics$prevalence = private$raw.info$prevalence
			} else {
				binary.params = c("no.cases", "no.controls", "case.proportion", "prevalence")
				if (any(binary.params %in% names(private$raw.info))) private$error("%parameter% ", items.and=binary.params[binary.params %in% names(private$raw.info)], " %is% not valid for continuous phenotype")

				if (is.null(private$raw.info$sample.size)) {
					if (!private$data.interface$main$has.available("sample.size")) private$error("missing both column and global value for required parameter 'sample.size'")
					private$configuration$metrics = "sample.size"
				} else private$sample.metrics$sample.size = private$raw.info$sample.size
			}
		},

		get.printer = function() {
			printer = ObjectPrinter$new("PhenotypeInfo (raw input)")$parameter.list(name=self$get.name(), file=private$filename, type=self$get.traits())
			if (self$is("localized")) {
				units = ifelse(length(private$unit.info$values) > 0, paste0(private$unit.info$names, "=", private$unit.info$values, collapse=", "), paste(private$unit.info$names, collapse=", "))
				printer$add.parameter("localized units", units)
			}
			printer$add.parameter("sample size", private$sample.metrics$sample.size)

			if (self$is("binary")) {
				printer$add.line("case proportion: ", round(100*private$sample.metrics$case.proportion, 2), "%", add.space=F)
				printer$parameter.list(prevalence=private$sample.metrics$prevalence)
			}
			return(printer)
		}
	),
	public = list(
		initialize = function(...) {
			super$initialize(...)

			input.type = private$raw.info$input.type
			if (length(input.type) != 1 && input.type %in% c("continuous", "binary")) private$error("no valid input.type specified")
			private$set.trait(input.type)

			private$validate.files()
			if ("properties" %in% names(private$raw.info)) private$validate.properties(recognized = c("truncated"))
			if ("parameters" %in% names(private$raw.info)) private$validate.columns()
			if ("unit" %in% names(private$raw.info)) private$validate.localized()

			private$determine.columns()
		},

		## get.name creates full name with localized units included, get.basename only returns core phenotype name
		get.name = function() {return(ifelse(self$is("localized") && length(private$unit.info$values) > 0, paste0(c(private$name, private$unit.info$values), collapse="::"), private$name))},
		get.basename = function() {return(private$name)},
		abbreviate = function() {return(paste0(super$abbreviate(), " (", private$filename, ")"))},

		get.filename = function(chromosomes=NULL) {return(if(!is.null(chromosomes)) gsub("[CHR]", chromosome, private$filename, fixed=T) else private$filename)},

		get.global = function(name) {out = private$sample.metrics[[name]]; return(if (!is.null(out)) out else NA)},
		get.unit.info = function() {return(private$unit.info)},

		## return list of available parameter names per category; convert to data.frame with corresponding column name if add.columns=T
		get.configuration = function(add.columns=F) {
			if (add.columns) {
				out = data.frame(category=unlist(lapply(names(private$configuration), function(c) {rep(c, length(private$configuration[[c]]))})), stringsAsFactors=F)
				out$parameter = unlist(private$configuration)
				out$column = unlist(private$data.interface$main$get.map(out$parameter, use="first", drop.missing=F))
				return(out)
			} else return(private$configuration)
		},

		get.data = function(chromosomes=NULL) {
			if (!is.null(chromosomes)) {
				if (!self$by.chromosome()) private$error("cannot load by chromosome for whole genome input")
				chromosomes = validate.chromosomes(chromosomes, condense=F)
				available = self$has.chromosomes(chromosomes, aggregate=F)
				if (!all(available)) private$error("cannot load data for %chromosome% ", items.and=chromosomes[!available])
				interfaces = private$data.interface$chromosomes[chromosomes]
			} else if (is.null(chromosomes)) {
				if (self$by.chromosome()) interfaces = private$data.interface$chromosomes[!unlist(lapply(private$data.interface$chromosomes, is.null))]
				interfaces = list(private$data.interface$main)
			}

			if (length(interfaces) == 0) private$error("no data to load")
			for (i in seq_along(interfaces)) {
				curr.input = interfaces[[i]]$load.input(input.type = "summary statistics", error.tag = paste0("phenotype: ", self$get.name(), "; file: ", interfaces[[i]]$get.filename()))
				curr.data = curr.input$get.subset(unlist(private$configuration))
				if (i > 1) {
					if (!identical(names(input.data), names(curr.data))) private$error("cannot merge chromosome input data, headers are incompatible", filename=private$filename)
					else input.data = rbind(input.data, curr.data)
				} else input.data = curr.data
			}

			return(input.data)
		},

		by.chromosome = function() {return(grepl("[CHR]", private$filename, fixed=T))},
		has.chromosomes = function(chromosomes, aggregate=F) {
			if ("all" %in% chromosomes) aggregate = T
			all.avail = which(sapply(1:23, function(i) {length(private$data.interface$chromosomes) >= i && !is.null(private$data.interface$chromosomes[[i]])}))
			available = validate.chromosomes(chromosomes, condense=F) %in% all.avail
			if (aggregate) return(all(available))
			else return(available)
		},

		## return PhenotypeInfo object with localized units set to specific value
		instantiate.units = function(...) {
			values = flatten.arglist(...)
			if (length(values) != length(private$unit.info$names) || any(sapply(values, length) != 1) || !(all(private$unit.info$names %in% names(values)))) private$error("cannot instantiate Phenotype object, unit specification is invalid")
			out = self$clone()
			out$.__enclos_env__$private$unit.info$values = values[private$unit.info$names]
			return(out)
		}
	)
)



## container class for ad hoc combination of other phenotypes
# CompositePhenotype = R6::R6Class("CompositePhenotype",
# 	inherit = PhenotypeInfo,
# 	private = list(
# 		components = NULL,
# 		get.printer = function() {return(ObjectPrinter$new("PhenotypeInfo (composite)")$parameter.list(name=self$get.name(), type=self$get.traits()))}
# 	),
# 	public = list(
# 		initialize = function(name, components) {
# 			super$initialize(phenotype=name)
# 			private$set.trait("composite")
# 			private$components = components
# 		},
#
# 		get.components = function() {return(components)}
# 	)
# )








