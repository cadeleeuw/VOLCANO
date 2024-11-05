AnalysisProcessor = R6::R6Class("AnalysisProcessor",
	private = list(
		locus.data = NULL, univariate = NULL, results = NULL,
		settings = list(), settings.restore = list(),

		init.settings = function(..., default.settings=NULL, import.settings=NULL) {
			private$settings = list(
				univ.threshold = 1,
				add.pval = T,
				add.ci = F,
				adaptive.threshold = globals$adaptive.pval.threshold,
				signif.decimals = globals$analysis.signif.decimals
			)
			if (!is.null(default.settings))	for (par in names(default.settings)) private$settings[[par]] = default.settings[[par]]
			private$settings.restore$default = private$settings

			if (!is.null(import.settings)) for (par in names(import.settings)) private$settings[[par]] = import.settings[[par]]
			private$process.settings(...)
			private$settings.restore$initial = private$settings
		},

		process.settings = function(...) {
			args = flatten.arglist(...)

			valid.params = c(names(private$settings), private$valid.params)
			invalid = !(names(args) %in% valid.params)
			if (any(invalid)) private$error("unknown %parameter% ", items.and=names(args)[invalid])
			for (par in names(args)) private$settings[[par]] = args[[par]]
		},

		## if lines is a named list, will print "[name] ([value])"
		format.condition = function(type.func, ..., phenotype=NULL, lines=NULL) {
			if (!is.null(lines) && is.list(lines)) lines = paste0(names(lines), " (", unlist(lines), ")")
			label = list(phenotype=phenotype, "locus ID"=private$locus.data$get.data()$get.locus()$get.name())
			type.func(..., .label=label, .lines=lines)
		},

		warning = function(..., phenotype=NULL, lines.list=NULL) {private$format.condition(throw.warning, ..., phenotype=phenotype, lines=lines.list)},
		error = function(..., phenotype=NULL, lines.list=NULL) {private$format.condition(input.error, ..., phenotype=phenotype, lines=lines.list)},
		failure = function(..., phenotype=NULL, lines.list=NULL) {private$format.condition(data.error, ..., phenotype=phenotype, lines=lines.list)},

		check.status = function(phenotypes, check.univ=T) {
			if (!is.null(phenotypes)) {
				if (length(phenotypes) == 0) return(list(available = character(0), failed = list()))
				info = private$univariate[match(phenotypes, private$univariate$phenotype),]
			} else info = private$univariate

			failed = list()
			for (i in 1:nrow(info)) {
				if (info$status[i] == "failed") failed[[info$phenotype[i]]] = "could not be processed in locus"
				else if (info$status[i] == "negative") failed[[info$phenotype[i]]] = "local genetic variance is negative"
				else if (check.univ && !is.null(private$settings$univ.threshold)) {
				 	if (!is.na(info$p.univariate[i]) && info$p.univariate[i] > private$settings$univ.threshold) failed[[info$phenotype[i]]] = "insufficient genetic signal"
				}
			}
			return(list(available = info$phenotype[!(info$phenotype %in% names(failed))], failed = failed))
		},

		## match partial names, check and throw error for ambiguous and unknown, remove duplicates, expand wildcards
		## if fill.null=T, expands NULL phenotypes argument to list of all phenotypes
		## otherwise, returns NULL if input is NULL, or character(0) if input is empty
		check.phenotypes = function(phenotypes, fill.null=F) {
			if (length(phenotypes) > 0) {
				res = tryCatch(check.phenotypes(phenotypes, private$univariate$phenotype, label="analysis input"), error=function(err) {private$error(err)})
				return(res)
			} else {
				if (is.null(phenotypes)) {
					if (fill.null) return(private$univariate$phenotypes)
					else return(NULL)
				} else return(character(0))
			}
		},

		## set named arguments with 'column=bounds', where bounds is considered symmetric around zero if only a single value
		## truncates columns with given name to specified bounds and rounds off to settings$signif.decimals with signif() (set bound to Inf to round only)
		## ignores columns if not present, rounds .lower and .upper columns to same level if present
		format.results = function(results, ...) {
			params = lapply(list(...), function(x) {return(if(length(x) == 1) c(-x,x) else x[1:2])})
			for (param in names(params)) {
				for (col in paste0(param, c("", ".lower", ".upper"))) {
					if (col %in% names(results)) {
						curr = results[[col]]
						curr[!is.na(curr) & curr < params[[param]][1]] = params[[param]][1]
						curr[!is.na(curr) & curr > params[[param]][2]] = params[[param]][2]
						results[[col]] = signif(curr, private$settings$signif.decimals)
					}
				}
			}

			return(results)
		}
	),
	public = list(
		initialize = function(locus.data, ..., default.settings=NULL, import.settings=NULL) {
			check.types(locus.data="ProcessedLocus", .class=class(self))
			private$locus.data = locus.data

			private$init.settings(..., default.settings=default.settings, import.settings=import.settings)
			if (private$locus.data$no.pheno() == 0) private$failure("input locus data contains no phenotypes to analyse")

			private$univariate = private$locus.data$get.marginal()
		},

		no.pheno = function() {return(private$locus.data$no.pheno())},
		phenotypes = function() {return(private$locus.data$phenotypes())},

		get.univariate = function(phenotypes=NULL) {
			if (length(phenotypes) > 0) {
				out = private$univariate[match(private$check.phenotypes(phenotypes), private$univariate$phenotype),]
				rownames(out) = NULL
				return(out)
			} else return(private$univariate)
		},

		## subset to specified rows/columns, use all if NULL
		get.omega = function(pheno.cols=NULL, pheno.rows=NULL, scale=c("component", "delta")) {
			omega = private$locus.data$get.omega(match.arg(scale))
			return(matrix.subset(omega, cols=private$check.phenotypes(pheno.cols), rows=private$check.phenotypes(pheno.rows)))
		},

		get.rg = function(pheno.cols=NULL, pheno.rows=NULL, rg.limit=globals$rg.invalid.bounds, verbose=T) {
			rg = matrix.subset(suppressWarnings(cov2cor(private$locus.data$get.omega())), cols=private$check.phenotypes(pheno.cols), rows=private$check.phenotypes(pheno.rows))
			invalid = !is.na(rg) & abs(rg) > rg.limit & rownames(rg)[row(rg)] != colnames(rg)[col(rg)]
			if (any(invalid)) {
				pairs = cbind(rownames(rg)[row(rg)[invalid]], colnames(rg)[col(rg)[invalid]])
				pairs = matrix(pairs[!duplicated(apply(pairs, 1, function(x) {paste(sort(x), collapse=" ")})),], ncol=2)
				if (nrow(pairs) > 0) {
					if (verbose) private$warning("correlation estimates too far out of bounds (+/-", rg.limit, ") for following phenotype pair(s); values will be set to NA", lines.list=paste(pairs[,1], pairs[,2], sep=" ~ "))
					rg[invalid] = NA
				}
			}
			return(rg)
			return(matrix.subset(rg, cols=private$check.phenotypes(pheno.cols), rows=private$check.phenotypes(pheno.rows)))
		},

		get.results = function() {return(private$results)},

		get.settings = function() {private$settings[substr(names(private$settings), start=1, stop=1) != "."]},

		set = function(..., init=F) {
			if (init) self$restore.settings("initial")
			private$process.settings(...)
			invisible(self)
		},
		restore.settings = function(restore=c("initial", "default")) {private$settings = private$settings.restore[[match.arg(restore)]]; invisible(self)}
	)
)



