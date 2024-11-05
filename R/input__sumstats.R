## reads in summary statistics and aligns them to the LD reference, filtering out variants not present in reference and entries with NA values
## provided counts are:
##    rows: number of variants read from input file
##    valid: number remaining after filtering on missing summary statistics, other missing values in relevant columns, zero sample size, and duplicate SNP IDs
##    available: number remaining after filtering those not in reference
##    aligned: number remaining after filtering those that could not be aligned to reference (allele mismatch, ambiguous AT/CG)
SummaryStatistics = R6::R6Class("SummaryStatistics",
	private = list(
		phenotype = NULL, filename = NULL, pheno.info = NULL,
		ld.reference = NULL, sum.stats = NULL,
		chromosomes = NULL, counts = list(rows=NA, valid=NA, available=NA, aligned=NA),

		warning = function(..., filename=NULL) {throw.warning(..., .label=list(phenotype=private$phenotype, file=ifelse(is.null(filename), private$filename, filename)))},
		error = function(..., filename=NULL) {input.error(..., .label=list(phenotype=private$phenotype, file=ifelse(is.null(filename), private$filename, filename)))},

		get.printer = function() {
			printer = ObjectPrinter$new(class(self)[1])$parameter.list(phenotype=private$phenotype, file=private$filename, type=self$get.type())
			printer$add.line("variants:", private$counts$rows, "in file,", private$counts$aligned, "valid and aligned to reference")
			return(printer)
		},

		load.data = function(chromosomes) {
			if (private$pheno.info$by.chromosome()) {
				chromosomes = chromosome.files(private$filename, chromosomes, prune.missing=T)$chromosome
				if (length(chromosomes) > 0) {
					for (i in seq_along(chromosomes)) {
						curr = private$load.file(chromosomes[i])
						if (i > 1) {
							if (identical(names(private$sum.stats), names(curr))) private$sum.stats = rbind(private$sum.stats, curr)
							else private$error("cannot merge chromosomes, headers are incompatible", filename=private$pheno.info$get.filename(chromosomes[i]))
						} else private$sum.stats = curr
					}
				} else private$error("specified input files do not exist")
			} else private$sum.stats = private$load.file()

			private$counts$rows = nrow(private$sum.stats)
			private$post.process()
		},

		## to keep consistent counts load.file function shouldn't filter rows, put this in post.process()
		## NB: configuration assumed to have been validated in pheno.info object
		load.file = function(chromosome=NULL) {
			config = private$pheno.info$get.configuration()
			sum.stats = private$pheno.info$get.data(chromosome)
			if (is.null(sum.stats)) private$error("unable to load input data")

			sum.stats$snp.id = tolower(sum.stats$snp.id)

			if (!("statistic" %in% names(sum.stats))) {
				pvals = suppressWarnings(as.numeric(sum.stats$p.value))
				pvals[!is.na(pvals) & pvals < globals$min.pvalue] = globals$min.pvalue

				effect.param = config$statistics[2]
				if (effect.param == "odds.ratio") {
					direction = sign(sum.stats[[effect.param]] - 1)
					if (any(!is.na(direction) & direction < 0)) private$warning("input contains negative odds ratio values; please check that you did not provide log odds or betas", filename=private$pheno.info$get.filename(chromosome))
				}	else direction = sign(sum.stats[[effect.param]])

				avg.direction = mean(direction == 1)
				if (avg.direction < globals$stat.direction.avg.threshold || avg.direction > (1-globals$stat.direction.avg.threshold)) private$warning("almost all ", effect.param, " values have the same direction; please verify input is correct", filename=private$pheno.info$get.filename(chromosome))
				sum.stats$statistic = -qnorm(pvals/2) * direction

				sum.stats[config$statistics] = NULL
			} else sum.stats$statistic = suppressWarnings(as.numeric(sum.stats$statistic))

			if (private$pheno.info$is("binary")) {
				binary.columns = config$metrics
				if (length(binary.columns) >= 2) {
					if (!all(c("sample.size", "case.proportion") %in% binary.columns)) {
						if (all(c("sample.size", "no.cases") %in% binary.columns)) sum.stats$case.proportion = suppressWarnings(as.numeric(sum.stats$no.cases)) / sum.stats$sample.size
						else if (all(c("no.cases", "no.controls") %in% binary.columns)) sum.stats$case.proportion = 1 - suppressWarnings(as.numeric(sum.stats$no.controls)) / sum.stats$sample.size
						else {
							no.cases = suppressWarnings(as.numeric(sum.stats$no.cases)); no.controls = suppressWarnings(as.numeric(sum.stats$no.controls))
							sum.stats$sample.size = no.cases + no.controls
							sum.stats$case.proportion = no.cases / sum.stats$sample.size
						}
					} else sum.stats$case.proportion = suppressWarnings(as.numeric(sum.stats$case.proportion))

					if (!(all(c("sample.size", "case.proportion") %in% names(sum.stats)))) private$error("invalid parameter configuration for binary phenotype")
				}
			}

			if (nrow(sum.stats) == 0) private$error("file is empty", filename=private$pheno.info$get.filename(chromosome))
			return(sum.stats)
		},

		## dupl.columns sets additional columns to use to identify rows, in addition to snp.id
		post.process = function(dupl.columns=NULL) {
			config = private$pheno.info$get.configuration()
			private$counts$rows = nrow(private$sum.stats)

			## filter duplicate SNP IDs
			dupl.id = paste.columns(private$sum.stats, c(config$core[1], dupl.columns), sep=" ")
			duplicate = dupl.id %in% unique(dupl.id[duplicated(dupl.id)])
			if (any(duplicate)) {
				private$warning("found ", sum(duplicate), " lines with duplicate SNP IDs, discarding these from input")
				private$sum.stats = private$sum.stats[!duplicate,]
			}

			## convert infinite valued test statistics and filter out invalid entries (missing values in any used column, invalid sample size)
			is.infinite = which(abs(private$sum.stats$statistic) > abs(qnorm(globals$min.pvalue/2)))
			private$sum.stats$statistic[is.infinite] = abs(qnorm(globals$min.pvalue/2)) * sign(private$sum.stats$statistic[is.infinite])
			keep = column.none(is.na(private$sum.stats))
			if ("sample.size" %in% names(private$sum.stats)) keep = keep & private$sum.stats$sample.size > 0
			private$sum.stats = private$sum.stats[keep,]
			private$counts$valid = nrow(private$sum.stats)

			## map to reference data, filtering out unmatched SNPs (adds direction column)
			private$sum.stats = private$ld.reference$align.sumstats(private$sum.stats)
			private$counts$available = nrow(private$sum.stats)

			## align summary statistics, filtering out unalignable SNPs
			private$sum.stats$statistic = private$sum.stats$statistic * private$sum.stats$direction
			private$sum.stats = private$sum.stats[!is.na(private$sum.stats$statistic),]
			private$sum.stats$direction = NULL
			private$counts$aligned = nrow(private$sum.stats)

			sq.stat.avg = mean(private$sum.stats$statistic^2)
			if (sqrt(sq.stat.avg) < globals$stat.deflation.threshold) private$warning("SNP summary statistics appear severely deflated; please check that input is correct")
		},

		process.statistics = function(stats, select.ids, add.info) {
			if (!is.null(select.ids)) {
				if (has.type(select.ids, "SNPindex")) select.ids = select.ids$snps()
				stats = stats[stats$internal.id %in% select.ids,]
			}
			if (add.info) stats = cbind(private$ld.reference$get.snp.info(stats$internal.id)[,c("snp.id", "chromosome", "position", "allele1", "allele2")], subset(stats, select=-internal.id))
			return(stats)
		}
	),
	public = list(
		initialize = function(pheno.info, ld.reference, chromosomes="all") {
			private$phenotype = pheno.info$get.name(); private$filename = pheno.info$get.filename();
			private$pheno.info = pheno.info; private$ld.reference = ld.reference;

			if (pheno.info$is("localized") && !has.type(self, "LocalizedSummaryStatistics")) fatal.error("localized phenotype summary statistics require LocalizedSummaryStatistics object")

			private$load.data(chromosomes)
		},

		print = function(...) {cat(private$get.printer()$to.string())},

		get.statistics = function(select.ids=NULL, add.info=F) {return(private$process.statistics(private$sum.stats, select.ids, add.info))},
		get.reference = function() {return(private$ld.reference)},

		get.phenotype = function() {return(private$phenotype)},
		get.type = function() {return(private$pheno.info$get.traits())},

		is.binary = function() {return(private$pheno.info$is("binary"))},
		is.localized = function() {return(private$pheno.info$is("localized"))},

		get.counts = function() {return(private$counts)},
		get.info = function() {return(private$pheno.info)}
	)
)


LocalizedSummaryStatistics = R6::R6Class("LocalizedSummaryStatistics",
	inherit = SummaryStatistics,
	private = list(
		units = NULL,

		get.printer = function() {
			printer = super$get.printer()
			for (unit in private$units) printer$add.line("localized units: ", private$counts$units[[unit]], " (", unit, ")", add.space=F)
			return(printer)
		},

		check.units = function(unit.names) {
			matched = unit.names %in% private$units
			if (!all(matched)) private$error("invalid unit specification for localized phenotype, unknown %unit% ", items.and=names(args)[!matched])
		},

		unit.index = function(...) {
			args = flatten.arglist(...)
			private$check.units(names(args))

			index = rep(T, nrow(private$sum.stats))
			for (unit in names(args)) {
				if (length(args[[unit]]) > 1) index = index & private$sum.stats[[unit]] %in% args[[unit]]
				else index = index & private$sum.stats[[unit]] == args[[unit]]
			}
			return(index)
		},

		post.process = function() {
			private$units = private$pheno.info$get.configuration()$units
			for (unit in private$units) {
				is.repeat = private$sum.stats[[unit]] == "."; is.repeat[1] = F
				index = which(!is.repeat); count = c(index[-1], length(is.repeat)+1) - index
				private$sum.stats[[unit]] = private$sum.stats[[unit]][rep(index, times=count)]
			}
			failed = apply(private$sum.stats[private$units] == ".", 1, any)
			if (any(failed)) private$sum.stats = private$sum.stats[!failed,]

			super$post.process(dupl.columns=private$units)

			private$counts$units = list()
			for (unit in private$units) private$counts$units[[unit]] = length(unique(private$sum.stats[[unit]]))
		}
	),
	public = list(
		initialize = function(pheno.info, ld.reference, chromosomes="all") {
			super$initialize(pheno.info, ld.reference, chromosomes=chromosomes)
			if (!pheno.info$is("localized")) fatal.error("input to LocalizedSummaryStatistics object is not localized")
		},

		get.units = function() {return(private$units)},

		get.unit.values = function(units=NULL) {
			if (is.null(units)) units = private$units
			else private$check.units(units)

			values = data.frame(unique(private$sum.stats[units]), stringsAsFactors=F)
			names(values) = units; rownames(values) = NULL
			return(values)
		},

		## ... argument specifies unit.name=values pairs (can be NULL, values can be vector of multiple)
		get.statistics = function(..., select.ids=NULL, add.info=F) {return(private$process.statistics(private$sum.stats[private$unit.index(...),], select.ids, add.info))},

		create.annotation = function(unit.name) {
			if (length(unit.name) > 1) private$error("cannot create multi-unit annotation")
			private$check.units(unit.name)

			locus.names = unique(private$sum.stats[[unit.name]])
			location = private$ld.reference$get.snp.info(private$sum.stats$internal.id)[,c("chromosome", "position")]
			chr = aggregate(location$chromosome, list(unit=private$sum.stats[[unit.name]]), function(x) {unique(x)})
			if (any(sapply(chr$x, length) != 1)) private$error("cannot create annotation for unit ", unit.name, ", some instances span multiple chromosomes")
			ranges = aggregate(location$position, list(unit=private$sum.stats[[unit.name]]), range)

			annot = data.frame(locus.names, chr[match(locus.names, chr$unit),2], ranges[match(locus.names, ranges$unit),2], stringsAsFactors=F)
			names(annot) = c("locus", "chromosome", "start", "stop")
			return(RangeAnnotation$new(annot, source=private$pheno.info$get.filename(), type=unit.name))
		}
	)
)




