## composite of LD reference and any number of summary statistics inputs
## summary statistics are loaded from file according to input specification, and are aligned to the LD reference
## phenotypes are loaded on demand, or when data is explicitly accessed; phenotype info objects for all are stored in 'available'
## sampling correlations are assumed to be zero unless explicitly specified using set.sampling.correlation function
##   if partial sampling correlation matrix is set, correlations with excluded phenotypes are assumed to be zero
##   input matrix needs to have phenotype names set as column names
## data can be set to load only a subset of chromosomes
##   this only affects data loading for phenotypes with input files split by chromosome, whole-genome input files are loaded in their entirety
DataSet = R6::R6Class("DataSet",
	private = list(
		ld.reference = NULL, sum.stats = list(), sampling.correlation = NULL,
		available = list(),	chromosomes = "all",

		## check validity of phenotypes, match partial phenotype names against available phenotypes, and remove duplicates/expand wildcards
		check.phenotypes = function(phenotypes, force.load=F) {
			if (length(phenotypes) > 0) {
				phenotypes = check.phenotypes(phenotypes, names(private$available), label="input data")
				if (force.load) self$load(phenotypes)
				return(phenotypes)
			} else {
				if (is.null(phenotypes)) return(NULL)
				else return(character(0))
			}
		},

		load.phenotypes = function(phenotypes) {
			for (phenotype in private$check.phenotypes(phenotypes)) {
				if (!(phenotype %in% names(private$sum.stats))) {
					pheno.info = private$available[[phenotype]]

					log.message("reading summary statistics for phenotype '", pheno.info$get.name(), "' (", toupper(pheno.info$get.traits()), ") from file ",  pheno.info$get.filename())

					if (pheno.info$by.chromosome()) {
						if (!private$chromosomes[1] == "all") log.message("loading input for chromosome(s) ", print.chromosomes(private$chromosomes), .indent=2)
						chromosomes = validate.chromosomes(private$chromosomes, condense=F)
						available = pheno.info$has.chromosomes(chromosomes, aggregate=F)
						if (!all(available)) {
							if (!any(available)) input.error("input files for specified chromosomes do not exist")
							if (private$chromosomes[1] != "all" || any(chromosomes[!available] != 23)) throw.warning("missing input files for chromosome(s) ", print.chromosomes(chromosomes[!available]))
						}
					}

					config = pheno.info$get.configuration(add.columns=T)
					log.message("using columns ", paste0(config$column, " (", config$parameter, ")", collapse=", "), .indent=2)

					InputType = if (pheno.info$is("localized")) LocalizedSummaryStatistics else SummaryStatistics
					statistics = InputType$new(pheno.info, private$ld.reference, chromosomes=private$chromosomes)
					private$check.counts(statistics$get.counts())
					if (pheno.info$is("localized")) {
						unit.counts = statistics$get.counts()$units
						log.message("localized unit counts: ", paste0(unlist(unit.counts), " of type '", names(unit.counts), "'", collapse=", "), .indent=2)
					}

					private$sum.stats[[pheno.info$get.name()]] = statistics
				}
			}
		},

		check.counts = function(...) {
			counts = flatten.arglist(...); prop.thresh = globals$variant.match.warning.proportion
			if (counts$available / counts$valid < prop.thresh) {
				if (counts$available == 0) input.error("no variants could be matched to reference data; please make sure the SNP ID format matches the reference data")
				else throw.warning("less than ", 100*prop.thresh, "% of variants matched to reference data; please make sure the SNP ID format matches the reference data")
			}
			if (counts$aligned / counts$available < prop.thresh) {
				if (counts$aligned == 0) input.error("none of the variants matched to reference data could be aligned; please make sure the allele columns are correct")
				else throw.warning("less than ", 100*prop.thresh, "% of variants matched to reference data could be aligned; please make sure the allele columns are correct")
			}
			log.message("read ", counts$rows, " variants, of which ", counts$aligned, " valid and aligned to reference", .indent=2)
		}
	),
	public = list(
		initialize = function(ld) {private$ld.reference = ld},

		print = function(...) {
			printer = ObjectPrinter$new("DataSet")
			printer$add.parameter("LD reference", private$ld.reference$abbreviate())
			printer$add.parameter("loaded phenotypes", length(private$sum.stats))
			for (ph in names(private$sum.stats)) printer$add.line(private$available[[ph]]$abbreviate(), indent=2)

			unloaded = names(private$available)[!(names(private$available) %in% names(private$sum.stats))]
			if (length(unloaded) > 0) {
				printer$add.parameter("available phenotypes", length(unloaded))
				for (ph in unloaded) printer$add.line(private$available[[ph]]$abbreviate(), indent=2)
			}
			cat(printer$to.string())
		},

		no.pheno = function() {return(list(available=length(private$available), loaded=length(private$sum.stats)))},
		phenotypes = function(loaded.only=T) {return(if (loaded.only) names(private$sum.stats) else names(private$available))},

		get.reference = function() {return(private$ld.reference)},
		get.sumstats = function(phenotypes=NULL) {
			if (!is.null(phenotypes)) {
				phenotypes = private$check.phenotypes(phenotypes, force.load=T)
				return(private$sum.stats[phenotypes])
			} else return(private$sum.stats)
		},
		get.sampling.correlation = function(phenotypes=NULL) { ## fills out to full matrix matching phenotypes argument, or full list of phenotypes in data if NULL
			if (is.null(phenotypes)) phenotypes = names(private$sum.stats)
			else phenotypes = private$check.phenotypes(phenotypes)

			if (length(phenotypes) > 0) {
				out = diag(length(phenotypes))
				rownames(out) = colnames(out) = phenotypes
				if (!is.null(private$sampling.correlation)) {
					mapping = match(phenotypes, colnames(private$sampling.correlation))
					index.out = which(!is.na(mapping)); index.in = mapping[!is.na(mapping)]
					if (length(index.out) > 0) out[index.out,index.out] = private$sampling.correlation[index.in,index.in]
				}
				return(out)
			} else return(matrix(1, 0, 0))
		},

		get.interface = function(phenotypes=NULL, settings.override=list()) {
			if (is.null(phenotypes)) {
				if (length(private$sum.stats) == 0) input.error("no phenotypes have been loaded")
				phenotypes = names(private$sum.stats)
			} else phenotypes = private$check.phenotypes(phenotypes, force.load=T)
			if (length(phenotypes) == 0) input.error("list of phenotypes is empty")
			return(DataInterface$new(self, phenotypes, settings.override=settings.override))
		},

		add.phenotype = function(pheno.info) {
			if (pheno.info$get.name() %in% names(private$available)) input.error("duplicate phenotype label '", pheno.info$get.name(), "'; please rename")
			private$available[[pheno.info$get.name()]] = pheno.info
			invisible(self)
		},

		## chromosomes argument is one or a vector of numeric chromosome codes (chromosome X can be denoted as 'X' or 23), or 'all'
		set.chromosomes = function(chromosomes="all") {
			chromosomes = validate.chromosomes(chromosomes)
			if (length(chromosomes) != length(private$chromosomes) || !all(chromosomes == private$chromosomes)) {
				reload = private$chromosomes[1] != "all" && !all(chromosomes %in% private$chromosomes)
				private$chromosomes = chromosomes
				if (reload) {
					for (ph in names(private$sum.stats)) {
						if (private$available[[ph]]$by.chromosome()) {
							private$sum.stats[[ph]] = NULL
							private$load.phenotypes(ph)
						}
					}
				}
			}
			invisible(self)
		},

		set.sampling.correlation = function(covar) {
			private$sampling.correlation = cov2cor(check.covariance(covar, allow.NA=T, label="sampling correlation matrix"))
			invisible(self)
		},

		## load specified phenotypes, or all if none specified
		load = function(phenotypes=NULL) {
			if (is.null(phenotypes)) phenotypes = names(private$available)
			else phenotypes = private$check.phenotypes(phenotypes)

			private$load.phenotypes(phenotypes)
			invisible(self)
		},

		## unload specified phenotypes, or all if none specified
		## NB: memory may not be freed immediately if object is still referenced elsewhere
		unload = function(phenotypes=NULL) {
			if (!is.null(phenotypes)) {
				phenotypes = private$check.phenotypes(phenotypes)
				private$sum.stats[phenotypes] = NULL
			} else private$sum.stats = list()
			invisible(self)
		},

		check.overlap = function() {
			pattern.matrix = NULL
			if (length(private$sum.stats) > 1) {
				pattern.matrix = matrix(0, nrow=private$ld.reference$no.snps(), ncol=length(private$sum.stats))
				for (i in 1:length(private$sum.stats)) pattern.matrix[private$sum.stats[[i]]$get.statistics()$internal.id,i] = 1
				colnames(pattern.matrix) = names(private$sum.stats)
				if (private$chromosomes[1] != "all") pattern.matrix = pattern.matrix[private$ld.reference$get.snp.info()$chromosome %in% private$chromosomes,]

				variants.used = sum(apply(pattern.matrix == 1, 1, any))
				overlap.global = sum(apply(pattern.matrix == 1, 1, all))
				overlap.proportion = overlap.global/variants.used
				pairwise.overlap = t(pattern.matrix) %*% pattern.matrix
				log.message("data contains ", variants.used, " variants aligned for at least one phenotype, of which ", overlap.global, " variants shared across all phenotypes (", round(100*overlap.proportion, 2), "%)")

				prop.thresh = globals$variant.overlap.warning.proportion
				if (ncol(pattern.matrix) > 2) {
					pairwise.proportion = pairwise.overlap / (outer(diag(pairwise.overlap), diag(pairwise.overlap), FUN="+") - pairwise.overlap)

					proportion.range = range(pairwise.proportion[lower.tri(pairwise.proportion)])
					log.message("pairwise overlap proportions are between ", round(100*proportion.range[1], 2), "% and ", round(100*proportion.range[2],2), "% ", .indent=2)
					if (proportion.range[1] < prop.thresh) throw.warning(ifelse(proportion.range[2] < prop.thresh, "all", "some"), " phenotype pairs have variant overlap proportion of less than ", 100*prop.thresh, "%; analysis may be unstable")
				} else {
					if (overlap.proportion < prop.thresh) throw.warning("variant overlap proportion is less than ", 100*prop.thresh, "%; analysis may be unstable")
				}
				if (max(pairwise.overlap[lower.tri(pairwise.overlap)]) <= globals$variant.overlap.failure.count) input.error("variant overlap between phenotypes is negligible")
			}

			invisible(pattern.matrix)
		}
	)
)

## provides an interface to retrieve LD and summary statistics for specific phenotypes and locus
## phenotypes are set on initialization, the locus and phenotype subunits are set and updated later
DataInterface = R6::R6Class("DataInterface",
	private = list(
		settings = list(
			prune.threshold = globals$pca.prune.threshold,
			minimum.snps = globals$locus.minimum.snps,
			minimum.components = globals$locus.minimum.components,
			correlation.bound = globals$sampling.correlation.truncation,
			load.individuals = FALSE
		),
		ld.reference = NULL, ## LDreferenceData object
		sum.stats = NULL, ## named list of SummaryStatistics object per phenotype
		sampling.correlation = NULL, ## row/col named matrix of sampling correlations
		units = list(
			unit.index = list(), ## named list of phenotypes per unit
			pheno.index = list(), ## named list of units per phenotype
			values = list() ## named list of available unit values per unit name
		),
		current = list(
			annotation = NULL, ## LocusAnnotation object
			locus = NULL, ## LocusDefinition object
			units = list(), ## named list of specified unit values per unit name (if NULL, use all available in locus)
			data = list() ## list containing cached summary statistics loaded by fetch.statistics
		),

		validate.settings = function(...) {
			private$settings = validate.settings(flatten.arglist(...), private$settings, self, merge=T)
			for (param in c("minimum.snps", "minimum.components")) {
				if (param %in% names(private$settings) && private$settings[[param]] < globals$locus.minimum.lower.bound) input.error("parameter '", param, "' cannot be lower than ", globals$locus.minimum.lower.bound)
			}
		},

		fetch.statistics = function() {
			if (is.null(private$current$locus)) input.error("no locus has been set")

			if (length(private$current$data) == 0) {
				statistics = list(); counts = list()
				snp.index = private$current$locus$snp.index(private$ld.reference$get.snp.info())
				for (ph in names(private$sum.stats)) {
					if (private$sum.stats[[ph]]$is.localized()) {
						unit.values = lapply(private$units$pheno.index[[ph]], function(u) {private$current$units[[u]]})
						names(unit.values) = private$units$pheno.index[[ph]]
						locus.stats = private$sum.stats[[ph]]$get.statistics(unit.values, select.ids=snp.index)
						if (nrow(locus.stats) > 0) {
							for (u in names(unit.values)[sapply(unit.values, is.null)]) unit.values[[u]] = unique(locus.stats[[u]])

							value.sets = combinations(unit.values); no.snps = rep(0, nrow(value.sets)); curr.stats = list()
							for (i in 1:nrow(value.sets)) {
								select = rep(T, nrow(locus.stats))
								for (j in 1:ncol(value.sets)) select = select & locus.stats[[names(value.sets)[j]]] == value.sets[i,j]
								no.snps[i] = sum(select)
								if (sum(select) >= private$settings$minimum.snps) {
									curr.info = private$sum.stats[[ph]]$get.info()$instantiate.units(value.sets[i,,drop=F])
									statistics[[curr.info$get.name()]] = list(phenotype=ph, info=curr.info, data=locus.stats[select,!(names(locus.stats) %in% names(value.sets))])
								}
							}
							counts[[ph]] = cbind(value.sets, no.snps = no.snps)
						}
					} else {
						curr.stats = private$sum.stats[[ph]]$get.statistics(select.ids=snp.index)
						if (nrow(curr.stats) >= private$settings$minimum.snps) statistics[[ph]] = list(phenotype=ph, info=private$sum.stats[[ph]]$get.info(), data=curr.stats)
					}
				}
				private$current$data = list(statistics=statistics, snp.index=snp.index, counts=counts)
			}
			return(private$current$data)
		}
	),
	public = list(
		initialize = function(data, phenotypes=NULL, settings.override=list()) {
			private$ld.reference = data$get.reference()
			private$sum.stats = data$get.sumstats(phenotypes)
			private$sampling.correlation = data$get.sampling.correlation(phenotypes)
			if (length(settings.override) > 0) private$validate.settings(settings.override)

			if (any(duplicated(names(private$sum.stats)))) input.error("duplicate phenotypes selected in DataInterface")

			localized = sapply(private$sum.stats, function(ss) {ss$is.localized()})
			if (any(localized)) {
				for (ph in names(private$sum.stats[localized])) {
					curr = private$sum.stats[[ph]]$get.units()
					private$units$pheno.index[[ph]] = curr
					for (u in curr) {
						private$units$unit.index[[u]] = c(private$units$unit.index[[u]], ph)
						private$units$values[[u]] = unique(c(private$units$values[[u]], private$sum.stats[[ph]]$get.unit.values(u)[,1]))
					}
				}
			}
		},

		print = function(...) {
			printer = ObjectPrinter$new("DataInterface")
			printer$add.parameter("LD reference", private$ld.reference$abbreviate())$add.parameter("phenotypes", length(private$sum.stats))
			for (i in 1:length(private$sum.stats)) printer$add.line(private$sum.stats[[i]]$get.info()$abbreviate(), indent=2)

			if (!is.null(private$current$locus)) printer$add.line("current locus:", private$current$locus$get.name())
			else printer$add.line("current locus: NA")
			if (length(names(private$units$unit.index)) > 0) {
				printer$add.line("current units:")
			 	for (unit in names(private$units$unit.index)) printer$add.parameter(unit, ifelse(is.null(private$current$units[[unit]]), "[ALL]", paste(private$current$units[[unit]], collapse=", ")), indent=2)
			}

			if (length(private$settings) > 0) {
				printer$add.line("settings:")
				for (s in names(private$settings)) printer$add.parameter(s, private$settings[s], indent=2)
			}
			cat(printer$to.string())
		},

		set.annotation = function(annotation=NULL) {
			if (!is.null(annotation)) {
				check.types(annotation="LocusAnnotation")
				self$set.locus(NULL)
			}
			private$current$annotation = annotation
			invisible(self)
		},
		create.annotation = function(unit.name, phenotypes=NULL, merge.mode=c("union", "intersection")) {merge.mode = match.arg(merge.mode)
			if (length(unit.name) > 1) input.error("cannot create multi-unit annotation")
			if (!(unit.name %in% names(private$units$unit.index))) input.error("unknown unit ", unit.name)

			if (length(phenotypes) > 0) {
				phenotypes = check.phenotypes(phenotypes, names(private$sum.stats), label="data interface")
				matched = phenotypes %in% names(private$units$unit.index)
				if (!all(matched)) input.error("unit '", unit.name, "' not available for %phenotype% ", items.or=phenotypes[!matched])
			} else phenotypes = private$units$unit.index[[unit.name]]

			annot = private$sum.stats[[phenotypes[1]]]$create.annotation(unit.name)
			if (length(phenotypes) > 1) {
				for (p in 2:length(phenotypes)) annot = annot$merge(private$sum.stats[[phenotypes[p]]]$create.annotation(unit.name), merge.mode=merge.mode)
			}

			self$set.annotation(annot)
			invisible(self)
		},

		## advance to next locus in internal annotation, return FALSE if no further remaining
		next.locus = function() {
			if (is.null(private$current$annotation)) input.error("no annotation has been specified for data interface")
			next.id = ifelse(!is.null(private$current$locus), private$current$locus$get.id() + 1, 1)
			if (next.id <= private$current$annotation$size()) {self$set.locus(next.id); invisible(TRUE)	}
			else {self$set.locus(NULL); invisible(FALSE)}
		},

		## provide LocusDefinition object obtained from a LocusAnnotation (output of read.loci() function)
		## if internal annotation is set, locus.info can also be string or numeric index, to retrieve corresponding locus from internal annotation
		## if detect.units is true, will automatically fill in unit values for units with name matching the LocusDefinition type
		set.locus = function(locus.info, detect.units=T) {
			if (!is.null(locus.info) && !has.type(locus.info, "LocusDefinition")) {
				if (!is.null(private$current$annotation)) {
					locus.info = private$current$annotation$get.locus(locus.info)
				} else input.error("argument 'locus.info' must have type 'LocusDefinition' unless annotation has been specified for data interface")
			}
			private$current$locus = locus.info
			private$current$units = list()
			private$current$data = list()

			if (!is.null(locus.info) && detect.units && length(names(private$units$unit.index)) > 0 && !("unknown" %in% locus.info$get.type())) {
				values = list()
				for (type in locus.info$get.type()[locus.info$get.type() %in% names(private$units$unit.index)]) values[[type]] = locus.info$get.name()
				if (length(values) > 0) self$set.units(values)
			}

			invisible(self)
		},

		## specify unit values as (list of) named arguments (can be vectors of multiple values)
		set.units = function(...) {
			units = flatten.arglist(...)
			if (is.null(names(units))) input.error("unit names not specified")
			unknown = !(names(units) %in% names(private$units$unit.index))
			if (any(unknown)) input.error("unknown %unit% ", items.and=names(units)[unknown])

			private$current$data = list()  ## clear cached statistics
			for (u in names(units)) {
				values = unique(units[[u]])
				valid = values %in% private$units$values[[u]]
				if (!all(valid)) input.error("unknown %value% ", items.and=values[!valid], " for unit '", u, "'")
				private$current$units[[u]] = values
			}
			invisible(self)
		},

		available.units = function() {
			units = names(private$units$unit.index); out = NULL
			if (length(units) > 0) {
				counts = lapply(private$fetch.statistics()$counts, function(curr) {curr[units[!(units %in% names(curr))]] = NA; return(curr)})
				for (ph in names(counts)) out = rbind(out, data.frame(phenotype=ph, counts[[ph]][,c(units, "no.snps")], stringsAsFactors=F))
			}
			return(out)
		},

		## retrieves and processes LD and summary statistics, returns LocusData object
		## if discard.empty=T, discard phenotypes with insufficient SNPs in input and return reduced data object, rather than terminating
		## if compute.marginal=T, process immediately into ProcessedLocus object, otherwise return LocusData object
		process = function(discard.empty=T, compute.marginal=T) {
			data = private$fetch.statistics()

			## create SNP index and LD block, and verify configuration
			truncated = sapply(data$statistics, function(s) {s$info$is("truncated")})
			available = matrix(F, nrow=data$snp.index$size(), ncol=length(data$statistics))
			for (i in seq_along(data$statistics)) available[,i] = data$snp.index$snps() %in% data$statistics[[i]]$data$internal.id
			if (any(truncated)) {
				if (!all(truncated)) {
					if (sum(!truncated) > 1) keep = apply(available[,!truncated], 1, all)
					else keep = available[,!truncated]
					data$statistics = data$statistics[apply(available[keep,], 2, sum) >= private$settings$minimum.snps]
				} else keep = apply(available, 1, any)
			} else keep = apply(available, 1, all)
			data$snp.index = SNPindex$new(data$snp.index$snps()[keep])

			empty.pheno = !(names(private$sum.stats) %in% sapply(data$statistics, function(s) {s$phenotype}))
			if (all(empty.pheno) || data$snp.index$size() < private$settings$minimum.snps) data.error("current configuration contains fewer than ", private$settings$minimum.snps, " SNPs, cannot load data")
			if (any(empty.pheno)) {
				msg = list("locus contains summary statistics for fewer than ", private$settings$minimum.snps, " SNPs for %phenotype% ", items.and=names(private$sum.stats)[empty.pheno])
				if (discard.empty) {msg = c(msg, "; discarding %phenotype%"); throw.warning(msg) } else data.error(msg)
			}

			sample.sizes = unlist(lapply(data$statistics, function(s) {ifelse("sample.size" %in% names(s$data), mean(s$data$sample.size), s$info$get.global("sample.size"))}))
			if (any(is.na(sample.sizes))) fatal.error("sample size information is not available for all phenotypes")
			max.components = max(globals$local.minimum.components, globals$component.sample.ratio * min(sample.sizes))

			ld.block = private$ld.reference$get.block(data$snp.index, load.individuals=private$settings$load.individuals, settings.override=list(prune.threshold=private$settings$prune.threshold, max.components=max.components))
			if (ld.block$no.components() < private$settings$minimum.components) data.error("fewer than ", private$settings$minimum.components, " genetic principal components in locus genotype data")


			## filter statistics to SNP selection, and create SumstatBlocks
			stat.blocks = list(); truncate.failed = c()
			for (ph.out in names(data$statistics)) {
				curr.stats = data$statistics[[ph.out]]$data; curr.ld = ld.block
				if (data$statistics[[ph.out]]$info$is("truncated")) {
					curr.ld = ld.block$intersection(curr.stats$internal.id)
					if (curr.ld$no.components() < private$settings$minimum.components) {truncate.failed = c(truncate.failed, ph.out); next}
					curr.stats = tryCatch(curr.ld$align.data(curr.stats, mode="truncate"), error=function(err) {private$fatal("invalid input to SumstatBlock object, ", err)})
					if (nrow(curr.stats) < private$settings$minimum.snps) {truncate.failed = c(truncate.failed, ph.out); next}
				} else curr.stats = tryCatch(curr.ld$align.data(curr.stats, mode="truncate"), error=function(err) {private$fatal("invalid input to SumstatBlock object, ", err)})

				BlockType = if (data$statistics[[ph.out]]$info$is("binary")) BinarySumstatBlock else ContinuousSumstatBlock
				stat.blocks[[ph.out]] = BlockType$new(curr.ld, curr.stats, data$statistics[[ph.out]]$info)
			}

			if (!is.null(truncate.failed)) {
				failed.index = names(data$statistics) %in% truncate.failed
				failed.pheno = unique(sapply(data$statistics[failed.index], function(s) {s$phenotype}))
				empty.pheno = !(failed.pheno %in% sapply(data$statistics[!failed.index], function(s) {s$phenotype}))
				if (any(empty.pheno)) {
					msg = list("not enough SNPs and/or genetic principal components in locus genotype data for truncated %phenotype% ", items.and=failed.pheno[empty.pheno])
					if (discard.empty) {msg = c(msg, "; discarding %phenotype%"); throw.warning(msg) } else data.error(msg)

				}
				data$statistics = data$statistics[!failed.index]
				if (length(data$statistics) == 0) data.error("after filtering on minimum number of SNPs and genetic principal components, no phenotypes left to analyse")
			}

			## create sampling correlation matrix
			margin = sapply(data$statistics, function(s) {s$phenotype})
			if (length(margin) > 1) {
				corr.matrix = private$sampling.correlation[margin, margin]
				corr.matrix[outer(margin, margin, FUN="==") & row(corr.matrix) != col(corr.matrix)] = NA
				truncate = !is.na(corr.matrix) & row(corr.matrix) != col(corr.matrix) & abs(corr.matrix) < private$settings$correlation.bound
				corr.matrix[truncate] = 0
			} else corr.matrix = matrix(1,1,1)
			colnames(corr.matrix) = rownames(corr.matrix) = names(stat.blocks)

			locus.data = LocusData$new(private$current$locus, ld.block, stat.blocks, corr.matrix, settings.override=private$settings)
			if (compute.marginal) return(locus.data$process())
			else return(locus.data)
		},

		phenotypes = function() {return(names(private$sum.stats))},

		get.settings = function() {return(private$settings)},
		get.input = function() {return(list(ld.reference=private$ld.reference, sum.stats=private$sum.stats, sampling.correlation=private$sampling.correlation))},

		get.annotation = function() {return(private$current$annotation)},
		get.locus = function() {return(private$current$locus)},
		get.units = function() {return(list(units=private$units$unit.index, values=private$current$units))}
	)
)




