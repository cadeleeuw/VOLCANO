## LD and decomposition for block of SNPs
LDblock = R6::R6Class("LDblock",
	private = list(
		source.info = list(label=NULL),
		block.id = NULL,
		snp.index = NULL,  ## SNPindex object
		snp.info = NULL, ## data.frame with columns internal.id, snp.id, chromosome, position, allele1, allele2
		decomposition = list(no.components=NULL, vectors=NULL, values=NULL, components=NULL),  ## pruned decomposition of LD matrix
		settings = list(
			prune.threshold = globals$pca.prune.threshold,  ## proportion of variance to retain after decomposition
			max.components = Inf
		),

		get.printer = function() {
			printer = ObjectPrinter$new("LD Block")$parameter.list(source=private$source.info$label)
			no.comp = self$no.components()
			printer$add.parameter("no. variants", self$no.snps())$add.parameter("no. components", ifelse(is.na(no.comp), "unknown", no.comp))
			return(printer)
		},

		validate.settings = function(settings) {
			private$settings = validate.settings(settings, private$settings, self, merge=T)
			if (private$settings$prune.threshold > 1) private$settings$prune.threshold = private$settings$prune.threshold/100
			if (private$settings$max.components < globals$locus.minimum.lower.bound) private$settings$max.components = globals$locus.minimum.lower.bound
		},

		## subsets block to only those SNPs specified in internal.ids argument (discarding any not found in block)
		slice = function(internal.ids) {undefined.error("slice", class(self)[1])}
	),
	public = list(
		initialize = function(source, settings.override=list()) {
			private$source.info$label = source$abbreviate()
			if (length(settings.override) > 0) private$validate.settings(settings.override)
		},
		print = function(...) {cat(private$get.printer()$to.string())},

		id = function() {return(private$block.id)},
		equals = function(other) {return(class(self)[1] == class(other)[1] && self$id() == other$id())},

		no.snps = function() {return(private$snp.index$size())},
		no.components = function(force.compute=T) {return(ifelse(is.null(private$decomposition$no.components), NA, private$decomposition$no.components))},

		get.snp.index = function() {return(private$snp.index)},
		get.snp.info = function() {return(private$snp.info)},

		get.vectors = function() {return(private$decomposition$vectors)}, # Q
		get.values = function() {return(private$decomposition$values)}, # L
		get.root = function() {return(private$decomposition$vectors %*% diag(sqrt(private$decomposition$values)))}, ## Q %*% diag(sqrt(L)
		get.inverse.root = function() {return(private$decomposition$vectors %*% diag(sqrt(1/private$decomposition$values)))}, ## Q %*% diag(sqrt(1/L)

		## create LD block containing subset of SNPs based on internal.ids (either SNPindex object or vector of internal IDs)
		subset = function(internal.ids) {
			if (!has.type(internal.ids, "SNPindex")) internal.ids = SNPindex$new(internal.ids)
			if (!private$snp.index$contains(internal.ids)) input.error("LD block does not contain all requested SNPs, cannot retrieve subblock")
			return(self$clone()$.__enclos_env__$private$slice(internal.ids))
		},

		## create LD block containing subset of SNPs shared with internal.ids argument
		intersection = function(internal.ids) {return(self$clone()$.__enclos_env__$private$slice(internal.ids))},

		## check/align data in data.frame with internal.id column to LD block,
		## - mode=truncate checks that all SNPs are in data, discards excess in data and reorders if needed
		## - mode=sort checks that SNPs are the same, but reorders data if needed
		## - mode=strict checks that SNPs are the same and already in same order, error if not
		## if the function returns without error, the input data will be aligned to the LD block
		align.data = function(data, mode=c("truncate", "sort", "strict")) {mode = match.arg(mode)
			if (!is.data.frame(data) || !("internal.id" %in% names(data))) fatal.error("invalid input to align.data()")
			if (mode == "truncate") data = data[data$internal.id %in% private$snp.index$snps(),]
			if (mode %in% c("sort", "truncate")) data = data[order(data$internal.id),]
			if (!private$snp.index$equals(data$internal.id)) input.error("data is not aligned with LD block")
			return(data)
		}
	)
)

## LD block based on raw genotype input data
LDblockRaw = R6::R6Class("LDblockRaw",
	inherit = LDblock,
	private = list(
		data = NULL,  ## data.frame with standardized genotype data (missing values imputed to 0)
		indiv.info = NULL,  ## data.frame with indiv.id, family.id, sex, phenotype (if loaded)

		process.data = function(internal.ids, data) {
			if (ncol(data) == 0) input.error("loaded genotype data contains zero SNPs")
			if (nrow(data) < globals$ref.minimum.N) input.error("sample size for genotype reference data is smaller than minimum of ", globals$ref.minimum.N)
			if (ncol(data) != length(internal.ids)) input.error("size of genotype data does not match internal ID vector")

			data = sweep(data.matrix(data), 2, apply(data, 2, mean, na.rm=T), FUN="-")
			data[is.na(data)] = 0
			sd = apply(data, 2, sd, na.rm=T)
			if (any(sd == 0)) {
				data = data[,sd > 0]
				internal.ids = internal.ids[sd > 0]
				sd = sd[sd > 0]
			}
			private$data = sweep(data, 2, sd, FUN="/")

			private$snp.index = SNPindex$new(internal.ids)
			if (!private$snp.index$equals(internal.ids)) private$data = private$data[,private$snp.index$extract.index(internal.ids)] ## make sure the SNPindex and data are fully aligned
			if (private$snp.index$size() != ncol(data)) fatal.error("unable to align SNP index to input data in LDblockRaw")
		},

		compute.decomposition = function() {
			if (ncol(private$data) > 0) {
				decomp = try(svd(private$data), silent=T)
				if (class(decomp) == "try-error") decomp = try(decomp(private$data), silent=T)
				if (class(decomp) == "try-error") {
					decomp = eigen(cov(private$data))
					vectors = decomp$vectors
					values = decomp$values
				} else {
					vectors = decomp$v
					values = decomp$d * decomp$d / (nrow(private$data) - 1) * sign(decomp$d)
					values = values[decomp$d > 0]
				}

				req.comp = which(cumsum(values / sum(values)) >= private$settings$prune.threshold)
				private$decomposition$no.components = min(req.comp, length(values), private$settings$max.components)
				private$decomposition$vectors = vectors[,1:private$decomposition$no.components]
				private$decomposition$values = values[1:private$decomposition$no.components]
			} else private$decomposition = list(no.components = 0, vectors = matrix(0, 0, 0), values = numeric(0))

			private$block.id = rlang::hash(c(private$snp.index$id(), private$decomposition$values))
		},

		## modifies object, should only be called directly after cloning
		## subsets block to only those SNPs specified in internal.ids argument (discarding any not found in block)
		slice = function(internal.ids) {
			if (!has.type(internal.ids, "SNPindex")) internal.ids = SNPindex$new(internal.ids)
			if (!private$snp.index$equals(internal.ids)) {
				subset = internal.ids$extract.index(private$snp.index)
				private$snp.info = private$snp.info[subset,]
				private$data = private$data[,subset]
				private$snp.index = SNPindex$new(private$snp.index$snps()[subset])
				private$compute.decomposition()
			}
			return(self)
		},

		get.printer = function() {return(super$get.printer()$set.name("LD Block (raw data)"))}
	),
	public = list(
		initialize = function(source, internal.ids, data, indiv.info=NULL, settings.override=list()) {
			super$initialize(source, settings.override)

			private$process.data(internal.ids, data)
			private$compute.decomposition()  ## also sets ID

			private$snp.info = subset(source$get.snp.info(private$snp.index), select=-internal.id)
			private$indiv.info = indiv.info
		},

		no.indiv = function() {return(nrow(private$data))},
		get.indiv.info = function() {return(private$indiv.info)},

		get.data = function() {return(private$data)},

		get.components = function() {
			if (is.null(private$decomposition$components)) private$decomposition$components = private$data %*% self$get.inverse.root()
			return(private$decomposition$components)
		},



#TODO add source ID check here; + move to external?
		## compute correlations of genetic components with those in other LDblock
		component.correlations = function(other) {
			if (other$no.indiv() != self$no.indiv()) input.error("LD blocks are incompatible, cannot compute correlations")
			return(t(self$get.components()) %*% other$get.components() / (self$no.indiv() - 1))
		},

		trim.components = function(drop=NULL, truncate=NULL) {
			curr.pcs = self$no.components()
			if (!is.null(drop)) {
				if (!is.null(truncate)) input.error("for function trim.components(), cannot use 'truncate' argument if 'drop' argument is specified")
				if (min(drop) <= 0 || max(drop) > curr.pcs) input.error("for function trim.components(), elements of 'drop' argument are out of bounds")
				drop = 1:curr.pcs %in% drop
			} else drop = 1:curr.pcs > truncate
			if (!is.null(drop)) {
				private$decomposition$vectors = private$decomposition$vectors[,!drop]
				private$decomposition$values = private$decomposition$values[!drop]
				private$decomposition$no.components = length(private$decomposition$values)
				if (!is.null(private$decomposition$components)) private$decomposition$components = private$decomposition$components[,!drop]
			}
			invisible(self)
		}
	)
)
