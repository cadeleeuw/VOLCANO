## alignment helper object, translates ordered alelle pairs into numeric code indicating type (absolute value) and strand (sign)
## stores allele codes for reference to match other input
AlleleAligner = R6::R6Class("AlleleAligner",
	private = list(
		alleles = c("A", "C", "G", "T", "a", "c", "t", "g"), allele.pairs = list(AA=NA, CC=NA, GG=NA, TT=NA, 	AC=1, CA=-1, TG=1, GT=-1,		AG=2, GA=-2, TC=2, CT=-2,		AT=11, TA=11, CG=12, GC=12),
		pair.index = NULL, ## data.frame with columns allele1, allele2, code,
		ref.alignment = NULL,  ## allele code vector for input, to compare against

		## create allele pair index from allele string codes
		map.alleles = function(allele1.str, allele2.str) {
			if (length(allele1.str) != length(allele2.str)) input.error("allele string vectors do not have the same length")
			allele1.index = match(allele1.str, private$alleles)
			allele2.index = match(allele2.str, private$alleles)

			return(private$pair.index$code[allele1.index + (allele2.index-1)*length(private$alleles)])
		},

		compare.pairs = function(allele.pair1, allele.pair2) {
			res = rep(NA, length(allele.pair1))
			res[allele.pair1 == -allele.pair2] = -1
			res[allele.pair1 == allele.pair2] = 1
			return(res)
		}
	),
	public = list(
		initialize = function(allele1.str, allele2.str) {
			private$pair.index = data.frame(allele1=rep(private$alleles, times=length(private$alleles)), allele2=rep(private$alleles, each=length(private$alleles)), stringsAsFactors=F)

			private$pair.index$code = unlist(private$allele.pairs[paste0(toupper(private$pair.index$allele1), toupper(private$pair.index$allele2))])
			private$pair.index$code[private$pair.index$code > 10] = NA	# setting strand ambiguous alleles to NA

			private$ref.alignment = private$map.alleles(allele1.str, allele2.str)
		},

		## return values: NA denotes mismatch (including invalid allele pair codes), -1 denotes matching but flipped, 1 denotes matching and same direction
		compare = function(allele1.str, allele2.str, index=NULL) {
			ref.pairs = if (!is.null(index)) private$ref.alignment[index] else private$ref.alignment
			if (length(allele1.str) != length(ref.pairs)) input.error("inconsistent input to function compare()")
			return(private$compare.pairs(ref.pairs, private$map.alleles(allele1.str, allele2.str)))
		}
	)
)


## wrapper for internal IDs for LD reference data; IDs are a numeric index based on the order SNPs are stored in the reference data
## NB: indices are kept unique and sorted
SNPindex = R6::R6Class("SNPindex",
	private = list(internal.ids = c(), hash.code = NULL),
	public = list(
		initialize = function(internal.ids) {
			if (length(internal.ids) > 0) {
				if (is.numeric(internal.ids)) private$internal.ids = sort(unique(internal.ids))
				if (length(private$internal.ids) == 0 || private$internal.ids[1] < 0) fatal.error("invalid input to SNPindex object")
			}
			private$hash.code = rlang::hash(private$internal.ids)
		},

		id = function() {return(private$hash.code)},
		snps = function() {return(private$internal.ids)},
		size = function() {return(length(private$internal.ids))},

		equals = function(internal.ids) {
			if (!has.type(internal.ids, "SNPindex")) {
				return(self$size() == length(internal.ids) && all(self$snps() == internal.ids))
			} else return(self$id() == internal.ids$id())
		},

		contains = function(internal.ids) {
			if (has.type(internal.ids, "SNPindex")) internal.ids = internal.ids$snps()
			return(length(internal.ids) <= self$size() && all(internal.ids %in% self$snps()))
		},

		intersection = function(internal.ids) {return(SNPindex$new(private$internal.ids[self$extract.index(internal.ids)]))},

		## returns numeric index of where this SNPindex matches input, to be used to subset (objects aligned with) internal.ids to SNPs in this SNPindex, in the same order as in this SNPindex
		extract.index = function(internal.ids) {
			if (has.type(internal.ids, "SNPindex")) out = match(self$snps(), internal.ids$snps())
			else out = match(self$snps(), internal.ids)
			return(out[!is.na(out)])
		}
	)
)



#########################################



## base class for different LD reference input formats
LDreferenceData = R6::R6Class("LDreferenceData",
	private = list(
		aligner = NULL, ## AlleleAligner object
		alignment.index = NULL, ## vector of allele pair code for each SNP
		snp.info = NULL, ## data.frame with columns internal.id, snp.id, chromosome, position, allele1, allele2
		max.id = NA,

		## at present, internal IDs match numeric index into snp.info, as codified here
		from.index = function(index) {return(index)},  ## convert snp.info index to internal ID
		to.index = function(internal.ids) {return(if (has.type(internal.ids, "SNPindex")) internal.ids$snps() else internal.ids)}  ## convert internal ID to snp.info index
	),
	public = list(
		initialize = function() {fatal.error("cannot initialize base LDreference objects")},

		no.snps = function() {return(nrow(private$snp.info))},
		get.snp.info = function(select.ids=NULL) {return(if (!is.null(select.ids)) private$snp.info[private$to.index(select.ids),] else private$snp.info)},

		## sum.stats: data.frame with snp.id, allele1 and allele2 columns
		## updates sum.stats data.frame: filtered to SNPs in reference data, in same order as reference, and adds direction column
		## - adds direction column multiply with statistics, with a value of -1, 0 or 1; set to NA if SNP cannot be aligned
		## - removes SNPs not in reference, reorders remaining to match order in reference
		## - removes the snp.id, allele1 and allele2 columns, and adds an internal.id column corresponding to rows of the reference data SNP info
		align.sumstats = function(sum.stats) {
			sum.stats$index = match(sum.stats$snp.id, private$snp.info$snp.id)
			sum.stats = sum.stats[!is.na(sum.stats$index),]

			if (nrow(sum.stats) > 0) {
				sum.stats$direction = private$aligner$compare(sum.stats$allele1, sum.stats$allele2, index=sum.stats$index)
				if (nrow(sum.stats) > 1) sum.stats = sum.stats[order(sum.stats$index),]
				sum.stats = subset(sum.stats, select=-c(snp.id, allele1, allele2))
			}

			sum.stats = sum.stats[c("index", names(sum.stats)[names(sum.stats) != "index"])]
			sum.stats$index = private$from.index(sum.stats$index); names(sum.stats)[1] = "internal.id"
			rownames(sum.stats) = NULL
			return(sum.stats)
		},

		get.block = function(subset) {undefined.error("get.block", class(self)[1])}
	)
)


## interface to PLINK data file set; loads .bim file on initialization, SNP IDs are stored in lower case
## call load.data to obtain raw genotype data for specified selection of SNPs
PlinkData = R6::R6Class("PlinkData",
	inherit = LDreferenceData,
	private = list(
		data.prefix = NULL,
		indiv.info = NULL,  ## data.frame with columns indiv.id, family.id, sex, phenotype

		error = function(...) {input.error(..., .label=list("PLINK data prefix"=private$data.prefix))},
		data.file = function(ext) {return(paste0(private$data.prefix, ".", ext))},

		## internal ID represents index into .bed file
		read.bim = function() {
			snp.info = data.table::fread(private$data.file("bim"), data.table=F, showProgress=F)
			if (ncol(snp.info) != 6) private$error("invalid .bim file, incorrect number of columns")

			snp.info = cbind(1:nrow(snp.info), snp.info[,c(2,1,4:6)])
			names(snp.info) = c("internal.id", "snp.id", "chromosome", "position", "allele1", "allele2")
			snp.info$snp.id = tolower(snp.info$snp.id)

			if (is.character(snp.info$chromosome)) {
				snp.info$chromosome[snp.info$chromosome == "x" | snp.info$chromosome == "X"] = 23
				snp.info$chromosome = suppressWarnings(as.numeric(snp.info$chromosome))
				snp.info$chromosome[is.na(snp.info$chromosome)] = -1
			}

			if (any(duplicated(snp.info$snp.id))) private$error("duplicate SNP IDs in .bim file")

			private$snp.info = snp.info; private$max.id = nrow(snp.info)
			private$aligner = AlleleAligner$new(snp.info$allele1, snp.info$allele2)
		},

		read.fam = function() {
			indiv = data.table::fread(private$data.file("fam"), data.table=F, showProgress=F)
			if (ncol(indiv) != 6) private$error("invalid format for .fam file")
			if (nrow(indiv) < globals$ref.minimum.N) private$error("sample size for genotype reference data is smaller than minimum of ", globals$ref.minimum.N)
			indiv = indiv[,c(2,1,5,6)]; names(indiv) = c("indiv.id", "family.id", "sex", "phenotype")

			indiv$sex = suppressWarnings(as.numeric(indiv$sex))
			indiv$phenotype = suppressWarnings(as.numeric(indiv$phenotype))

			indiv$sex[!(indiv$sex %in% 1:2)] = NA
			if (all(indiv$phenotype[!is.na(indiv$phenotype)] %in% c(-9,0,1,2))) indiv$phenotype[!(indiv$phenotype %in% 1:2)] = NA  ## set missing values for binary phenotype

			private$indiv.info = indiv
		},

		read.bed = function(subset=NULL) {
			if (!is.null(subset)) {
				if (has.type(subset, "SNPindex")) internal.ids = subset$snps()
				else internal.ids = private$from.index(which(private$snp.info$snp.id %in% tolower(subset)))

				if (length(internal.ids) == 0) return(list(data=data.frame(1:nrow(private$indiv.info))[,-1], internal.id=numeric(0)))
				if (tail(internal.ids, 1) > private$max.id) fatal.error("encountered invalid internal IDs", .label=list("PLINK data prefix"=private$data.prefix))
			} else internal.ids = NULL
			indiv.ids = private$indiv.info$indiv.id
			snp.ids = self$get.snp.info(internal.ids)$snp.id

			invisible(getNamespace("snpStats"))
    	out = as(.Call("readbed", private$data.file("bed"), indiv.ids, snp.ids, NULL, internal.ids, PACKAGE = "snpStats"), "numeric")
    	rownames(out) = NULL; colnames(out) = NULL
    	return(list(data = out, internal.ids = internal.ids))
		}
	),
	public = list(
		initialize = function(prefix)	{
			private$data.prefix = prefix

			ext = c("bed", "bim", "fam")
			exists = file.exists(private$data.file(ext))
			if (!all(exists)) private$error("missing ", paste0(".", ext[!exists], collapse=", "), " input file(s)")

			private$read.bim()
			private$read.fam()
		},

		abbreviate = function() {return(paste0("PLINK genotype data (", private$data.prefix, ")"))},
		print = function(...) {
			printer = ObjectPrinter$new("LD Reference (PlinkData)")$add.line("data prefix:", private$data.prefix)
			printer$add.line("data size:", self$no.indiv(), "individuals,", self$no.snps(), "variants")
			cat(printer$to.string())
		},

		no.indiv = function() {return(nrow(private$indiv.info))},
		get.indiv.info = function() {return(private$indiv.info)},

		## subset: list of SNP IDs or a SNPindex object; if NULL, returns all SNPs in data
		get.block = function(subset, load.individuals=F, settings.override=list()) {
			data = private$read.bed(subset)
			return(LDblockRaw$new(self, data$internal.id, data$data, indiv.info=if (load.individuals) private$indiv.info, settings.override=settings.override))
		}
	)
)



