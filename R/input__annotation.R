## access interface for locus definitions
## either Range (based on genomic position, chr/start/stop) or List (list of SNP IDs)
## emits LocusDefinition objects for individual loci, to pass around

LocusAnnotation = R6::R6Class("LocusAnnotation",
	private = list(
		locus.names = NULL, source = NULL,
		type = NULL,  ## type of locus, if any (eg. gene)

		create.locus = function(index) {undefined.error("create.locus", class(self)[1])}
	),
	public = list(
		initialize = function(locus.names, source, type) {
			private$locus.names = locus.names; private$source = source; private$type = if (!("unknown" %in% NULL)) type
			if (any(duplicated(locus.names))) input.error("locus annotation input contains duplicate locus names")
		},
		print = function(...) {
			type.tag = switch(class(self)[1], RangeAnnotation="genomic coordinates", ListAnnotation="SNP lists", "UNKNOWN TYPE")
			printer = ObjectPrinter$new("Locus Annotation (", type.tag, ")")
			printer$parameter.list(source=private$source)
			type.tag = if (!is.null(private$type)) paste0("(", private$type, ")")
			printer$add.line("number of loci:", self$size(), type.tag)
			cat(printer$to.string())
		},

		size = function() {return(length(private$locus.names))},

		get.locus = function(id) {
			if (!is.numeric(id)) {
				index = match(id, private$locus.names)
				if (is.na(index)) input.error("locus ID '", id, "' not found")
				return(private$create.locus(index))
			} else {
				if (id <= 0 || id > length(private$locus.names)) input.error("locus index ", id, " is out of bounds")
				return(private$create.locus(id))
			}
		},
		get.names = function() {return(private$locus.names)},
		get.contents = function(index) {undefined.error("get.contents", class(self)[1])},

		get.type = function() {return(if (!is.null(private$type)) private$type else "unknown")},
		has.type = function(type) {
			if (is.null(type) || "unknown" %in% type) return(is.null(private$type))
			else return(all(type %in% private$type))
		},

		## output merge with other RangeAnnotation
		## - intersection: keep only shared loci, reduce loci to overlapping regions (Range) or shared SNPs (List)
		## - union: keep all loci, expand loci to combination of regions (Range; includes any area between regions, if disjoint!) or union of SNPs (List)
		merge = function(other, merge.mode=c("union", "intersection")) {
			check.types(other="LocusAnnotation", .func="merge", .class=class(self))
			if (!has.type(other, class(self)[1])) input.error("cannot perform annotation merge, annotation object types are incompatible")
			if (!(other$has.type(self$get.type()) && self$has.type(other$get.type()))) input.error("cannot merge annotations, locus types do not match")
			return(private$merge.internal(other, match.arg(merge.mode)))
		}
	)
)

## annotation specified in genomic coordinates
## input should be a data.frame with columns: locus, chromosome, start, stop
RangeAnnotation = R6::R6Class("RangeAnnotation",
	inherit = LocusAnnotation,
	private = list(
		ranges = NULL,

		create.locus = function(index) {return(RangeLocus$new(index, private$locus.names[index], private$ranges[index,], type=private$type, source=private$source))},

		merge.internal = function(other, merge.mode) {
			data1 = self$get.contents(); data2 = other$get.contents()

			merged = merge(data1, data2, by="locus", all=F); names(merged)[2:4] = c("chromosome", "start", "stop")
			if (nrow(merged) > 0) {
				if (any(merged$chromosome != merged$chromosome.y)) input.error("chromosomes do not match for all loci")
				start.higher = merged$start.y > merged$start; stop.higher = merged$stop.y > merged$stop
				if (merge.mode == "union") {merged$start[!start.higher] = merged$start.y[!start.higher]; merged$stop[stop.higher] = merged$stop.y[stop.higher]}
				if (merge.mode == "intersection") {merged$start[start.higher] = merged$start.y[start.higher]; merged$stop[!stop.higher] = merged$stop.y[!stop.higher]}
			}
			merged = merged[,1:4]

			if (merge.mode == "union") {
				merged = rbind(merged, data1[!(data1$locus %in% merged$locus),], data2[!(data2$locus %in% merged$locus),])
				merged = merged[order(merged$chromosome, merged$start, merged$stop),]
			}
			return(RangeAnnotation$new(merged, source=paste0("merged (", merge.mode, ")"), type=private$type))
		}
	),
	public = list(
		initialize = function(input, source=NULL, type=NULL) {
			names(input) = tolower(names(input))
			required = c("locus", "chromosome", "start", "stop"); present = required %in% names(input)
			if (!all(present)) input.error("missing %column% ", items.and=required[!present])
			super$initialize(input$locus, source, type)
			rownames(input) = NULL
			private$ranges = input[required[-1]]
		},

		get.contents = function() {return(cbind(locus=private$locus.names, private$ranges))}
	)
)


## annotation specified in lists of SNPs
## input can be provided as either a data.frame with LOC and SNPS column, with SNPS column entries semi-colon separated strings
## or as a list with a LOC and SNPS element of equal length, with SNPS a list of SNP ID vectors
ListAnnotation = R6::R6Class("ListAnnotation",
	inherit = LocusAnnotation,
	private = list(
		snp.ids = NULL,

		create.locus = function(index) {return(ListLocus$new(index, private$locus.names[index], private$snp.ids[[index]], type=private$type, source=private$source))},

		merge.internal = function(other, merge.mode) {
			data1 = self$get.contents(); data2 = other$get.contents()

			if (merge.mode == "union") locus.names = unique(c(data1$locus, data2$locus))
			else locus.names = intersect(data1$locus, data2$locus)

			snps = list()
			for (l in locus.names) {
				if (merge.mode == "union") curr = unique(c(data1$snp.ids[[l]], data2$snp.ids[[l]]))
				else curr = intersect(data1$snp.ids[[l]], data2$snp.ids[[l]])
				if (length(curr) > 0) snps[[l]] = curr
			}
			return(ListAnnotation$new(snps, source=paste0("merged (", merge.mode, ")"), type=private$type))
		}
	),
	public = list(
		initialize = function(input, sep=";", source=NULL, type=NULL) {
			if (is.list(input) && !is.null(names(input))) {
				found = c("locus", "snps") %in% tolower(names(input))
				if (all(found)) names(input) = tolower(names(input))  ## data.frame or list with 'locus' and 'snps' elements
				else if (!any(found) && !is.data.frame(input)) input = list(locus = names(input), snps=input)  ## named list of SNP ID vectors
				else input = NULL
			}
			if (!is.list(input)) input.error("invalid input for ListAnnotation object")

			super$initialize(input$locus, source, type)
			if (!is.data.frame(input)) {
				if (length(input$snps) != length(input$locus)) input.error("'locus' and 'snps' elements are not of same length")
				private$snp.ids = input$snps
			} else private$snp.ids = lapply(strsplit(input$snps, sep), trimws)
		},

		get.contents = function() {return(list(locus=private$locus.names, snp.ids=private$snp.ids))}
	)
)



## type can be a vector
LocusDefinition = R6::R6Class("LocusDefinition",
	private = list(
		id = NULL, name = NULL, type = NULL, source = NULL,

		match.index = function(snp.info) {undefined.error("match.index", class(self)[1])}
	),
	public = list(
		initialize = function(id, name, type, source) {private$id = id; private$name = name; private$type = if (!("unknown" %in% NULL)) type; private$source = source},

		get.id = function() {return(private$id)},
		get.name = function() {return(private$name)},
		get.type = function() {return(if (!is.null(private$type)) private$type else "unknown")},
		has.type = function(type) {
			if (is.null(type) || "unknown" %in% type) return(is.null(private$type))
			else return(all(type %in% private$type))
		},

		## snp.info: data.frame with internal.id, snp.id (for ListLocus), and chromosome/position (for RangeLocus) columns
		## outputs SNPindex object for variants in the locus
		snp.index = function(snp.info) {return(SNPindex$new(snp.info$internal.id[private$match.index(snp.info)]))}
	)
)


RangeLocus = R6::R6Class("RangeLocus",
	inherit = LocusDefinition,
	private = list(
		chromosome = NULL, start = NULL, stop = NULL,

		match.index = function(snp.info) {return(which(snp.info$chromosome == private$chromosome & snp.info$position >= private$start & snp.info$position <= private$stop))}
	),
	public = list(
		initialize = function(id, name, ..., type=NULL, source=NULL) {
			args = flatten.arglist(...)
			super$initialize(id, name, type, source)
			private$chromosome = args$chromosome; private$start = args$start; private$stop = args$stop
		},

		print = function(...) {
			printer = ObjectPrinter$new("Locus (genomic coordinates)")$parameter.list(locus=private$name, source=private$source, type=private$type, chromosome=private$chromosome, start=private$start, stop=private$stop)
			cat(printer$to.string())
		},

		get.coordinates = function() {return(list(chromosome=private$chromosome, start=private$start, stop=private$stop))}
	)
)


ListLocus = R6::R6Class("ListLocus",
	inherit = LocusDefinition,
	private = list(
		snp.ids = NULL,

		match.index = function(snp.info) {return(which(snp.info$snp.id %in% private$snp.ids))}
	),
	public = list(
		initialize = function(id, name, snp.ids, type=NULL, source=NULL) {
			super$initialize(id, name, type, source)
			private$snp.ids = tolower(snp.ids)
		},

		print = function(...) {
			printer = ObjectPrinter$new("Locus (SNP list)")
			printer$parameter.list(locus=private$name, source=private$source, type=private$type)
			snp.ids = head(private$snp.ids, 5)
			if (length(private$snp.ids) > length(snp.ids)) snp.ids = c(snp.ids, "...")
			printer$add.line("SNPs (", length(private$snp.ids), " total): ", paste(snp.ids, collapse=", "), add.space=F)
			cat(printer$to.string())
		},

		get.snps = function() {return(private$snp.ids)}
	)
)













