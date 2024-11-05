

InputInterface = R6::R6Class("InputInterface",
	private = list(
		parameter.index = NULL, input.names = NULL, input.type = NULL,
		parameter.map = NULL, parameter.names = NULL, match.counts = NULL,

		create.map = function() {
			duplicates = duplicated(unlist(private$parameter.index))
			if (any(duplicates)) input.error("invalid index, %", private$input.type, "% ", items.and=unique(unlist(private$parameter.index)[duplicates]), " %is% %[each]% mapped to multiple parameters")

			private$parameter.map = list()
			for (h in names(private$parameter.index)) private$parameter.map[[h]] = unique(private$input.names[tolower(private$input.names) %in% private$parameter.index[[h]]])

			private$parameter.names = names(private$parameter.map)
			private$match.counts = unlist(lapply(private$parameter.map, length))
		}
	),
	public = list(
		initialize = function(input.names, parameter.index, input.type="variable", include.param.names=T) {
			private$input.names = input.names; private$input.type = input.type

			if (include.param.names) {for (h in names(parameter.index)) parameter.index[[h]] = c(parameter.index[[h]], h)}
			private$parameter.index = lapply(parameter.index, function(cols) {unique(tolower(cols))})

			private$create.map()
		},

		## add additional parameter mappings, adding new parameters and extending mapping of existing parameters (if append.existing=T)
		## returns a list of parameters included in the input mapping that already existed
		update.parameters = function(add.index, append.existing=T) {
			existing = c()
			for (param in names(add.index)) {
				if (param %in% names(private$parameter.index)) {
					existing = c(existing, param)
					if (!append.existing) next
				}
				private$parameter.index[[param]] = unique(c(private$parameter.index[[param]], tolower(add.index[[param]])))
			}
			private$create.map()
			return(existing)
		},

		get.variables = function() {return(private$input.names)},

		get.available = function() {return(private$parameter.names[private$match.counts > 0])},
		get.duplicates = function() {return(private$parameter.names[private$match.counts > 1])},
		get.missing = function() {return(private$parameter.names[private$match.counts == 0])},

		## return map to first/last mapped input element, trimming away unmapped parameters unless trim.empty=F
		## if parameters=NULL, return full map of all parameters, otherwise subset to selected parameters
		## if drop.missing=T, unknown parameters are discarded, otherwise they are included with an NA value
		get.map = function(parameters=NULL, use=c("first", "last"), drop.missing=T, trim.empty=T) {
			select = if (match.arg(use) == "first") head else tail
			map = lapply(private$parameter.map, select, 1)
			if (trim.empty) map = map[private$match.counts > 0]
			if (!is.null(parameters)) {
				map = map[names(map) %in% parameters]
				if (!drop.missing) map[parameters[!(parameters %in% names(map))]] = NA
				map = map[order(match(names(map), parameters))]
			}
			return(map)
		},

		has.parameters = function(parameters, aggregate=F) {status = parameters %in% private$parameter.names; return(if (aggregate) all(status) else status)},
		has.available = function(parameters, aggregate=F) {status = parameters %in% self$get.available(); return(if (aggregate) all(status) else status)}
	)
)


FileInputInterface = R6::R6Class("FileInputInterface",
	inherit = InputInterface,
	private = list(
		filename = NULL, file.header = NULL, global.parameters = list(),

		read.header = function() {undefined.error("read.header", class(self)[1])}
	),
	public = list(
		initialize = function(filename, parameter.index, include.param.names=T) {
			private$filename = filename
			private$read.header()
			super$initialize(private$file.header, parameter.index, input.type="column",	include.param.names=include.param.names)
		},

		get.filename = function() {return(private$filename)},
		get.globals = function() {return(private$global.parameters)}
	)
)


TextFileInputInterface = R6::R6Class("TextFileInputInterface",
	inherit = FileInputInterface,
	private = list(
		comment.lines = 0,

		read.header = function() {
			comments = NULL; curr.step = 1; increment = 100
			while (length(private$file.header) == 0) {
				current = scan(private$filename, sep="\n", what="", quiet=T, strip.white=T, blank.lines.skip=F, skip=(curr.step-1)*increment, nlines=increment)
				if (length(current) == 0) break

				input.lines = grep("^[^#]", current) ## match a non-empty line not starting with a # (NB: leading/trailing whitespaces are stripped by scan, so empty lines are "")
				if (length(input.lines) > 0) {
					if (input.lines[1] > 1) {
						comments = c(comments, current[1:(input.lines[1]-1)])
						private$comment.lines = length(comments)
					}
					private$file.header = strsplit(current[input.lines[1]], "\\s+")[[1]]
					break
				} else {
					comments = c(comments, input.lines)
					curr.step = curr.step + 1
				}
			}
			if (length(private$file.header) == 0) input.error("input file ", private$filename, " contains no header")

			comments = trimws(gsub("^#", "", comments))
			parameters = grep("^\\S+\\s*:\\s*\\S+", comments[comments != ""], value=T)   ## select parameter comments (non-whitespace string, (optional whitespace), colon, (optional whitespace), any amount of non-whitespace

			if (length(parameters) > 0) {
				param.parts = lapply(strsplit(parameters, ":"), function(v) {trimws(c(v[1], paste(v[-1], collapse=":")))})  ## split on first c
				for (i in 1:length(param.parts)) private$global.parameters[[param.parts[[i]][1]]] = preprocess.specifier(param.parts[[i]][2])
			}
		}
	),
	public = list(
		load.input = function(input.type=NULL, error.tag=NULL, defer.loading=F) {return(TextFileInput$new(self, input.type=input.type, error.tag=error.tag, defer.loading=defer.loading))},
		load.data = function() {
			data = data.table::fread(private$filename, data.table=F, showProgress=F, skip=private$comment.lines)
			if (length(private$file.header) != ncol(data)) input.error("invalid header for file ", private$filename)
			return(data)
		}
	)
)


InputData = R6::R6Class("InputData",
	private = list(
		input.type = NULL, error.tag = NULL,
		data = NULL, param.interface = NULL,

		error = function(...) {input.error(..., .label=private$error.tag)},

		match.parameters = function(..., target, as.bool=F, check.unknown=T) {parameters = c(...)
			unknown = !private$param.interface$has.parameters(parameters)
			if (check.unknown && any(unknown)) private$error("%parameter% ", items.and=parameters[unknown], " %is% unknown")
			if (as.bool)	return(parameters %in% target)
			else return(parameters[parameters %in% target])
		},

		check.required = function(...) {
			missing = private$match.parameters(..., target=private$param.interface$get.missing())
			if (length(missing) > 0) private$error("no input %column% %was% found for required %parameter% ", items.and=missing)
		},

		check.duplicates = function(...) {
			duplicates = private$match.parameters(..., target=private$param.interface$get.duplicates())
			if (length(duplicates) > 0) throw.warning("matched multiple columns for %parameter% ", items.and=duplicates, ", using first matched column %[for each]%", .label=private$error.tag)
		},

		extract.parameters = function(...) {
			private$check.required(...)
			private$check.duplicates(...)

			parameters = c(...)
			index = private$param.interface$get.map(parameters, use="first", drop.missing=T)

			out = data.frame(private$data[,unlist(index)], stringsAsFactors=F)
			names(out) = names(index)
			return(out)
		}
	),
	public = list(
		initialize = function(data, interface, input.type=NULL, error.tag=NULL) {
			private$data = data; private$param.interface = interface
			private$error.tag = error.tag; private$input.type = input.type

			duplicate = duplicated(tolower(interface$get.variables()))
			if (any(duplicate)) private$error("duplicate column %name% ", items.and=interface$get.variables()[duplicate])
		},

		has.parameters = function(..., aggregate=T) {
			available = private$match.parameters(..., target=private$param.interface$get.available(), as.bool=T, check.unknown=F)
			return(if (aggregate) all(available) else available)
		},

		get.data = function() {return(private$extract.parameters(private$param.interface$get.available()))},
		get.subset = function(..., optional=NULL) {
			parameters = c(...)
			if (!is.null(optional)) parameters = c(parameters, optional[private$param.interface$has.available(optional)])
			return(private$extract.parameters(parameters))
		},
		get.parameter = function(param.name) {return(private$extract.parameters(param.name)[,1])}
	)
)


DataframeInput = R6::R6Class("DataFrameInput",
	inherit = InputData,
	public = list(
		initialize = function(data, header.index, input.type=NULL, error.tag=NULL) {
			if (is.null(error.tag)) {
				if (is.null(input.type)) error.tag = "in data.frame input"
				else error.tag = paste0("in ", input.type, " data.frame input")
			}
			if (!is.data.frame(data)) fatal.error("input to DataframeInput object is not a data frame", .label=private$error.tag)

			interface = InputInterface$new(names(data), header.index, input.type="column", include.param.names=T)
			super$initialize(data, interface, input.type=input.type, error.tag=error.tag)
		}
	)
)



TextFileInput = R6::R6Class("TextFileInput",
	inherit = InputData,
	private = list(
		filename = NULL,

		extract.parameters = function(...) {
			if (is.null(private$data)) self$load.data()
			return(super$extract.parameters(...))
		}
	),
	public = list(
		initialize = function(input, header.index, input.type=NULL, error.tag=NULL, defer.loading=F) {
			if (has.type(input, "TextFileInputInterface")) interface = input
			else interface = TextFileInputInterface$new(input, header.index, include.param.names=T)
			private$filename = interface$get.filename()

			if (is.null(error.tag)) {
				if (is.null(input.type)) error.tag = paste0("in file ", private$filename)
				else error.tag = paste0("in ", input.type, " file ", private$filename)
			}

			super$initialize(NULL, interface, input.type=input.type, error.tag=error.tag)
			if (!defer.loading) self$load.data()
		},

		get.globals = function() {return(private$param.interface$get.globals())},

		load.data = function() {if (is.null(private$data)) private$data = private$param.interface$load.data()}
	)
)







