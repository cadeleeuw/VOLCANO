has.type = function(object, types, R6=T) {return(all(c(types, if (R6) "R6") %in% class(object)))}


## input is named arguments (or named list) of format [variable]=[types] (can be multiple types), the variable is checked to see if it has all those types
## variables are pulled from the environment where check.types() is called; if called (directly) inside a function this will be the internal environment of the function, including the input arguments of that function
check.types = function(..., .func=NULL, .class=NULL) {
	env = parent.frame(); args = list(...); label = NULL
	if (!is.null(.func)) label = paste0("in function ", if (!is.null(.class)) paste0(.class[1], "$", collapse=""), .func)
	else if (!is.null(.class)) label = paste0("initializing ", .class[1], " object")

	for (obj in names(args)) {
		if (!(obj %in% names(env))) fatal.error("argument '", obj, "' does not exist in input", .label=label)
		if (!has.type(env[[obj]], args[[obj]], R6=F)) fatal.error("argument '", obj, "' must have %type% ", items.and=args[[obj]], .label=label)
	}
}


## validates named list of settings 'input' against a reference of possible settings, completing partial names if possible and checking for duplicates
## if merge=T, values in input are copied to reference and reference is returned; otherwise, input is returned
validate.settings = function(input, reference, source, merge=F) {
	if (!is.character(source)) source = paste(class(source)[1], "object")
	if (length(input) == 0) input = list()
	else if (!is.list(input) || is.null(names(input))) input.error("invalid settings specification in ", source)
	else input = input[!sapply(input, is.null)]

	names(input) = complete.names(names(input), names(reference))
	unknown = !(names(input) %in% names(reference))
	if (any(unknown)) input.error("unknown %parameter% ", items.and=names(input)[unknown], " in ", source)

	multiple = sapply(input, length) > 1
	if (any(multiple)) input.error("multiple values provided for %parameter% ", items.and=names(input)[unknown], " in ", source)

	if (merge) {
		for (param in names(input)) reference[[param]] = input[[param]]
		return(reference)
	} else return(input)
}


## returns named list of input values: unnamed list inputs are recursively flattened to extract named elements
## subsequently, NULL and unnamed elements are removed (unless keep.unnamed=T)
flatten.arglist = function(..., filter.NA=F, keep.unnamed=F) {
	input = list(...); output = list()
	while (length(input) > 0) {
		if (length(input[[1]]) > 0) {
			is.named = !is.null(names(input)[1]) && names(input)[1] != ""
			if (!is.list(input[[1]]) || is.named) {
				if (is.named || keep.unnamed) {
					output[[length(output)+1]] = input[[1]]
					names(output)[length(output)] = names(input)[1]
				}
				input = input[-1]
			} else input = c(input[[1]], input[-1])
		} else input = input[-1]
	}
	if (length(output) > 0) {
		drop = unlist(lapply(output, function(o) {is.null(o) || (filter.NA && length(o) == 1 && is.na(o))}))
		output = output[!drop]
	}

	if (!is.null(names(output))) names(output)[is.na(names(output))] = ""
	return(output)
}


## completes elements in 'targets' to match element in 'range' (ie. exact match, or matches start of single element in 'range')
## if no match, returns original value, or NA if failed.NA=T
complete.names = function(targets, range, failed.NA=F) {
	exact = targets %in% range
	for (i in which(!exact)) {
		partial = targets[i] == substr(range, 1, nchar(targets[i]))
		if (sum(partial) == 1) targets[i] = range[partial]
		else if (failed.NA) targets[i] = NA
	}
	return(targets)
}


## performs partial match of targets to range, returns list with valid matches, as well as ambiguous and unknown targets
## 'matches' output value is in same order as 'targets' input, and same length if none failed and all entries in 'targets' are unique and no wildcards were expanded
## wildcard * matches any number of characters, wildcard ? matches a single character
partial.match = function(targets, range, expand.wildcard=F) {
	matches = list(); is.wild = c()
	for (t in targets) {
		if (!any(range == t)) {
			curr = unique(grep(t, range, value=T, fixed=T))
			matches[[t]] = curr[substr(curr, 1, nchar(t)) == t]
		} else matches[[t]] = t

		if (expand.wildcard && length(matches[[t]]) == 0 && grepl("\\*|\\?", t)) {
			pattern = gsub("\\?", ".", gsub("\\*", ".*", gsub("\\.", "\\\\.", t))) ## escape periods, convert wildcards to regular expression
			pattern = paste0("^", pattern, "$")
			res = try(grep(pattern, range, value=T), silent=T)
			if (class(res) != "try-error")	matches[[t]] = res
			if (length(matches[[t]]) > 0) is.wild = c(is.wild, t)
		}
	}

	count = sapply(matches, length)
	if (length(is.wild) > 0) count[count > 1 & names(matches) %in% is.wild] = 1
	return(list(
		matches = as.vector(unlist(matches[count == 1])),
		ambiguous = targets[count > 1],
		unknown = targets[count == 0],
		failure = any(count != 1)
	))
}


## scrub whitespace around commas and equality signs, and split by whitespace / semi-colons, removing empty elements
## also removes leading/trailing commas and commas directly before/after equality signs
preprocess.specifier = function(input.str) {
	input.str = gsub(",+", ",", gsub("\\s*,\\s*", ",", gsub("\\s*=\\s*", "=", trimws(input.str))))   ## remove whitespace around commas and equality signs, and collapse multiple commas
	elems = strsplit(trimws(input.str), "(\\s*;\\s*)|(\\s+)")[[1]]   ## split by whitespace / semicolons
	elems = gsub(",*=,*", "=", gsub("(^,)|(,$)", "", elems))   ## remove extraneous commas
	return(elems[elems != ""])
}

## parse input.string of the format "name1=value1,value2,... name2=value3,value4,...; name3=value5 value6" etc.
## whitespace within primary name/value elements around commas and equality signs is ignored
## primary are separated by whitespace and/or semi-colons
## if infer.names=T, unnamed elements like 'value6' are treated as 'value6=value6'
parse.value.string = function(input.str, allow.unnamed=T, infer.names=allow.unnamed) {
	raw = preprocess.specifier(input.str)
	pairs = strsplit(raw, "="); counts = sapply(pairs, length)

	invalid = counts > 2 | grepl("(^=)|(=$)", raw)  ## more than one '=', or '=' at beginning or end of element
	if (any(invalid)) input.error("invalid specifier %token% ", items.and=raw[invalid])

	output = strsplit(sapply(pairs, tail, 1), ",")
	names(output) = sapply(pairs, function(spec) {ifelse(length(spec) == 2, spec[1], "")})

	if (infer.names) {
		infer = names(output) == "" & sapply(output, length) == 1
		if (any(infer)) names(output)[infer] = unlist(output[infer])
	}

	if (!allow.unnamed && any(names(output) == "")) input.error("specification string '", input.str, "' contains unnamed ", if (infer.names) "multi-value ", "elements")
	return(output)
}



## check list of (partial) phenotype names against list of available phenotypes
## expands * wildcard to match multiple characters if not already matched as literal *, and same for expanding ? wildcard to match single character
## throws error if any are ambiguous or unmatched, returns vector with corresponding full phenotype names otherwise (removing any duplicates)
## returns NULL if phenotypes argument is NULL and character(0) if phenotypes argument is empty
check.phenotypes = function(phenotypes, available, label=NULL, discard.ambiguous=F, discard.unknown=F) {
	if (length(phenotypes) > 0) {
		index = partial.match(phenotypes, available, expand.wildcard=T)

		if (!is.null(label)) label = paste0(" in ", label)
		if (!discard.ambiguous && length(index$ambiguous) > 0) input.error("ambiguous partial phenotype %name% ", items=index$ambiguous, " %matches% multiple phenotypes", label)
		if (!discard.unknown && length(index$unknown) > 0) input.error("%phenotype% ", items.and=index$unknown, " %is% not present", label)
		return(unique(index$matches))
	} else {
		if (is.null(phenotypes)) return(NULL)
		else return(character(0))
	}
}



## if resolve.chr=T and filename contains [CHR] token, considered to exist if at least one chromosome-specific file matching the filename template is found
check.files.exist = function(..., resolve.chr=F) {
	file.names = c(...)
	if (length(file.names) > 0) {
		exist = !is.na(file.names) & file.exists(file.names)
		if (resolve.chr) {
			tpl = !is.na(file.names) & grepl("[CHR]", file.names, fixed=T)
			if (any(tpl))	exist[tpl] = unlist(lapply(file.names[tpl], function(f) {nrow(chromosome.files(f)) > 0}))
		}
		if (!all(exist)) input.error("missing input %file% ", items.and=file.names[!exist], "; please ensure correct file names and paths have been provided")
	}
}

check.dirs.exist = function(...) {
	dirs = c(...)
	if (length(dirs) > 0) {
		valid = !is.na(dirs) & dir.exists(dirs)
		if (!all(valid)) input.error("missing or invalid %directory% ", items=dirs[!valid], "; please ensure correct paths have been provided")
	}
}

chromosome.files = function(file.tpl, chromosomes="all", prune.missing=T) {
	if (grepl("[CHR]", file.tpl, fixed=T)) {
		files = sapply(c(1:23, "X", "x"), function(c) {gsub("[CHR]", c, file.tpl, fixed=T)})
		files[!file.exists(files)] = NA
		chr.x = files[23:25]
		if (any(!is.na(chr.x))) files[23] = chr.x[!is.na(chr.x)][1]
		out = data.frame(chromosome=1:23, file=files[1:23], stringsAsFactors=F)

		chromosomes = validate.chromosomes(chromosomes)
		if (chromosomes[1] != "all") out = out[out$chromosome %in% chromosomes,]
		if (prune.missing) out = out[!is.na(out$file),]
		return(out)
	} else return(NULL)
}

## validate list of chromosome codes, and condense to 'all' if all chromosomes included
validate.chromosomes = function(chromosomes, condense=T) {
	if (!is.null(chromosomes) && !any(chromosomes == "all")) {
		chromosomes[tolower(chromosomes) == "x"] = 23
		unknown = !(chromosomes %in% 1:23)
		if (any(unknown)) input.error("unknown chromosome %code% ", items.and=chromosomes[unknown])
		chromosomes = sort(unique(as.numeric(chromosomes)))
	} else chromosomes = 1:23

	if (condense && length(chromosomes) == 23 && all(chromosomes == 1:23)) return("all")
	else return(chromosomes)
}


## check if covariance matrix is valid
## allow.NA=T setting only allows NA covariance values, not variances
## if check.names=T, will mirror row/column names if one is NULL
check.covariance = function(M, allow.NA=F, check.names=T, label="covariance matrix", tolerance=1e-10) {
	if (is.null(dim(M)) || nrow(M) != ncol(M)) input.error(label, " is not a square matrix")
	if (!allow.NA && any(is.na(M))) input.error(label, " contains missing values")
	if (!isSymmetric(M, tol=tolerance, check.attributes=F)) input.error(label, " is not symmetrical")
	if (any(is.na(diag(M)) | diag(M) <= 0)) input.error(label, " contains missing, negative or zero variance values")
	if (any(abs(cov2cor(M)) > 1, na.rm=T)) input.error(label, " contains correlations exceeding 1 or -1")

	if (check.names) {
		if (is.null(rownames(M))) rownames(M) = colnames(M)
		else if (is.null(colnames(M))) colnames(M) = rownames(M)
		if (!all(rownames(M) == colnames(M))) input.error("row and column names of ", label, " do not match")
	}
	return(M)
}








