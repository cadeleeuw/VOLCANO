## merge list of data frames, aligning column order and filling in missing columns
fill.rowmerge = function(data) {
	data = data[!sapply(data, is.null)]
	columns = unique(unlist(lapply(data, names)))
	if (!all(sapply(data, ncol) == length(columns)) || !all(sapply(data, function(d) {all(names(d) == columns)}))) {
		for (i in seq_along(data)) {
			missing = !(columns %in% names(data[[i]]))
			if (any(missing)) data[[i]][columns[missing]] = NA
			data[[i]] = data[[i]][columns]
		}
	}
	return(do.call(rbind, data))
}


## subset matrix by row/colnames; cols/rows values of NULL return all rows/columns
## maintains row/colnames in output and correct dimensions even if size of one or both dimensions reduced to zero
matrix.subset = function(M, cols=NULL, rows=NULL) {
	if (!is.null(rows)) {
		if (is.null(rownames(M))) input.error("in matrix.subset(), cannot select rows if rownames is NULL")
		if (!all(rows %in% rownames(M))) input.error("in matrix.subset(), unknown rows selected")
		row.index = match(rows, rownames(M))
	} else row.index = 1:nrow(M)
	if (!is.null(cols)) {
		if (is.null(colnames(M))) input.error("in matrix.subset(), cannot select cols if colnames is NULL")
		if (!all(cols %in% colnames(M))) input.error("in matrix.subset(), unknown cols selected")
		col.index = match(cols, colnames(M))
	} else col.index = 1:ncol(M)

	return(matrix(M[row.index,col.index], nrow=length(row.index), ncol=length(col.index), dimnames=list(rownames(M)[row.index], colnames(M)[col.index])))
}

## add extra row/column to matrix M, with row/colname as specified
matrix.expand = function(M, name, values) {
	N = length(values)
	if (ncol(M) == N && !is.null(colnames(M))) { ## add as column
		M = cbind(M, values); colnames(M)[ncol(M)] = name
	} else if (nrow(M) == N && !is.null(rownames(M))) { ## add as row
		M = rbind(M, values); rownames(M)[nrow(M)] = name
	} else if (ncol(M) == N-1 && ncol(M) == nrow(M) && !is.null(colnames(M)) && !is.null(rownames(M))) { ## add as row and column
		M = rbind(cbind(M, values[-N]), values)
		colnames(M)[N] = rownames(M)[N] = name
	} else input.error("incompatible dimensions in matrix.expand()")
	return(M)
}


## create data.frame with all combinations of named input vectors or input named list
combinations = function(...) {
	args = flatten.arglist(...); values = NULL
	for (param in names(args)) {
		if (!is.null(values)) values = cbind(values[rep(1:nrow(values), each=length(args[[param]])),], args[[param]])
		else values = matrix(args[[param]], ncol=1)
	}
	values = as.data.frame(values); names(values) = names(args)
	return(values)
}


## create diagonal matrix, optionally set row/column names
## argument 'diagonal' is always treated as a vector, even if length is one
diag.matrix = function(diagonal, add.names=NULL) {
	M = diag(diagonal, nrow=length(diagonal))
	if (!is.null(add.names)) {
		if (length(add.names) != length(diagonal)) input.error("add.names argument for diag.matrix() does not match size of specified matrix")
		colnames(M) = rownames(M) = add.names
	}
	return(M)
}


## mirror square matrix values and names across the diagonal, and optionally set new row/column names
make.symmetric = function(M, add.names=NULL, use=c("lower", "upper")) {
	if (is.null(dim(M)) || nrow(M) != ncol(M)) input.error("matrix input for make.symmetric() is not square")
	if (match.arg(use) == "upper") M = t(M)

	if (!is.null(add.names)) {
		if (length(add.names) != nrow(M)) input.error("add.names argument for make.symmetric() does not match size of matrix input")
		rownames(M) = add.names
	}

	colnames(M) = rownames(M)
	M[upper.tri(M)] = t(M)[upper.tri(M)]
	return(M)
}


## invert using eigendecomposition, checking invertibility
## NB: symmetry of M is assumed
eigen.invert = function(M, threshold=globals$decomposition.eigen.threshold) {
	decomp = eigen(M); no.vars = ifelse(!is.null(dim(M)), nrow(M), 1)
	if (any(decomp$values / sum(decomp$values) < threshold / no.vars)) data.error("matrix is not invertible")
	if (no.vars > 1) return(decomp$vectors %*% diag(1/decomp$values) %*% t(decomp$vectors))
	else return(matrix(M, 1, 1))
}

## sum columns of input matrix (considerably faster than apply(data, 1, sum) for large inputs)
column.sum = function(data) {return((data %*% rep(1, ncol(data)))[,1])}

## apply logical function to boolean matrix (considerably faster than apply(bool.data, 1, op) for large inputs)
column.none = function(bool.data) {return(column.sum(bool.data) == 0)}
column.any = function(bool.data) {return(column.sum(bool.data) > 0)}
column.all = function(bool.data) {return(column.sum(bool.data) == ncol(bool.data))}



## paste columns of input data.frame together
paste.columns = function(data, columns, sep=" ") {
	if (length(columns) > 1) {
		if (is.numeric(columns)) columns = names(data)[columns]
		expr.str = paste0("paste(", paste0("data[[\"", columns, "\"]]", collapse=", "), ", sep=\"", sep, "\")")
		return(eval(parse(text=expr.str)))
	} else return(data[[columns]])
}



