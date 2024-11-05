
min.snps = 100
min.pcs = 20
univ.thresh = 1e-5  

chromosome = 22


ref.prefix.base = "g1000" # prefix of reference data

info.file.name = ...  
# info file must contain 'phenotype' and 'filename' columns
# may contain 'cases', 'controls', 'unit', 'properties'

annot.file.name = ... # file with definitions to be analysed (columns: LOC, CHR, START, STOP)



add.root = function(file, root=".") {paste0(root, "/", file)}
pkg.dir = add.root("scripts")             ## directory containing the 'R' folder and DESCRIPTION / NAMESPACE files
data.dir = add.root("data")
out.dir = add.root("output")


ref.prefix = add.root(ref.prefix.base, root=data.dir)
info.file = add.root(info.file.name, root=data.dir)
annot.file = add.root(annot.file.name, root=data.dir)

pkgload::load_all(pkg.dir)


#########################

## load in and validate data files for specified chromosome, create DataSet object; input.dir is appended as prefix to file names specified in info file
data = process.input(info.file, sample.overlap.file=NULL, ref.prefix=ref.prefix, input.dir=data.dir, chromosomes=chromosome)

## create DataInterface object for all phenotypes included in 'data'
input = data$get.interface(settings=list(minimum.snps=min.snps, minimum.components=min.pcs))

# read in loci and set as annotation on input object
annot = read.loci(annot.file)
input$set.annotation(annot)


univariate = NULL; bivariate = NULL
while (input$next.locus()) { ## advances to next locus in internal annotation, returns FALSE and ends loop if no more remaining
  ## equivalent of process.locus(), any data failures are converted to warnings, with locus.data becoming NULL (errors will still halt the script)
  locus.data = catchErrors(input$process()$filter()$filter.pval(univ.thresh), all.failures=T)
  
  if (!is.null(locus.data)) {
    univariate = fill.rowmerge(list(univariate, locus.data$summarize()))
    if (locus.data$no.pheno() > 1) {
      res = catchErrors(compute.rg(locus.data), all.failures=T)
      bivariate = fill.rowmerge(list(bivariate, res))
    }
  }
}

write.table(univariate, file=paste0(out.dir, "/univariate_chr", chromosome, ".res"), row.names=F, quote=F, sep="\t")
write.table(bivariate, file=paste0(out.dir, "/bivariate_chr", chromosome, ".res"), row.names=F, quote=F, sep="\t")



