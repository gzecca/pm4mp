#'
#'@title Phylocom's 'nodesigl' analysis from R
#'
#'@description This function is a simple wrapper for running Phylocom's
#'  'nodesigl' analysis (Webb et al., 2008) directly from R.
#'
#'@details This function calls the software Phylocom and performs multiple runs
#'  of the 'nodesigl' analysis directly from R. It assumes that input files (see
#'  Parameters) and a copy of the Phylocom binaries are in the current working
#'  directory. Please refer to the Phylocom manual for OS details.
#'
#'@param start,stop Two integers, corresponding to the first and last
#'  independent runs of the 'nodesigl' analysis to be performed. The number of
#'  independent runs is equal to: length(start:stop). For example, setting
#'  start=100 and stop=300 will run two hundred independent analysis replicates.
#'
#'@param s Character string, indicating the name of the sample file. Typically,
#'  the output of 'sample4nodesigl' function.
#'
#'@param f Character string, indicating the name of the phylogeny file.
#'  The output of 'tree4nodesigl' function NOT including branch lengths must be used.
#'  Empty files are obtained if branch lengths are included!
#'
#'@param r Integer, indicating the number of randomizations to use in each run
#'  of the analysis (default = 999).
#'
#'@return N text files, named "Nodesigl_start"..."Nodesigl_stop", are saved in
#'  the current working directory, N being the number of independent runs (i.e.,
#'  N = length(start:stop)). Each file contain the output of a single 'nodesigl'
#'  analysis based on 'r' randomizations. Please see the Phylocom manual for a
#'  description of the contents of the files. Since each run is independent, the
#'  output of different call of the'nodesiglR' function can be easily combined
#'  by specifying different values for the 'start' and 'stop' parameters.
#'
#'@author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#'
#'@references
#' - Webb, C. O.; Ackerly, D. D. & Kembel, S. W. (2008) Phylocom: software for
#'   the analysis of phylogenetic community structure and trait evolution.
#'   Bioinformatics 24: 2098-2100.
#' - Zecca, G., Toini, E., Labra, M, Grassi, F. (2025) Accelerating the
#'   identification and the prioritisation of new plants with medicinal
#'   potential: the  pm4mp R package.\cr
#'
#'@seealso [tree4nodesigl()], [sample4nodesigl()]
#'
#'@examplesIf interactive()
#'# WARNING: to run the example a copy of the Phylocom binaries MUST be in the
#'# current temp directory. DO NOT RUN otherwise!
#'# Use tempdir() to find the current temp directory.
#'
#'  library(ape)
#'  # Generate a random tree with 10 tips and save it in the current temp folder.
#'  tr<-rtree(10)
#'  tree<-write.tree(tr)
#'  treefile <-tempfile("tree", fileext = ".tree")
#'  cat(tree, file = treefile, sep = "\n")
#'
#'  #Create a sample file
#'  sample<-c("disease\t1\tt1
#'  disease\t1\tt3
#'  disease\t1\tt5
#'  disease\t1\tt4")
#'  samplefile<- tempfile("sample", fileext ="")
#'  cat(sample, file = samplefile, sep = "\n")
#'
#'  #Get the temporary directory path and set it as current working directory.
#'  path<-paste0(normalizePath(tempdir(), winslash = "/"),"/")
#'  setwd(path)
#'
#'  #Test.
#'  tree4nodesigl(treefile, "random_tree")
#'  # This is a toy example. Many more repetitions are normally needed
#'  nodesiglR(start=1,stop=20,s=samplefile, f="phylo_nobl_random_tree")
#'
#'  #Tidy up. Remove the output files from the current temporary directory.
#'  toremove<-list.files(path, pattern = "tree|phylo_|Nodesigl_|sample")
#'  file.remove(paste0(path,toremove))
#'
#'@export
#'

nodesiglR <-function(start, stop, s, f, r= 999){
  for (i in start:stop){
    system2(command="phylocom", args= c("nodesigl", paste0("-r ", r),"-m 2", paste0("-s ",s), paste0("-f ",f)), stdout=paste0("Nodesigl_",as.character(i),".txt"))}
}
