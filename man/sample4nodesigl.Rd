% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample4nodesigl.R
\name{sample4nodesigl}
\alias{sample4nodesigl}
\title{Preparing the input 'sample file' for Phylocom}
\usage{
sample4nodesigl(fdisease, tree, disout)
}
\arguments{
\item{fdisease}{A csv file of medicinal plants with the taxa assumed to have
the correct names. It is assumed that the csv file is formatted like the
"diseasename_cleanedup_exmatch.csv" file returned by the 'clean_and_match'
function, possibly updated by  user with the taxa included in the
"diseasename_cleanedup_qumatch_threshold.csv" file. See the
'clean_and_match' function for details and options.}

\item{tree}{A user-provided tree in Newick format whose tips are assumed to be
correctly labelled. See 'clean_and_match' function for details and renaming
options. Output files of 'tree4nodesigl' function can be used (assuming that
tips are correctly labelled).}

\item{disout}{Character string specifying the name of the disease of interest
to be added as a suffix to the end of the name of the generated "sample"
file. It can be the same name as in the "fdisease" files or any other
convenient name chosen by the user.}
}
\value{
A tab delimited file is saved in the current working directory.
According to the Phylocom manual, one row per taxon with three columns
("Sample name", "Abundance" and "Species code") is used. By default "Sample
name" and "Abundance" columns are set to "disease" and "1", respectively,
for all taxa. Column "Species code" includes only medicinal plant names that
give an exact match to the taxa present in the input tree.
}
\description{
This function prepares the "sample" file to be used with the
software Phylocom (Webb et al. 2008), starting from a tree file and a
medicinal plants file provided by the user.
}
\details{
This function prepares the "sample" file to be used in "nodesigl"
Phylocom analysis and it is intended for use with the output files of the
'clean_and_match' and 'tree4nodesigl' functions. Custom files, such as
subsamples of medicinal plants or any subtree of interest, can also be used
as long as they are formatted correctly (please double check column
headers!) and the included taxa labeled correctly (see Parameters for
details). For this reason the function performs an additional check and only
taxa included in both input files are retained in the final "sample" file.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
 # Use tempdir() to find the current temp directory.

 library(ape)

 # Generate a random tree with 10 tips and save it in the current temp folder.
 tr<-rtree(10)
 tree<-write.tree(tr)
 treefile <- tempfile("tree",fileext = ".tree")
 cat(tree, file = treefile, sep = "\n")

 #Create a medicinal plants list including also taxa not present in the tree
 # and save it in the current temp folder.
 dis<-data.frame(Plant.Name=paste0("t",c(2:7,12:15)))
 disfile<-tempfile("disease", fileext = ".csv")
 write.csv(dis, file= disfile, quote =FALSE,row.names =FALSE)

 # Get path to the current temp folder and set it as the current working directory.
 path<-paste0(normalizePath(tempdir(), winslash = "/"),"/")
 setwd(path)

 #Test.
 sample4nodesigl(fdisease=disfile,tree=treefile,fileout="disease")

 #Tidy up. Remove the output files from the current temporary directory.
 toremove<-list.files(path, pattern = "disease|.tree|sample")
 file.remove(paste0(path,toremove))
\dontshow{\}) # examplesIf}
}
\references{
\itemize{
\item Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics
and evolutionary analyses in R.” Bioinformatics, 35, 526-528.
doi:10.1093/bioinformatics/bty633.\cr
\item Webb, C. O.; Ackerly, D. D. & Kembel, S. W. (2008) Phylocom: software for the
analysis of phylogenetic community structure and trait evolution.
Bioinformatics 24: 2098-2100.
\item Zecca, G., Toini, E., Labra, M, Grassi, F. (2025) Accelerating the
identification and the prioritisation of new plants with medicinal
potential: the  pm4mp R package.\cr
}
}
\seealso{
\code{\link[=clean_and_match]{clean_and_match()}},\code{\link[=tree4nodesigl]{tree4nodesigl()}}
}
\author{
Giovanni Zecca, \email{giovanni.zecca@unimib.it}
}
