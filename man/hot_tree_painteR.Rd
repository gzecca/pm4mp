% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hot_tree_painteR.R
\name{hot_tree_painteR}
\alias{hot_tree_painteR}
\title{Draw an "independent hot tree" with a graphical summary of multiple
'nodesigl' analysis}
\usage{
hot_tree_painteR(
  inputTrees,
  inputNodes,
  import.probs = FALSE,
  inputProbs,
  indtr,
  tftype = "i",
  ncex = 0.6,
  br = TRUE,
  br.ncex = 0.8,
  br.pch = 8,
  uom = TRUE,
  branch.unit = "Ma",
  autolabel = TRUE,
  ttype = NULL,
  tfsize = NULL,
  arrow = TRUE,
  par.arrow = NULL,
  par.legend = NULL,
  ...
)
}
\arguments{
\item{inputTrees}{The "IndependentHotTreesNewick" file produced by the
'nodesigl_harvesteR' function, which includes the list of identified
"independent hot trees" written in Newick format.}

\item{inputNodes}{The "IndependentHotNodesByTrees" csv file produced by the
'nodesigl_harvesteR' function, which includes information summarised from
multiple 'nodesigl' analysis replicates.}

\item{import.probs}{Logical, specifying whether relative probabilities
computed by the 'hmpp' function must be attached to taxa labels
and displayed in the final plot (deafault = FALSE).}

\item{inputProbs}{The csv file produced by the 'hmpp' function,
which contain the relative probability of each taxon of being a medicinal
plant.}

\item{indtr}{Integer, specifying the tree to be drawn among those included in
the "IndependentHotTreesNewick" file . Trees are numbered in increasing
order from the first to the last row of the input file.}

\item{tftype}{Font type (default="i"). See the 'plotSimmap' function in the
package 'phytools' for details.}

\item{ncex}{Numeric, specifying the expansion applied to the symbols drawn at
generic "hot nodes"(i.e, black filled circles). Default value is ncex= 0.6.}

\item{br}{Logical, indicating whether best ranking "hot nodes" should be
marked with a different symbol or not. Default is TRUE.}

\item{br.ncex}{Numeric, specifying the expansion applied to the symbols drawn
at best ranking "hot nodes". Default value is ncex= 0.8.}

\item{br.pch}{Numeric, which specifies the symbol to be drawn corresponding to
the best ranking "hot nodes" when the parameter  'br' is set to TRUE.
Default value is br.pch= 8. See the R graphical parameter 'pch' for
available options.}

\item{uom}{Logical, specifying whether the unit of measurement of branch
lengths must be added to the scale bar or not (default= TRUE).}

\item{branch.unit}{Character string, specifying the unit used to measure
branch lengths (default="Ma"). Ignored when uom=FALSE.}

\item{autolabel}{Logical, specifying whether the function should automatically
set the \emph{'ttype'} and the \emph{'tfsize'} parameters or not. If TRUE ,
the 'ttype' parameter is set to "fan" if Ntip>50, to "phylogram" otherwise.
A rough attempt to find a appropriated value for the 'tfsize' parameter is
done and the selected value is printed on the screen. It should be regarded
as a reasonable starting point for a finer tuning of this parameter.}

\item{ttype}{Tree type. See the 'plotSimmap' function in the package
'phytools' for details. (Defaul= "phylogram" or "fan", depending on the tree
size).}

\item{tfsize}{Numeric, specifying the relative font size for tip labels.}

\item{arrow}{Logical, specifying whether known medicinal plants should be
indicated by arrows in the final graph or not (default= TRUE).}

\item{par.arrow}{List, including graphical parameters to modify the aspect of
the arrows pointing to known medicinal taxa. Ignored when 'arrow'=FALSE. See
the the 'add.arrow' function in 'phytools' package for details and options.
Please note that \emph{'tree'}, \emph{'tip'} and \emph{'col'} parameters are
already defined in the 'hot_tree_painteR' function and must not be modified
here. When 'arrow'=FALSE and 'par.arrow'=NULL the appearance of the arrows
is automatically determined by R with variable results depending on the size
of the tree (generally better results are obtained with small trees).}

\item{par.legend}{List, including graphical parameters internally passed to
the 'add.simmap.legend' function of the 'phytools' package. Please refer to
the 'add.simmap.legend' help for details and options.}

\item{...}{Additional arguments internally passed to the 'plotSimmap'
function. Please note that \emph{'tree'}, \emph{'colors'},\emph{'pts'},
\emph{'ftype'},\emph{'fsize'} and \emph{'type'} parameters are already
defined by the 'hot_tree_painteR' function and must not be modified here.
See the 'plotSimmap' function in the package 'phytools' for details and
available options.}
}
\value{
A tree with colored branches and symbols drawn at the "hot nodes" is
printed on the graphic device. The user is asked to interactively click
where to draw the the legend of colours, the scale bar and the unit used to
measure branch lengths (if requested).
}
\description{
This function draws a user-selected "independent hot tree",
graphically displaying the output of multiple runs of the 'nodesigl'
analysis as summarised by the 'nodesigl_harvesteR' function. Optionally, it
is also possible to display the relative probabilities calculated by the
'hmpp' function.
}
\details{
The output files of the 'nodesigl_harvesteR' function are assumed to
be in the current working directory. If the 'import.probs' parameter is set
to TRUE, the output of the 'hmpp' function must also be in the
current working directory. "Hot nodes" are grouped into nestendness levels
based on their level of subordination to other "hot nodes" (i.e., from how
many other "hot nodes" they descend from). The branches of the tree are
colored according to the nestendness level of the “hot node” from which they
descend, with colors indicating an increasing nestendness from gray to
bright red. A different color is used for branches leading to known
medicinal species. High levels of nestedness should be indicative of the
aggregation of "hot nodes" and known medicinal species. However, being
dependent on the number of "hot nodes" present in the tree (and therefore on
the number of known medicinal species included in it), branch colors are not
comparable between different trees.\cr
Different symbols are used to highlight the differences between generic and
best ranking "hot nodes" (see the 'nodesigl_harvesteR' function for best
ranking "hot nodes" definition).
For the former, solid black circles are used, for the latter, black stars
are used by default, but user can choose different symbols. Plotting
functions from the 'ape' and 'phytools' packages are internally called by
this function and different graphical parameters can be used to customize
the final output (see Parameters for further details).
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# WARNING: to run the example a copy of the Phylocom binaries MUST be in the
# current temp directory. DO NOT RUN otherwise!
# Use tempdir() to find the current temp directory.

 library(ape)
 library(phytools)
 library(grDevices)
 library(graphics)

 # Generate a random tree with 100 tips and save it in the current temp folder
 set.seed(1)
 tr<-rtree(150)

 #Add dummy genera
 tr$tip.label<-paste0(rep(paste0("Genus", LETTERS[1:6]),each=25),"_",tr$tip.label)

 #Plot the tree
 tipcol<-c(rep("black",150))
 tipcol[c(10,11,12,13,14,25,26,27,28,29,30,32,35,41,42,43,46,47,48,55,56,57,58,
 59,60,101,102,103,105,106,108)]<-"red"
 #medicinal taxa in red
 plot(tr, tip.color =tipcol, cex=0.7, underscore=TRUE)

 #Show medicinal taxa labels
 tr$tip.label[c(10,11,12,13,14,25,26,27,28,29,30,32,35,41,42,43,46,47,48,55,56,
 57,58,59,60,101,102,103,105,106,108)]

 #Write the tree in the current temp folder
 tree<-write.tree(tr)
 treefile <-tempfile("tree", fileext = ".tree")
 cat(tree, file = treefile, sep = "\n")

 #Create a sample file and write it in the current temp folder
 sample_disease<-c("disease\t1\tGenusA_t28
disease\t1\tGenusA_t116
disease\t1\tGenusA_t37
disease\t1\tGenusA_t61
disease\t1\tGenusA_t113
disease\t1\tGenusA_t56
disease\t1\tGenusB_t118
disease\t1\tGenusB_t130
disease\t1\tGenusB_t50
disease\t1\tGenusB_t119
disease\t1\tGenusB_t111
disease\t1\tGenusB_t7
disease\t1\tGenusB_t150
disease\t1\tGenusB_t2
disease\t1\tGenusB_t48
disease\t1\tGenusB_t65
disease\t1\tGenusB_t129
disease\t1\tGenusB_t77
disease\t1\tGenusB_t147
disease\t1\tGenusC_t31
disease\t1\tGenusC_t35
disease\t1\tGenusC_t46
disease\t1\tGenusC_t125
disease\t1\tGenusC_t97
disease\t1\tGenusC_t137
disease\t1\tGenusE_t30
disease\t1\tGenusE_t101
disease\t1\tGenusE_t104
disease\t1\tGenusE_t80
disease\t1\tGenusE_t69
disease\t1\tGenusE_t106")
 samplefile<-tempfile("sample_disease", fileext = "")
 cat(sample_disease, file =samplefile, sep = "\n")

 #Get the temporary directory path and set it as current working directory
 path<-paste0(normalizePath(tempdir(), winslash = "/"),"/")
 setwd(path)

 #Output files are in the current temp folder (type 'path')
 tree4nodesigl(treefile, "example")
 nodesiglR(start=1,stop=100, s=samplefile,f="phylo_nobl_example")
 nodesigl_harvesteR(sample_file=basename(samplefile),tree="phylo_nobl_example",fract=5)
 files<-list.files(path, pattern ="IndependentHot")

 #Test
 hot_tree_painteR(inputTrees=files[2],inputNodes=files[1],indtr=1,branch.unit="Ma",
     autolabel=TRUE, arrow=TRUE, par.arrow=list(arrl=0.1,hedl=0.05,angle=25),
     par.legend=list(fsize=0.7), offset=3, lwd=3)

 #Tidy up. Remove the output files from the current temporary directory.
 toremove<-list.files(path, pattern = ".tree|phylo_|disease|Nodesigl_")
 file.remove(paste0(path,toremove))

\dontshow{\}) # examplesIf}
}
\references{
\itemize{
\item Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics
and evolutionary analyses in R.” Bioinformatics, 35, 526-528.
doi:10.1093/bioinformatics/bty633.\cr
\item Revell L (2024). “phytools 2.0: an updated R ecosystem for phylogenetic
comparative methods (and other things).” PeerJ, 12, e16505.
doi:10.7717/peerj.16505.\cr
\item Webb, C. O.; Ackerly, D. D. & Kembel, S. W. (2008) Phylocom: software for
the analysis of phylogenetic community structure and trait evolution.
Bioinformatics 24: 2098-2100.
\item Zecca, G., Toini, E., Labra, M, Grassi, F. (2025) Accelerating the
identification and the prioritisation of new plants with medicinal
potential: the  pm4mp R package.\cr
}
}
\seealso{
\code{\link[=nodesigl_harvesteR]{nodesigl_harvesteR()}}, \code{\link[=hmpp]{hmpp()}}
}
\author{
Giovanni Zecca, \email{giovanni.zecca@unimib.it} with the
contribution of Elisa Toini, \email{e.toini@campus.unimib.it}
}
