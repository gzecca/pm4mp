#'
#'@title Calculate the Hot Ancetry Score (HAS) for all taxa in the selected hot trees
#'
#'@description This function computes the Hot Ancestry Score (HAS) for all taxa
#'  included in the "hot trees" selected by the user defined 'cut_off'.
#'  The function uses the output files produced by the function
#'  'nodesigl_harvesteR()' as input for calculations.
#'
#'
#'@details The HAS approach is based on the "hot node" concept (Saslis-Lagoudakis et. al., 2012)
#'  and on the assumption that the regions of a hot tree where hot nodes are
#'  phylogenetically nested can be viewed as the clades in which the likelihood of finding
#'  new medicinal plants is highest.
#'  Cosequently,the greater the number of hierarchically nested hot nodes a species
#'  descends from, the higher its probability of being a medicinal taxon. The function
#'  calculates the number of ancestral hot nodes for each terminal taxon in the "hot trees"
#'  selected based on the parameter 'cut_off'. For each taxon found in a “hot tree”,
#'  the ratio between its HAS and the maximum HAS calculated in that tree is also
#'  calculated. Since this function does not need to keep track of branch descending
#'  from "hot nodes",it uses a different calculation method than the one implemented
#'  in the "hot_tree_painteR" function. The values ​​obtained from the two functions are
#'  however identical.
#'
#'@param sample_file Character (file name). A three columns, tab delimited file
#'  containing the output of the 'sample4nodesigl()' function. The same file
#'  used with the 'nodesiglR()' and 'nodesigl_harvesteR()' functions must be used.
#'
#'@param IndepTrees Character (file name). The text file named "IndependentHotTreesNewick_*"
#' outputted by the 'nodesigl_harvesteR()'function. This file includes the list of
#'  identified "independent hot trees" written in Newick format.
#'
#'@param IndepHotNodes Character (file name). The tab delimited file named
#' "IndependentHotNodesByTrees_*", produced by the 'nodesigl_harvesteR()'function.
#'
#'@param cut_off Integer. Defines the size that a tree must have to be included
#'  in the analysis. Only trees with the number of tips > cut_off are considered
#'  (default: cut_off = NULL).
#'
#'@return A csv file with four columns is saved in the current working
#'  "has_t\emph{i}_input_file_name.csv", where "t\emph{i}" is the
#'  "hot tree" analysed (with \emph{i} corresponding to the \emph{i}-th tree in
#'  the "IndepTrees" file), and "input_file_name" is the base name of the
#'  IndepHotNodes" input file, after that the prefix
#'  "IndependentHotNodesByTrees_" has been removed (i.e., it is the "*"
#'  part of the file name). Assuming that the ‘disout’ parameter of the function
#'  ‘sample4nodesigl()’ has been set to take this information into account,
#   the"input_file_name" tipically includes: the name of the disease of interest,
#'  the name of the reference phylogeny and the value of the "fract" parameter
#'  used by the function 'nodesigl_harvesteR'.
#'
#'  The column headings and related information displayed are as follows:\cr
#'  \emph{- species}, names of the species included in the analysed "hot
#'  tree";\cr
#'  \emph{- known_med_plants}, logical values. The column indicates
#'  whether a plant is an already known medicinal plant;\cr
#'  \emph{- HAS}, for each taxa included in the "hot tree", the computed Hot
#'  Ancestry Score is shown;\cr
#'  \emph{- HAS_ratio}, for each species of the ‘hot tree’, the ratio of the HAS
#'  calculated for that species to the maximum HAS calculated within the tree
#' is shown;\cr
#'
#'@author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#'
#'@references
#' - Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics
#'   and evolutionary analyses in R.” Bioinformatics, 35, 526-528.
#'   doi:10.1093/bioinformatics/bty633.
#' - Revell L.J.(2024). phytools 2.0: an updated R ecosystem for phylogenetic
#'   comparative methods (and other things). PeerJ 12:e16505
#'   https://doi.org/10.7717/peerj.16505
#' - Saslis-Lagoudakis CH, Savolainen V, Williamson EM, Forest F, Wagstaff SJ,
#'   Baral SR, Watson MF, Pendry CA, Hawkins JA. Phylogenies reveal predictive power
#'   of traditional medicine in bioprospecting. Proc Natl Acad Sci U S A. 2012
#'   Sep 25;109(39):15835-40. doi: 10.1073/pnas.1202242109.
#' - Webb, C. O.; Ackerly, D. D. & Kembel, S. W. (2008) Phylocom: software for the
#'   analysis of phylogenetic community structure and trait evolution.
#'   Bioinformatics 24: 2098-2100.
#' - Zecca, G., Toini, E., Labra, M, Grassi, F. (2025) Accelerating the
#'   identification and the prioritisation of new plants with medicinal
#'   potential: the  pm4mp R package.\cr
#'
#'@seealso [nodesigl_harvesteR()], [sample4nodesigl()], [hmpp()]
#'
#'@examplesIf interactive()
#' # WARNING: to run the example a copy of the Phylocom (Webb et al. 2008) binaries MUST be in the
#' # current temp directory. # DO NOT RUN otherwise!
#' # Use tempdir() to find the current temp directory.
#'
#' library(ape)
#' library(castor)
#' library(phytools)
#'
#' # Generate a random tree with 450 tips and save it in the current temp folder
#' set.seed(1) # do not change it
#' tr<-rtree(450)
#'
#' # Add dummy genera
#' tr$tip.label<-paste0(rep(paste0("Genus", LETTERS[1:18]),each=25),"_",tr$tip.label)
#'
#' # Plot the tree
#' tipcol<-c(rep("black",450))
#' tipcol[c(10,11,12,13,14,25,26,27,28,29,30,32,35,41,42,43,46,47,48,55,56,57,
#' 58,59,60,101,102,103,105,106,108,109,110,309, 313, 318,319,320,325,350,351,
#' 352,353,354,355,356,357, 358,362,364,368,374,380,399,400,403,407,408,411)]<-"red"
#' #known medicinal taxa in red
#' plot(tr,tip.color =tipcol, cex=0.7, underscore=TRUE)
#'
#' # Show medicinal taxa labels
#' tr$tip.label[c(10,11,12,13,14,25,26,27,28,29,30,32,35,41,42,43,46,47,48,55,56,
#' 57,58,59,60,101,102,103,105,106,108,109,110,309, 313, 318,319,320,325,350,351,
#' 352,353,354,355,356,357, 358,362,364,368,374,380,399,400,403,407,408,411)]
#'
#' # Write the tree in the current temp folder
#' tree<-write.tree(tr)
#' treefile <-tempfile("tree", fileext = ".tree")
#' cat(tree, file = treefile, sep = "\n")
#'
#' # Create a sample file and write it in the current temp folder
#'sample_disease<-c("disease\t1\tGenusA_t156
#'disease\t1\tGenusA_t439
#'disease\t1\tGenusA_t197
#'disease\t1\tGenusA_t220
#'disease\t1\tGenusA_t421
#'disease\t1\tGenusA_t185
#'disease\t1\tGenusB_t103
#'disease\t1\tGenusB_t383
#'disease\t1\tGenusB_t358
#'disease\t1\tGenusB_t253
#'disease\t1\tGenusB_t180
#'disease\t1\tGenusB_t221
#'disease\t1\tGenusB_t167
#'disease\t1\tGenusB_t141
#'disease\t1\tGenusB_t247
#'disease\t1\tGenusB_t229
#'disease\t1\tGenusB_t406
#'disease\t1\tGenusB_t368
#'disease\t1\tGenusB_t420
#'disease\t1\tGenusC_t398
#'disease\t1\tGenusC_t203
#'disease\t1\tGenusC_t104
#'disease\t1\tGenusC_t3
#'disease\t1\tGenusC_t179
#'disease\t1\tGenusC_t147
#'disease\t1\tGenusE_t334
#'disease\t1\tGenusE_t71
#'disease\t1\tGenusE_t315
#'disease\t1\tGenusE_t405
#'disease\t1\tGenusE_t58
#'disease\t1\tGenusE_t114
#'disease\t1\tGenusE_t320
#'disease\t1\tGenusE_t34
#'disease\t1\tGenusM_t67
#'disease\t1\tGenusM_t342
#'disease\t1\tGenusM_t267
#'disease\t1\tGenusM_t93
#'disease\t1\tGenusM_t433
#'disease\t1\tGenusM_t441
#'disease\t1\tGenusN_t237
#'disease\t1\tGenusO_t28
#'disease\t1\tGenusO_t12
#'disease\t1\tGenusO_t94
#'disease\t1\tGenusO_t126
#'disease\t1\tGenusO_t224
#'disease\t1\tGenusO_t380
#'disease\t1\tGenusO_t194
#'disease\t1\tGenusO_t250
#'disease\t1\tGenusO_t437
#'disease\t1\tGenusO_t440
#'disease\t1\tGenusO_t76
#'disease\t1\tGenusO_t120
#'disease\t1\tGenusP_t204
#'disease\t1\tGenusP_t29
#'disease\t1\tGenusP_t402
#'disease\t1\tGenusQ_t244
#'disease\t1\tGenusQ_t151
#'disease\t1\tGenusQ_t214
#'disease\t1\tGenusQ_t419")
#'
#' samplefile<- tempfile("sample_disease", fileext = "")
#' cat(sample_disease,file = samplefile, sep = "\n")
#'
#' # Get the temporary directory path and set it as current working directory
#' path<-paste0(normalizePath(tempdir(), winslash = "/"),"/")
#' setwd(path)
#'
#' # Output files are in the current temp folder (type 'path')
#' tree4nodesigl(treefile, "example")
#' # it may take a few minutes to complete
#' nodesiglR(start=1,stop=500, s=samplefile, f="phylo_nobl_example")
#' nodesigl_harvesteR(sample_file=basename(samplefile),tree="phylo_nobl_example",fract=5)
#' files<-list.files(path, pattern ="IndependentHot")
#'
#' #Test. Output files are in the current temp folder
#'
#' has(IndepTrees=files[2],IndepHotNodes=files[1],sample_file=basename(samplefile),cut_off=100)
#'
#' # Tidy up. Remove the output files from the current temporary directory.
#' toremove<-list.files(path, pattern = ".tree|phylo_|disease|Nodesigl_|has_")
#' file.remove(paste0(path,toremove))
#'
#' # You can delete the Phylocom copy yourself if you want
#'
#' @export
#'

has <-function(IndepTrees,IndepHotNodes,sample_file, cut_off=NULL){
  #import data
  HotTrees<-ape::read.tree(file=IndepTrees)
  if(inherits(HotTrees,"phylo")){
    HotTrees<-phytools::as.multiPhylo(HotTrees)}
  names(HotTrees)<-paste0("t",seq_along(HotTrees))
  dfHn<-utils::read.csv(IndepHotNodes)
  Medtips<-unlist(utils::read.delim(sample_file,header=FALSE,colClasses=c("NULL","NULL",NA)),
                  use.names = FALSE)
  #Extract root labels
  root_labels<-as.list(unique(dfHn$HotRoot))	#root_labels1<-lapply(HotTrees, \(x) x$node.label[castor::find_root(x)-length(x$tip.label)])
  #Apply cut_off
  if (!is.null(cut_off)){
    cutoff<-vapply(HotTrees, \(x) ape::Ntip(x)>cut_off,logical(1))
    HotTrees<-HotTrees[cutoff]
    root_labels<-root_labels[cutoff]
  }
  if(length(HotTrees)==0){
    stop("No tree has passed the set 'cut_off' threshold.")
  }
  #Extract hot node labels
  hn_labels<-lapply(root_labels, \(x) dfHn[dfHn$HotRoot == x,"HotNodes"])
  #extract tips from hot trees
  ht_tips<-lapply(HotTrees, \(x) x$tip.label)
  #Find subtrees at hot nodes and their sets of labels
  subtrees<-mapply(\(x,y) castor::get_subtrees_at_nodes(x,y)$subtrees, x=HotTrees, y=hn_labels, SIMPLIFY=FALSE)
  subtips<-vector("list", length(subtrees))
  for (i in seq_along(subtrees)){
    subtips[[i]]<- lapply(subtrees[[i]], \(x) x$tip.label)
  }
  #Check for presence of hot trees tips withini labels of hot subtrees i.e., computed the degree of nestendness of each tip
  presabs<-mapply(\(x,y) as.data.frame(lapply(x, \(z) y %in% z)), x=subtips, y=ht_tips,SIMPLIFY=FALSE)
  #absolute ranking
  has<-lapply(presabs, \(x) rowSums(x))
  has_ratio<-lapply(has, \(x) x/max(x))
  #Output
  kmp<-lapply(ht_tips, \(x) x %in% Medtips)
  df<- mapply(cbind, species=ht_tips,known_med_plants=kmp, HAS=has, HAS_ratio=has_ratio, SIMPLIFY=FALSE)
  dis_phy<- sub("IndependentHotNodesByTrees_", "", basename(IndepHotNodes))
  invisible(lapply(seq_along(df),
                   \(x) utils::write.csv(df[[x]], file=paste0("has_",names(HotTrees)[x],
                                                              "_", dis_phy, ".csv"),row.names = FALSE)))
}
