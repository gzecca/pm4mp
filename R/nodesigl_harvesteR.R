#'
#' @title filteR()
#' @usage NULL
#' @description Internal function to filter "hot nodes" obtained from multiple
#'   runs of Phylocom's 'nodesigl' analysis.
#' @details "Hot nodes" are filtered based on three criteria: 1) significant
#'   overabundance of medicinal plants in the terminal taxa distal to it; 2)
#'   stability across multiple runs; 3) presence of at least a minimum
#'   percentage of medicinal plants among the terminal taxa distal to it (i.e.,
#'   'fract' parameter).
#' @param mytree,fract See the 'nodesigl_harvesteR' function documentation for
#'   details.
#' @param nodesigl A list produced by the 'nodesigl_harvesteR' function, which
#'   includes the outputs of multiple 'nodesigl' analysis.
#' @return A list of R objects.
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#' @examplesIf interactive()
#' #NOT RUN
#' library(castor)
#' filtered<-filteR(mytree=mytree, nodesigl=nodesigl, fract=fract)
#' rm(filtered)
#' @keywords internal
#' @noRd
#'

filteR<-function(mytree,nodesigl, fract){
  plus<-lapply(nodesigl, \(x) x[x$mark == "+",]) #1st filter
  plusname<-lapply(plus, \(x) x$node_name)
  stableplus<-Reduce(intersect, plusname)				 #2nd filter
  counter<-data.frame(Node=mytree$node.label, TotTaxa=castor::count_tips_per_node(mytree))
  DF<-counter[counter$Node %in% stableplus,]
  colnames(DF)[1]<-"HotNodes"
  DF$NMed<-plus[[1]]$ntaxa[plus[[1]]$node_name %in% stableplus]
  DF$Fraction <-round((DF$NMed/DF$TotTaxa)*100, digits=3)
  DF<-DF[DF$Fraction > fract,] 			#3rd filter, adjust to the appropriate %value
  filtered<-list(plus=plus, DF=DF)
}

#'
#' @title bestRank()
#' @usage NULL
#' @description Internal function to identify the "hot nodes" with the best
#'   ranking obtained from multiple runs of Phylocom's 'nodesigl' analysis.
#' @details The "hot nodes" that always obtain the highest possible ranking in
#'   all replications are identified and collected.
#' @param DF A dataframe created by the 'filteR' function, which is modified and
#'   passed to the 'nodesigl_harvesteR' function.
#' @param plus A list produced by the 'filteR' function, which includes all
#'   internal node with significantly more medicinal plants than chance.
#' @return A list of R objects.
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#' @examplesIf interactive()
#' #NOT RUN
#' DF<-bestRank(plus=filtered$plus,DF=filtered$DF)
#' rm(DF)
#' @keywords internal
#' @noRd
#'

bestRank<-function(plus, DF){
  KeepMe<-lapply(plus, \(x) which(x$node_name %in% DF$HotNodes))
  ranks<-as.data.frame(lapply(seq_along(plus), \(x) plus[[x]][unlist(KeepMe[[x]]),"rank"]))
  DF$bestrank<-apply(ranks, 1,\(x) min(x)== max(x))
  DF<-DF[c("HotNodes","bestrank", "TotTaxa", "NMed", "Fraction")]
}

#'
#'@title Combine and summarise the results of multiple "nodesigl" analysis.
#'
#'@description This function summaries the output of multiple runs of the
#'  'nodesigl' analysis, such as those obtained from the 'nodesiglR' function.
#'
#'@details The function assumes that the files returned from the 'nodesiglR'
#'  function are in the current working directory and are named "Nodesigl_1"...
#'  "Nodesigl_n", where "n" is the number of replicates. "Hot nodes" identified in
#'  different 'nodesigl' analysis replicates are filtered based on three
#'  criteria: 1) significant overabundance of medicinal plants in the terminal
#'  taxa distal to it (i.e., only nodes marked with the "+" symbol are
#'  retained); 2) stability across multiple replicates (i.e., only nodes present in
#'  all replicates are retained); 3) presence of at least a minimum percentage of
#'  medicinal plants among the terminal taxa distal to it (i.e., only nodes
#'  whose set of descending tips contains at least a minimum percentage of
#'  medicinal plants, as defined by the user, are retained). Additionally, “hot
#'  nodes” that received the highest possible ranking across all replicates are
#'  identified and printed in the output.\cr
#'
#'  Among all "hot nodes" retained after filtering, independent subsets of
#'  nested nodes are identified. For each independent subset the "hot node"
#'  which is the MRCA of all other nodes in the subset is taken as the root of a
#'  subtree of the original reference phylogeny, thus defining a new set of
#'  "independent hot trees". These "independent hot trees" contain distinct
#'  subsets of nested "hot nodes" and thus different sets of species. They
#'  represent clades of the original reference phylogeny where known medicinal
#'  plants and stable hot nodes are clumped, thus representing the most
#'  promising clades for the search for new medicinal species.
#'
#'@param sample_file A three columns, tab delimited file containing the output
#'  of the 'sample4nodesigl' function. The same file used with the 'nodesiglR'
#'  function must be used.
#'
#'@param tree A reference phylogeny in Newick format containing the output of the
#'  'tree4nodesigl' function. The same file used with the 'nodesiglR' function
#'  MUST be used. In particular, it is essential that the internal nodes of the
#'  tree are named in the same way as in the previous analysis.
#'
#'@param fract A number between 0 and 100, which defines the threshold
#'  percentage of medicinal plants among the terminal taxa distal to a “hot
#'  node", above which it is retained by the filtering process. Please note that
#'  only nodes with a percentage of medicinal plants strictly greater than
#'  'fract' are retained. Default value is 5. Increasing this value will result
#'  in smaller (and generally more numerous) independent trees, while decreasing
#'  the number will result in larger (and generally fewer) independent trees.
#'
#'@return Three files are saved in the current working directory:
#' \itemize{
#' \item "IndependentHotNodesByTrees": a csv file with eight columns containing
#' several information:\cr
#' \emph{- IndepTree}, a progressive naming of the "independent hot trees"
#' identified by the function (i.e., t1,t2,t3...);\cr
#' \emph{- HotRoot}, the names of the "hot nodes" taken as the roots of the
#' "independent hot trees";\cr
#' \emph{- HotNodes},the names of the nested "hot nodes" included in each
#' "independent hot trees" (including the root of the tree);\cr
#' \emph{- bestrank}, a column of logic values with TRUE identifying those “hot nodes”
#' that always received the highest possible score (i.e., ranking) across all
#' replicates;\cr
#' \emph{- TotTaxa}, the number of the terminal taxa distal to each "hot node"
#' (i.e., the number of tips descending from it);\cr
#' \emph{- NMed},the number of the medicinal taxa distal to each "hot node";\cr
#' \emph{- Fraction}, the percentage (100*NMed/TotTaxa) of medicinal taxa descending
#' from each "hot node"(up to three decimal digits are calculated);\cr
#' \emph{- MedTaxa}, medicinal taxa included in each "independent hot tree" are
#' listed with names separated by a forward slash, preceded and followed by white
#' space. Taxa names are written only in the first row of each tree.\cr
#' \item "GeneraByTreeRoots": a csv file with a variable number of columns equal
#' to the number of "independent hot trees" identified by the function.
#' The first line contains the root names of the "independent hot trees". For
#' each "hot root" the list of the genera included in the corresponding
#' "independent hot tree" is given.\cr
#' \item "IndependentHotTreesNewick": a text file including the list of
#' "independent hot trees" identified by the function. Trees are written in Newick
#' format and are listed in the order given by the "IndepTree" column of the
#' "IndependentHotNodesByTrees" file.
#' }
#'  To prevent the results from being overwritten two suffixes are added to the
#'  end of the output file names. The first suffix indicating the disease of
#'  interest is taken directly from the name of the sample_file and it is equal
#'  to the 'disout' parameter of the 'sample4nodesigl' function. The second
#'  suffix is information about the date and time the files were created.
#'
#'@author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#'
#'@references
#' - Louca S, Doebeli M (2017). “Efficient comparative phylogenetics on large
#'   trees.” Bioinformatics. doi:10.1093/bioinformatics/btx701.
#' - Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics
#'   and evolutionary analyses in R.” Bioinformatics, 35, 526-528.
#'   doi:10.1093/bioinformatics/bty633.
#' - Webb, C. O.; Ackerly, D. D. & Kembel, S. W. (2008) Phylocom: software for the
#'   analysis of phylogenetic community structure and trait evolution.
#'   Bioinformatics 24: 2098-2100.
#'- Zecca, G., Toini, E., Labra, M, Grassi, F. (2025) Accelerating the
#'   identification and the prioritisation of new plants with medicinal
#'   potential: the  pm4mp R package.\cr
#'
#'@seealso [clean_and_match()], [sample4nodesigl()], [tree4nodesigl()]
#'
#'@examplesIf interactive()
#'# WARNING: to run the example a copy of the Phylocom binaries MUST be in the
#'# current temp directory. DO NOT RUN otherwise!
#'# Use tempdir() to find the current temp directory.
#'
#'  library(ape)
#'  library(castor)
#'
#'  # Generate a random tree with 100 tips and save it in the current temp
#'  # folder
#'  set.seed(1)
#'  tr<-rtree(100)
#'
#'  #Add dummy genera
#'  tr$tip.label<-paste0(c(rep("A",25),rep("B",15),rep("C",60)),"_",tr$tip.label)
#'
#'  #Plot the tree
#'  tipcol<-c(rep("black",100))
#'  tipcol[c(10,11,12,25,26,27,28,29,30,41,42,43)]<-"red"
#'
#'  #medicinal taxa in red
#'  plot(tr, tip.color =tipcol)
#'
#'  #Show medicinal taxa labels
#'  tr$tip.label[c(10,11,12,25,26,27,28,29,30,41,42,43)]
#'
#'  #Write the tree in the current temp folder
#'  tree<-write.tree(tr)
#'  treefile <-tempfile("tree", fileext = ".tree")
#'  cat(tree, file = treefile, sep = "\n")
#'
#'  #Create a sample file and write it in the current temp folder
#'  sample_disease<-c("disease\t1\tA_t99
#'disease\t1\tA_t33
#'disease\t1\tA_t45
#'disease\t1\tB_t59
#'disease\t1\tB_t26
#'disease\t1\tB_t15
#'disease\t1\tB_t58
#'disease\t1\tB_t95
#'disease\t1\tC_t19
#'disease\t1\tC_t74
#'disease\t1\tC_t35
#'disease\t1\tA_t43")
#'
#'  samplefile<- tempfile("sample_disease", fileext = "")
#'  cat(sample_disease, file = samplefile, sep = "\n")
#'
#'  #Get the temporary directory path and set it as current working directory
#'  path<-paste0(normalizePath(tempdir(), winslash = "/"),"/")
#'  setwd(path)
#'
#'  #Test: output files are in the current temp folder (type 'path')
#'  tree4nodesigl(treefile, "example")
#'  nodesiglR(start=1,stop=50, s=samplefile, f="phylo_nobl_example")
#'  nodesigl_harvesteR(sample_file=basename(samplefile),tree="phylo_nobl_example",fract=25)
#'
#'  #Tidy up. Remove the output files from the current temporary directory.
#'  toremove<-list.files(path, pattern = ".tree|phylo_|disease|Nodesigl_")
#'  file.remove(paste0(path,toremove))
#'
#'@export
#'

nodesigl_harvesteR<-function(sample_file, tree, fract=5){
  mytree<-ape::read.tree(tree)
  filenames<-list.files(path = ".", pattern ="Nodesigl_")
  nodesigl<- lapply(filenames, \(x) utils::read.delim(x,skip=1, header=FALSE, strip.white = TRUE,
                    col.names=c("plot","node","node_name","ntaxa","median","rank","sig","mark")))
  #Filtering results:keeping only stable "hot nodes" with appropriate number of med plants
  filtered<-filteR(mytree=mytree, nodesigl=nodesigl, fract=fract)
  #Best ranking nodes
  DF<-bestRank(plus=filtered$plus, DF=filtered$DF)
  #Preparing output files
  abovethreshold<-DF$HotNodes
  samp<-utils::read.table(sample_file, colClasses=c("NULL","NULL", NA))
  IndepTree<-c()
  HotRoot<- c()
  Genera<-list()
  MedTaxa<- c()
  Trees<-list()

#deb0<-list()

  n<-1
  while (length(abovethreshold) >0){
    tr<-castor::get_subtree_at_node(mytree, abovethreshold[1])
    gen<-sort(unique(gsub("_.*$", "",tr$subtree$tip.label)))

#deb1<-intersect(tr$subtree$tip.label, samp$V3)
#deb0[[n]]<-deb1

    DR<- paste(intersect(tr$subtree$tip.label, samp$V3), collapse=" / ")
    included<-abovethreshold %in% tr$subtree$node.label
    subtrue<-abovethreshold[included]
    IndepTree<-c(IndepTree, paste0("t",n),rep(".     ", (length(subtrue)-1)))
    HotRoot<-c(HotRoot,rep(abovethreshold[1], length(subtrue)))
    Genera[[abovethreshold[1]]]<-gen
    MedTaxa<-c(MedTaxa,DR,rep("", (length(subtrue)-1)))
    Trees[[abovethreshold[1]]]<-tr$subtree
    n<-n+1
    abovethreshold<- setdiff(abovethreshold, subtrue)
  }
  DF<-cbind(IndepTree,HotRoot,DF,MedTaxa)
  maxlength<-max(sapply(Genera, \(x) (length(x))))
  Genera<-as.data.frame(lapply(Genera, \(x) c(x, rep("", maxlength-length(x)))))
  #Writing output files
  disname<-sub("sample_","",sample_file)
  refphylo<-sub("phylo_bl_||phylo_nobl_","",tree)
  # stime<-sub(" ", "_", gsub("-|\\:","", round(Sys.time())))
  utils::write.csv(DF,
                   file = paste0("IndependentHotNodesByTrees_",disname,"_",refphylo,"_f",fract),
                   quote=FALSE, row.names=FALSE)
  utils::write.csv(Genera,
                   file = paste0("GeneraByTreeRoots_",disname,"_",refphylo,"_f",fract),
                   quote = FALSE, row.names = FALSE)
  invisible(lapply(Trees, \(x) ape::write.tree(x,
                  file =paste0("IndependentHotTreesNewick_",disname,"_",refphylo,"_f",fract),
                  append = TRUE, digits = 10)))
}

