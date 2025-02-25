#'@title Visual display of computed probabilities
#'
#'@description This function allows to graphically display the probabilities
#'  computed using the 'hmpp()' function. If different methods or
#'  settings were used to calculate probabilities for the same taxa, their
#'  results can be compared graphically.
#'
#'@details This function accepts the output of a single or multiple analysis of
#'  the 'hmpp()' function and displays the probability of "state 1"
#'  (i.e., of being a medicinal plant) computed for all the taxa of interest.
#'  Different colours are used to show results from different files. Different
#'  visualization options are available (see the parameter 'format' below).
#'
#'@param file_names Character. The name of a single file or a vector including
#'  different file names. Input file(s) must be the csv file(s) with eight
#'  columns produced by the 'hmpp()' function.
#'
#'@param legend_names Character. If only one file is loaded, the name to use in
#'  the legend for the analysis to be displayed. When multiple files are loaded
#'  at the same time, a vector containing the names to use with the different
#'  analyses. The length of 'legend_names' must be equal to the length of
#'  'file_names' parameter.
#'
#'@param format Character. Determines the format of the graphical display.
#'  Currently three options are available: \emph{h} to chose "histogram" as
#'  visual representation of data; \emph{kd} to chose "kernel density" as visual
#'  representation of data; \emph{comp} to combine both "histogram" and "kernel
#'  density" options into a single plot.
#'
#'@return A plot showing the distribution of the computed probabilities is
#'  displayed on the screen device.
#'
#'@author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#'
#'@references
#' - Wickham H. ggplot2: Elegant Graphics for Data Analysis.
#'   Springer-Verlag New York, 2016.
#' - Zecca, G., Toini, E., Labra, M, Grassi, F. (2025) Accelerating the
#'   identification and the prioritisation of new plants with medicinal
#'   potential: the  pm4mp R package.\cr
#'
#'@examplesIf interactive()
#' # WARNING: to run the example a copy of the Phylocom binaries MUST be in the
#' # current temp directory. DO NOT RUN otherwise!
#' # Use tempdir() to find the current temp directory.
#'
#'  library(ape)
#'  library(castor)
#'  library(factoextra)
#'  library(cluster)
#'  library(stats)
#'  library(phytools)
#'  library(ggplot2)
#'  library(ggpubr)
#'
#' # Generate a random tree with 100 tips and save it in the current temp folder.
#'  set.seed(1)
#'  tr<-rtree(450)
#'
#' # Add dummy genera
#'  tr$tip.label<-paste0(rep(paste0("Genus", LETTERS[1:18]),each=25),"_",tr$tip.label)
#'
#' # Plot the tree
#'  tipcol<-c(rep("black",450))
#'  tipcol[c(10,11,12,13,14,25,26,27,28,29,30,32,35,41,42,43,46,47,48,55,56,57,58,
#'  59,60,101,102,103,105,106,108,109,110,309, 313, 318,319,320,325,350,351,352,
#'  353,354,355,356,357, 358,362,364,368,374,380,399,400,403,407,408,411)]<-"red"
#'
#' # known medicinal taxa in red
#'  plot(tr,tip.color =tipcol, cex=0.7, underscore=TRUE)
#'
#' # Show medicinal taxa labels
#'  tr$tip.label[c(10,11,12,13,14,25,26,27,28,29,30,32,35,41,42,43,46,47,48,55,56,
#'  57,58,59,60,101,102,103,105,106,108,109,110,309, 313, 318,319,320,325,350,351,
#'  352,353,354,355,356,357, 358,362,364,368,374,380,399,400,403,407,408,411)]
#'
#' # Write the tree in the current temp folder
#'  tree<-write.tree(tr)
#'  treefile <-tempfile("tree", fileext = ".tree")
#'  cat(tree, file = treefile, sep = "\n")
#'
#' # Create a sample file and write it in the current temp folder
#'  sample_disease<-c("disease\t1\tGenusA_t156
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
#'  samplefile<- tempfile("sample_disease", fileext = "")
#'  cat(sample_disease,file = samplefile, sep = "\n")
#'
#' # Get the temporary directory path and set it as current working directory
#'  path<-paste0(normalizePath(tempdir(), winslash = "/"),"/")
#'  setwd(path)
#'
#' # Output files are in the current temp folder (type 'path')
#'  tree4nodesigl(treefile, "example")
#'  nodesiglR(start=1,stop=500, s=samplefile,f="phylo_nobl_example")
#' # It may take a few minutes to complete
#'  nodesigl_harvesteR(sample_file=basename(samplefile),tree="phylo_nobl_example", fract=5)
#'  files<-list.files(path, pattern ="IndependentHot")
#'
#' #Test. Output files are in the current temp folder
#' # Single threshold set at: 0.70
#'  hmpp(reftree="phylo_nobl_example",sample_file=basename(samplefile),
#'                 IndepTrees=files[2],method="single_thr", thr_level=0.7,
#'                 cut_off=100,min_revealed = 2,max_STE=0.25,seed=123)
#'
#' # Multiple (5) thresholds set at: 0.75, 0.80, 0.85, 0.90, 0.95
#'  hmpp(reftree="phylo_nobl_example",sample_file=basename(samplefile),
#'                 IndepTrees=files[2],method="multi_thr", thr_level=c(0.7,0.95),
#'                 by=0.05, cut_off=100,min_revealed= 2,max_STE =0.25,seed=123)
#'
#' # kmeans method testing a maximum of 5 clusters
#'  hmpp(reftree="phylo_nobl_example",sample_file=basename(samplefile),
#'                 IndepTrees=files[2],method="kmeans",cut_off=100,kmax=5,
#'                 min_revealed = 2,max_STE =0.25,seed=123)
#'
#' # PAM method testing a maximum of 20 clusters
#'  hmpp(reftree="phylo_nobl_example",sample_file=basename(samplefile),
#'                 IndepTrees=files[2],method="pam", cut_off=100,kmax =20,
#'                 min_revealed = 2,max_STE =0.25,seed=123)
#'
#' # Collect the outputs releted to the "hot tree" t1, obtained using
#' # different parameter settings
#'  t1<-list.files(path,pattern = "hspb_t1_")
#'
#' # Plot results
#'
#' #Histograms
#'  visual_comp(file_names= t1[1], legend_names="kmeans, k = 2",format="h")
#'
#'  visual_comp(file_names= c(t1[1], t1[2]), legend_names=c("kmeans, k =2",
#'            "multi_thr: 0.7 - 0.95, by= 0.05"), format="h")
#'
#' #Kernel density
#'  visual_comp(file_names= c(t1[1], t1[2]), legend_names=c("kmeans, k =2",
#'            "multi_thr: 0.7 - 0.95, by= 0.05"), format="kd")
#'
#' # Combined
#'  visual_comp(file_names= t1, legend_names=c("kmeans, k = 2",
#'              "multi_thr: 0.7 - 0.95, by= 0.05","pam, k=11",
#'               "single_thr: 0.7"), format="comb")
#'
#' # Tidy up. Remove the output files from the current temporary directory.
#'  toremove<-list.files(path, pattern = ".tree|phylo_|disease|Nodesigl_|hspb_")
#'  file.remove(paste0(path,toremove))
#'
#'@export
#'

visual_comp<-function(file_names, legend_names, format="h"){
  hsp_p1s<-sapply(file_names,
          \(x) utils::read.csv(x, colClasses =c("NULL","NULL",NA, "NULL","NULL","NULL","NULL","NULL")))
  method<-unlist(lapply(legend_names, \(x) rep(x, length(hsp_p1s[[1]]))))
  p<-unlist(hsp_p1s)
  datfr<-data.frame(p,method)
  tree<-sub("^.*_(t\\d+)_.*$","\\1", file_names[1])
  switch(format,
    h={ ggplot2::ggplot(datfr,  ggplot2::aes(x=p, fill=method, col=method)) +
       ggplot2::geom_histogram(alpha=0.4, binwidth = 0.05) +
       ggplot2::theme_minimal()+ ggplot2::ylim(0, length(hsp_p1s[[1]])) +
       ggplot2::facet_wrap(~method,strip.position = "top")+
       ggplot2::theme(strip.text = ggplot2::element_blank())+
       ggplot2::labs(x="P(state 1)", y="Count")},
    kd={ggplot2::ggplot(datfr,ggplot2::aes(x=p, fill=method)) +
       ggplot2::geom_density(alpha=0.4,adjust = 1/2) +
       ggplot2::theme_minimal()+
       ggplot2::labs(title=paste0("Kernel density plot of the estimated probability of state 1 - ",
                           "tree ",tree),x="Probability", y="Density")},
    comb={ggpubr::ggarrange(
       ggplot2::ggplot(datfr,  ggplot2::aes(x=p, fill=method, col=method))+
       ggplot2::geom_histogram(alpha=0.4, binwidth = 0.05) +
       ggplot2::theme_minimal()+ ggplot2::ylim(0, length(hsp_p1s[[1]]))+
       ggplot2::facet_wrap(~method,strip.position = "top")+
       ggplot2::theme(strip.text = ggplot2::element_blank())+
       ggplot2::labs(x="P(state 1)", y="Count"),

       ggplot2::ggplot(datfr,ggplot2::aes(x=p, fill=method))+
       ggplot2::geom_density(alpha=0.4,adjust = 1/2)+ggplot2::theme_minimal()+
       ggplot2::labs(title=paste0("Kernel density plot of the estimated probability of state 1 - ",
                   "tree ",tree),x="Probability", y="Density"),ncol = 1, nrow = 2)},
    stop("Invalid 'format' value")
  )
}

