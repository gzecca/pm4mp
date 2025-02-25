#'
#' @title to_root()
#' @usage NULL
#' @description Internal function used to root each one of the selected "hot
#'   trees" with an external taxon not belonging to the known medicinal plants.
#' @details the specified "hot tree"is rooted using an external taxa randomly
#'   chosen among one of the nearest non-medicinal taxa. Setting seeds (see
#'   'hmpp()' function) ensure reproducibility of the analysis.
#' @param HotTrees,Btree,Medtips internal parameters defined inside the
#'   'hmpp()' function.
#' @return A list of R objects.
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#' @examplesIf interactive()
#' #NOT RUN
#' library(castor)
#' to_root(HotTrees=HotTrees,Btree=Btree,Medtips=Medtips)
#' rm(to_root_out)
#' @keywords internal
#' @noRd
#'

to_root<-function(HotTrees,Btree,Medtips){
  backnodes<-lapply(HotTrees,
                    \(x) castor::get_ancestral_nodes(Btree, descendants= x$node.label[1], Nsplits=1))
  backtrees<-lapply(backnodes,
                    \(x) castor::get_subtree_at_node(tree= Btree, node=Btree$node.label[x])$subtree)
  NewTips<- lapply(seq_along(backtrees),
                   \(x) setdiff(backtrees[[x]]$tip.label, HotTrees[[x]]$tip.label))
  NoMedNewTips<-lapply(NewTips, \(x) setdiff(x, Medtips))
  Rs2Keep<-lapply(seq_along(NoMedNewTips),
                  \(x) which(NewTips[[x]] %in% sample(NoMedNewTips[[x]], 1)))
  outgrnames<-lapply(seq_along(Rs2Keep),\(x) NewTips[[x]][Rs2Keep[[x]]])
  Tips2omit<-lapply(seq_along(Rs2Keep),\(x) NewTips[[x]][- Rs2Keep[[x]]])
  RHotTrees<-lapply(seq_along(Tips2omit),
                    \(x) castor::get_subtree_with_tips(backtrees[[x]],
                                                       omit_tips = Tips2omit[[x]],
                                                       collapse_monofurcations = TRUE,
                                                       force_keep_root = TRUE)$subtree)
  to_root_out<-list(RHotTrees=RHotTrees, outgrnames=outgrnames)
}

#'
#' @title min_dist()
#' @usage NULL
#' @description Internal function used to compute the phylogenetic distance
#'   between each plant not belonging to the known medicinal plants and its
#'   nearest medicinal taxon.
#' @details all distances between taxa not belonging to medicinal plants and the
#'   known medicinal plants included in the selected "hot trees" are calculated.
#'   Then for each taxon not included in medicinal plants, the minimum distance
#'   is taken for further computations.
#' @param HotTrees,Btree,Medtips internal R object defined inside
#'   the'hmpp()' function.
#' @return A list of R objects.
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#' @examplesIf interactive()
#' #NOT RUN
#' library(castor)
#' to_root(HotTrees=HotTrees,Btree=Btree,Medtips=Medtips)
#' rm(min_dist_out)
#' @keywords internal
#' @noRd
#'

min_dist<-function(HotTrees,Btree,Medtips){
  TreesTipsList<-lapply(HotTrees, \(x) x$tip.label)
  TreesHotTipsList<- lapply (TreesTipsList, \(x) x[x %in% Medtips])
  TreesColdTipsList<- lapply (TreesTipsList, \(x) x[!x %in% Medtips])
  PhyDistAll<-castor::get_all_pairwise_distances(Btree,only_clades =unlist(TreesTipsList),
                                                 as_edge_counts=FALSE, check_input=TRUE)
  colnames(PhyDistAll)<-unlist(TreesTipsList)
  rownames(PhyDistAll)<-unlist(TreesTipsList)
  PhyDistList<-lapply(seq_along(HotTrees),
                      \(x) PhyDistAll[TreesHotTipsList[[x]], TreesColdTipsList[[x]]])
  MinPhyDistList<-lapply(PhyDistList, \(x) apply(x, 2, min))
  min_dist_out<-list(TreesHotTipsList=TreesHotTipsList,
                     TreesColdTipsList=TreesColdTipsList,
                     MinPhyDistList=MinPhyDistList)
}

#'
#' @title mthds()
#' @usage NULL
#' @description Internal function used to select the method applied to identify
#'   putative low medicinal potential plants based on a phylogenetic distance
#'   threshold from the nearest observed medicinal taxon.
#' @details The function uses the output of the 'min_dist()' function and
#'   applies one of four available methods to identify all tips that have a
#'   distance to the nearest medicinal taxa greater than a certain threshold(s).
#'   Thresholds are either specified by the user or automatically identified by
#'   the function.
#' @param method,thr_level,kmax See the 'hmpp()' function
#'   documentation for details.
#' @param MinPhyDistList,TreesColdTipsList Internal R objects defined inside the
#'   'min_dist()' function.
#' @return A list of R objects.
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#' @examplesIf interactive()
#'   #NOT RUN
#'   library(stats)
#'   library(cluster)
#'   library(factoextra)
#'   mthds(method,thr_level,kmax,MinPhyDistList,TreesColdTipsList)
#'   rm(mthds_out)
#' @keywords internal
#' @noRd
#'

mthds<- function(method,thr_level,kmax, MinPhyDistList, TreesColdTipsList){
  switch(method,
         # single number: the same threshold apply to all trees
         # numeric vector: different trees can have different threshold
         single_thr={
           thresholds<-mapply(\(x,y) stats::quantile(x, prob=y),
                              x=MinPhyDistList, y=thr_level, SIMPLIFY=FALSE)
           abvT<- mapply(\(x,y) x>y, x=MinPhyDistList, y=thresholds, SIMPLIFY=FALSE)
           TipsAboveT<-mapply(\(x,y) x[y], x=TreesColdTipsList, y=abvT, SIMPLIFY=FALSE)
           mthds_out<-list(TipsAboveT=TipsAboveT)},

         # 'thr_level' is a numeric vector. All Thresholds are applied to all trees
         multi_thr={
           thresholds<-lapply(thr_level,\(y) mapply(\(x,y) stats::quantile(x, prob=y),
                                                    x=MinPhyDistList, y=y, SIMPLIFY=TRUE))
           names(thresholds)<-thr_level
           abvT<- lapply(thresholds, \(y) mapply(\(x,y) x>y, x=MinPhyDistList, y=y,
                                                 SIMPLIFY=FALSE))
           TipsAboveT<- lapply(abvT, \(y) mapply(\(x,y) x[y], x=TreesColdTipsList, y=y,
                                                 SIMPLIFY=FALSE))
           mthds_out<-list(TipsAboveT=TipsAboveT)},

         kmeans={
           MinPhyScaled<-lapply(MinPhyDistList, \(x) as.data.frame(scale(x)))
           ks<-lapply(MinPhyScaled,
                      \(x) factoextra::fviz_nbclust(x, stats::kmeans, method = "silhouette",k.max = kmax))
           bestk<-lapply(ks, \(x) as.numeric(x$data[x$data$y == max(x$data$y), "clusters"]))
           clust<-mapply(\(x,y) stats::kmeans(x=x, centers=y, iter.max = 50, nstart = 100),
                         x=MinPhyScaled, y=bestk, SIMPLIFY=FALSE)
           max_cc<-lapply(clust, \(x) which(x$centers == max(x$centers)))
           TipsAboveT<-mapply(\(x,y) names(x$cluster[x$cluster == y]),
                              x=clust, y=max_cc, SIMPLIFY=FALSE)
           mthds_out<-list(TipsAboveT=TipsAboveT,bestk=bestk,ks=ks)},

         pam={
           MinPhyScaled<-lapply(MinPhyDistList, \(x) as.data.frame(scale(x)))
           ks<-lapply(MinPhyScaled,
                      \(x) factoextra::fviz_nbclust(x,  cluster::pam, method = "silhouette",k.max = kmax))
           bestk<-lapply(ks, \(x) as.numeric(x$data[x$data$y == max(x$data$y), "clusters"]))
           clust<-mapply(\(x,y) cluster::pam(x=x, k=y, variant = "original"),
                         x=MinPhyScaled, y=bestk, SIMPLIFY=FALSE)
           max_cm<-lapply(clust, \(x) which(x$medoids == max(x$medoids)))
           TipsAboveT<-mapply(\(x,y) names(x$clustering[x$clustering == y]),
                              x=clust, y=max_cm, SIMPLIFY=FALSE)
           mthds_out<-list(TipsAboveT=TipsAboveT,bestk=bestk,ks=ks)},
         stop("Invalid 'method' value."))
}

#'
#' @title one4all()
#' @usage NULL
#' @description Internal function used to compute tips probabilities when
#'   methods 'single_thr', 'kmeans' or 'pam' are used.
#' @details Prior states and prior probabilities of being 'state one' are
#'   assigned to each tip. Then tips probabilities are computed using the hidden
#'   state prediction approach implemented in the 'castor' R package.
#' @param MinPhyDistList,TipsAboveT,RHotTrees, R object defined within the
#'   'mthds()','min_dist()' and 'to_root()' functions.
#' @param min_revealed,max_STE parameters set by the user and internally used by
#'   the castor::hsp_binomial function.
#' @param EmCDFs,tip_states,st1_probs R object defined inside the main function
#'   'hmpp()'.
#' @param cf, correction factor used to correct the prior probability of measuring the
#'    'state 1' applied to each tip(default =0.1).Ignored if mthds is "kmeans" or "pam".
#' @return A list of R objects.
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#' @examplesIf interactive()
#' #NOT RUN
#' library(castor)
#' library(ape)
#' one4all(MinPhyDistList, TipsAboveT, EmCDFs, RHotTrees, tip_states,
#'         st1_probs, min_revealed, max_STE)
#' rm(one4all_out)
#' @keywords internal
#' @noRd
#'

one4all<-function(MinPhyDistList,TipsAboveT,EmCDFs, RHotTrees,tip_states,
                  st1_probs,min_revealed,max_STE, cf){
  MinDistAboveT <- mapply(\(x,y) x[y], x=MinPhyDistList, y=TipsAboveT, SIMPLIFY=FALSE)
  # Find prior probabilities above threshold & apply correction factor
  probsAboveT<-mapply(\(x,y) y(x)-cf, x=MinDistAboveT,  y=EmCDFs, SIMPLIFY=FALSE)
  for(i in seq_along(RHotTrees)){
    names(probsAboveT[[i]])<-names(MinDistAboveT [[i]])
    tip_states[[i]][TipsAboveT[[i]]]<-2					                # Complete tip states
    st1_probs[[i]][names(probsAboveT[[i]]),1]<-probsAboveT[[i]] # Complete 'state 1' pp
  }
  #Hidden state prediction & output preparation
  hspb<-lapply(seq_along(RHotTrees),
               \(x) castor::hsp_binomial(tree=RHotTrees[[x]],
                                         tip_states=tip_states[[x]],
                                         reveal_probs =NULL,
                                         state1_probs = st1_probs[[x]],
                                         min_revealed = min_revealed,
                                         max_STE = max_STE,
                                         check_input = TRUE))
  p1<- lapply(seq_along(hspb),
              \(x) round(hspb[[x]]$P1[1:ape::Ntip(RHotTrees[[x]])],digits=5))
  ste<-lapply(seq_along(hspb),
              \(x) round(hspb[[x]]$STE[1:ape::Ntip(RHotTrees[[x]])],digits=5))
  st1_probs<-lapply(st1_probs, \(x) apply(x,2, round, digits=5))
  st1_p<-lapply(st1_probs, \(x) paste0("c(", apply(x, 1, paste,collapse=","), ")"))
  one4all_out<-list(p1=p1, ste=ste, tip_states=tip_states, st1_p=st1_p)
}

#'
#' @title all4all()
#' @usage NULL
#' @description Internal function used to compute tips probabilities when th
#'   method 'multi_thr' is used.
#' @details Prior states and prior probabilities of being in 'state 1' are
#'   computed for each threshold specified by the user and are assigned to tips.
#'   Then tips probabilities are computed for each threshold using the hidden
#'   state prediction approach implemented in the 'castor' R package. Finally,
#'   the mean probability value is computed for each tip averaging over all
#'   values computed under different thresholds.
#' @param MinPhyDistList,TipsAboveT,RHotTrees, R object defined within the
#'   'mthds()','min_dist()' and 'to_root()' functions.
#' @param min_revealed,max_STE parameters set by the user and internally used by
#'   the castor::hsp_binomial function.
#' @param EmCDFs,tip_states,st1_probs R object defined inside the main function
#'   'hmpp()'.
#' @param thr_level parameter set by the user in the'hmpp()' function
#'   that define the thresholds to be used.
#' @param cf, correction factor used to correct the prior probability of measuring the
#'    'state 1' applied to each tip (default =0.1).Ignored if mthds is "kmeans" or "pam".
#' @return A list of R objects.
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#' @examplesIf interactive()
#' #NOT RUN
#' library(castor)
#' library(ape)
#' all4all(TipsAboveT, MinPhyDistList, EmCDFs, RHotTrees, thr_level,
#'         tip_states, st1_probs, min_revealed, max_STE)
#' rm(all4all_out)
#' @keywords internal
#' @noRd
#'

all4all<-function(TipsAboveT,MinPhyDistList,EmCDFs,RHotTrees,thr_level,tip_states,
                  st1_probs,min_revealed,max_STE, cf){
  MinDistAboveT <-lapply(TipsAboveT,
                         \(y) mapply(\(x,y) x[y], x=MinPhyDistList, y=y, SIMPLIFY=FALSE))
  probsAboveT<-lapply(MinDistAboveT,
                      \(x) mapply(\(x,y) y(x)-cf, x=x,  y=EmCDFs, SIMPLIFY=FALSE))
  for(i in seq_along(probsAboveT)) {
    for (k in seq_along(RHotTrees)){
      names(probsAboveT[[i]][[k]])<-names(MinDistAboveT [[i]][[k]])
    }
  }
  TS<-stats::setNames(vector("list", length(thr_level)), thr_level)
  tip_states<- lapply(TS, \(x) x<-tip_states)					# reorder
  for (i in seq_along(tip_states)){
    for (k in seq_along(RHotTrees)){
      tip_states[[i]][[k]][TipsAboveT[[i]][[k]]]<-2		# Complete tip states
    }
  }
  ST1<-stats::setNames(vector("list", length(thr_level)), thr_level)
  st1_probs<- lapply(ST1, \(x) x<-st1_probs)					# reorder
  # Complete state 1 pp
  for (i in seq_along(st1_probs)){
    for (k in seq_along(RHotTrees)){
      st1_probs[[i]][[k]][names(probsAboveT[[i]][[k]]),1]<-probsAboveT[[i]][[k]]
    }
  }
  #Hidden state prediction & output preparation
  hspb<-lapply(seq_along(thr_level),\(y) lapply(seq_along(RHotTrees),
                                                \(x) castor::hsp_binomial(tree=RHotTrees[[x]],
                                                                          tip_states=tip_states[[y]][[x]],
                                                                          reveal_probs =NULL,
                                                                          state1_probs = st1_probs[[y]][[x]],
                                                                          min_revealed = min_revealed,
                                                                          max_STE = max_STE,
                                                                          check_input = TRUE)))
  p1s<-lapply(seq_along(RHotTrees), \(x) lapply(seq_along(hspb),
                                                \(y) hspb[[y]][[x]]$P1[1:ape::Ntip(RHotTrees[[x]])]))
  p1<-lapply(seq_along(p1s),
             \(x) round(apply(mapply("c", p1s[[x]]),1, mean), digits=5))
  stes<-lapply(seq_along(RHotTrees),\(x) lapply(seq_along(hspb),
                                                \(y) hspb[[y]][[x]]$STE[1:ape::Ntip(RHotTrees[[x]])]))
  ste<-lapply(seq_along(stes),
              \(x) round(apply(mapply("c", stes[[x]]),1, mean), digits=5))
  st1p<-lapply(seq_along(RHotTrees),
               \(x) lapply(seq_along(thr_level), \(y) st1_probs[[y]][[x]][,1]))
  st1_probs<-lapply(seq_along(RHotTrees),
                    \(x) apply(round(mapply("c", st1p[[x]]),digits=5),1, paste, collapse="|"))
  st1_p<-lapply(st1_probs, \(x) paste0("c(", x, ",0)"))
  ordtip_st <- lapply(seq_along(RHotTrees),\(x) lapply(seq_along(tip_states),
                                                       \(y) tip_states[[y]][[x]]))
  tip_states<-lapply(seq_along(RHotTrees),
                     \(x) apply(mapply("c", ordtip_st[[x]]),1,paste, collapse="|"))
  all4all_out<-list(p1=p1, ste=ste, tip_states=tip_states, st1_p=st1_p)
}

#'
#'@title Calculate the (relative) probability of taxa to be medicinal species
#'  using HSP analysis.
#'
#'@description This function computes the (relative) probability of being a
#'  medicinal plant (i.e., 'state 1') for all taxa included in the "hot trees"
#'  with the number of tips greater than a user defined 'cut_off'. The function
#'  starts combining the information on the known medicinal plants with the
#'  computed phylogenetic distances to identify potential non-medicinal plants
#'  (i.e., 'state 2') according with an established distance threshold. Then the
#'  hidden state prediction (HSP) technique (Zaneveld and Thurber, 2014) for a
#'  binary trait (i.e., medicinal vs non-medicinal plants) based on the binomial
#'  distribution is applied to estimate the state probabilities of hidden
#'  states.
#'
#'@details The core of this function uses the 'hsp_binomial()' function
#'  available in the 'castor' R  package to calculate the (relative) probability
#'  of each taxon of being a medicinal plants. Others functions from 'ape',
#'  'castor', 'cluster', 'factoextra','ggplot2', 'phytools' and 'stats' R
#'  packages are also used in the computation. To be able to carry out this type
#'  of analysis it is necessary to assign some taxa to both states of the
#'  analyzed character. To this aim, the 'state 1' is assigned to all known
#'  medicinal plants (i.e., all species provided by the 'sample_file') assuming
#'  no error in measurement (i.e. setting the "state1_probs" parameter to
#'  c(1,0), see below). Less obvious is the assignment of taxa to 'state 2' for
#'  which the function operates as follows. Within each selected "hot tree"
#'  (i.e., all "hot trees" with Ntips > cut_off), the function computes the
#'  phylogenetic distance between known medicinal taxa and the remaining taxa.
#'  The minimum distances calculated for each taxon not belonging to the
#'  medicinal plant group are pooled together into a "minimum distance" sample.
#'  Then one of the following two approaches can be applied to assign species to
#'  'state 2' base on user preference:
#' \itemize{
#'  \item The score corresponding to the percentile of the "minimum distance"
#' sample specified by the user is calculated and applied as a reference distance threshold.
#' All taxa with a phylogenetic distance to the nearest known medicinal taxon greater
#' than this reference threshold are considered as species with putative low medicinal
#' potential and are assigned to 'state 2'.
#' Alternatively, multiple reference distance thresholds can be applied
#' simultaneously (see 'method' parameter).
#' \item  The "minimum distance" sample is automatically partitioned by a clustering
#' algorithm into an optimal number of clusters. The cluster including the greatest
#' distances is taken as a reference and the species corresponding to these distances
#' are considered as species with putative low medicinal potential and assigned to 'state 2'.}\cr
#'
#' The 'state1_probs' parameter available in the 'hsp_binomial()' function is
#' then used to account for potential state-measurement errors. For each
#' putative low medicinal potential taxon, its distance from the nearest known
#' medicinal plant is used to weight the uncertainty in measuring “state 1”.
#' To this aim, the percentile returned by the empirical cumulative distribution
#' function built on the "minimum distance" sample, optionally decreased by a
#' correction factor, is used to express probability of measuring the 'state 1'
#' conditional upon its true state actually being 'state 1' (and conditional upon
#' its state having been "measured"; see the 'castor' package manual).
#' In other words,the closer the distance to the defined threshold (s), the greater
#' the probability for a taxon that its 'state 1' has been erroneously "measured"
#' (and therefore that the 'status 2' was erroneously assigned). The correction
#' factor, if used, allows to avoid attributing 'state 2' with certainty to the
#' taxa that belong to a "hot tree" identified by the function "nodesigl_harvesteR".
#' To include at least one taxon with an error-free 'state 2', an additional step
#' is performed. The ancestor node of each selected "hot tree" is identified
#' within the “reference tree” from which it was obtained. Among the tips
#' descendant from this ancestor node and not included in the "hot tree", an
#' outgroup taxon is chosen randomly with the constraint that it does not
#' belong to the set of known medicinal plants. The selected taxon and the "hot
#' tree" are merged into a new tree rooted at their ancestor node, and 'state
#' 2' is assigned to the newly added outgroup assuming error-free measurement.
#' Seeds can be set, thus making taxon selection repeatable. The states of the
#' species not assigned to 'state 1' or 'state 2' are treated by the function
#' as hidden states.
#'
#'@param reftree Character (file name). The reference phylogeny from winch "hot
#'  trees" were extracted. It must be a tree file in Newick format with
#'  correctly labelled species (see 'clean_and_match()' function). The output of
#'  the 'tree4nodesigl()' function must be used (i.e. the same file used with
#'  the 'nodesiglR()' function).
#'
#'@param sample_file Character (file name). A three columns, tab delimited file
#'  containing the output of the 'sample4nodesigl()' function. The same file
#'  used with the 'nodesiglR()' function must be used.
#'
#'@param IndepTrees Character (file name). A text file including the list of
#'  "independent hot trees" to process written in Newick format. Typically, this
#'  is the "IndependentHotTreesNewick" produced by the 'nodesigl_harvesteR()'
#'  function.
#'
#'@param method Character. This parameter specifies the method to be applied to
#'  identify putative low medicinal potential plants based on a phylogenetic distance
#'  threshold from the nearest observed medicinal taxon. Currently, options are:
#'  \itemize{
#'  \item \emph{"single_thr"}: character. To identify species with putative low
#'  medicinal potential, a single user-defined distance threshold (see Details)
#'  is applied to each “hot tree” as a reference. The same reference threshold can
#'  be applied to all the "hot trees" analysed or different individual thresholds
#'   can be applied to different "hot trees" (see the 'thr_level'parameter).
#'   Final probabilities are computed based on a single phylogenetic distance
#'   threshold for each tree.
#'  \item \emph{"multi_thr"}: character. To identify putative low medicinal potential taxa,
#'  multiple distance thresholds are applied to each “hot tree”, defined by the user
#'  via the “thr_level” and “by” parameters. The same set of thresholds is applied
#'  to all trees analysed. For each tree the final probabilities are obtained by
#'  averaging all values computed under different distance thresholds.
#'  \item \emph{"kmeans"}: character. The K-means method implemented in the 'stats'
#'  R package is used to automatically identify the putative low medicinal potential plants.
#'  Phylogenetic distance included in the "minimum distance" sample (see Details)
#'  are grouped by the algorithm into a number of clusters ranging from 1 to 'kmax'
#'  (see 'kmax' parameter). The optimal number of clusters is determined using the
#'  "average silhouette" method implemented in the 'factoextra' R package.
#'  The species corresponding to the distances clustered around the largest centroid
#'  are considered as species with putative low medicinal potential and assigned to
#'  'state 2'.
#'  \item \emph{"pam"}: character.The PAM method implemented in the 'cluster' R
#'  package is used to automatically identify the putative low medicinal potential plants.
#'  Phylogenetic distance included in the "minimum distance" sample (see Details)
#'  are grouped by the algorithm into a number of clusters ranging from 1 to
#'  'kmax' (see 'kmax' parameter). The optimal number of clusters is determined
#'  using the "average silhouette" method implemented in the 'factoextra' R package.
#'  The species corresponding to the distances clustered around the largest medoid
#'  are considered as putative low medicinal potential species and assigned to
#'  'state 2'.
#'  PAM is sometimes considered more robust than the K-means method, because it is
#'  less  sensitive to outliers (i.e., the presence of few very long branches).
#'  Under certain circumstances the two methods can provide different results.
#' }
#'
#'@param thr_level A single numeric value or numeric vector. Values must be
#'  between 0 and 1. Four cases are possible:\cr
#' - If 'method' = "single_thr" and 'thr_level' is a single numeric value,
#' this value is taken as the percentile of the "minimum distance" sample to be
#' used to calculate the reference distance threshold. The same reference threshold
#' is applied to all the "hot trees" analysed.\cr
#' - If 'method' = "single_thr" and 'thr_level' is numeric vector, different values
#'  can be provided as vector elements so that different reference distance
#'  thresholds can be applied to different "hot trees". In this case the number
#'  of elements present in the vector must be equal to the number of "hot trees"
#'  analyzed and their order corresponds to that of the trees.
#'  A single threshold is applied to each tree, which however may differ from tree to tree.
#' - If 'method' = "multi_thr", then 'thr_level must be a numeric vector with
#' 2 elements . Their values correspond to the percentiles of the "minimum distance"
#' sample used to identify the minimum and maximum reference distance thresholds
#' to be used in the calculation. The number and the position of the intermediate
#' reference thresholds (i.e. those within the minimum and the maximum thresholds)
#' are determined via the 'by' parameter (see below). All the specified reference
#'  thresholds are applied to each tree.
#' - When the 'method' parameter is set to "kmeans" or "pam", then the parameter
#' 'thr_level' is ignored.
#'
#'@param by Numeric. The same as in 'seq()' function. It determines the increase
#'  in percentile between two subsequent reference thresholds when the selected
#'  'method' is "multi_thr". It must be set so that the difference between the
#'  maximum and minimum threshold is a multiple of the 'by' parameter.
#'
#'@param cut_off Integer, defines the size that a tree must have to be included
#'  in the analysis. Only trees with the number of tips > cut_off are considered
#'  (default: cut_off = 250).
#'
#'@param kmax Integer, sets the maximum number of clusters explored by the
#'  clustering algorithm when the 'method' parameter is set to "kmeans" or
#'  "pam". Ignored otherwise (default: kmax = 10).
#'
#'@param min_revealed The same as in the 'hsp_binomial()' function. See the
#'  'castor' package manual. Default: min_revealed = 2.
#'
#'@param max_STE The same as in the 'hsp_binomial()' function. See the 'castor'
#'  package manual. Default: max_STE = 0.25.
#'
#'@param seed Numeric, specifies seeds for reproducible results (see Details;
#'  default: seed=123).
#'
#'@param cf, correction factor used to correct the prior probability of measuring the
#'    'state 1' applied to each tip(default =0.1).Ignored if mthds is "kmeans" or "pam".
#'
#'@return A csv file with nine columns is saved in the current working
#'  directory. The file is named according to the format:
#'  "hmpp_t\emph{i}_method_thr_level_disease.csv", where "t\emph{i}" is the
#'  "hot tree" analysed (with \emph{i} corresponding to the \emph{i}-th tree in
#'  the "IndepTrees" file), "method" and "thr_level" are the method and the
#'  threshold(s) level applied and "disease" is the name of the disease of interest
#'  taken directly from the name of the sample_file (i.e.,it is equal to the
#'  'disout' parameter of the 'sample4nodesigl()' function).
#'
#'  The column headings and related information displayed are as follows:\cr
#'
#'  \emph{- species}, names of the species included in the analysed "hot
#'  tree";\cr
#'  \emph{- known_med_plants}, logical values. The column indicates
#'  whether a plant is an already known medicinal plant;\cr
#'  \emph{- hsp_p1},for each tip in the "hot tree", the estimated (relative)
#'  probability of being a medicinal plant (i.e., of being 'state 1') is shown;\cr
#'  \emph{max_est_p1}, logical values. The column indicates those taxa for which
#'  the estimated probability is equal to the maximum estimated probability in
#'  the "hot tree".Since for known medicinal plants the probability is fixed and
#'  not estimated, these are not taken into account when assessing the maximum
#'  probability among the tips;\cr
#'  \emph{hsp_p1_ratio},for each species in the "hot tree" the ratio between
#'  the estimated probability for that species and the estimated maximum probability
#'  in the tree is shown;\cr
#'  \emph{- hsp_p2},for each tip in the "hot tree", the column shows the (relative)
#'  probability of \emph{not} being a medicinal plant (i.e., of being 'state 2').
#'  It is computed as \emph{hsp_p2} = 1 - \emph{hsp_p1};\cr
#'  \emph{- hsp_STE_p1},for each tip in the "hot tree", the column shows the
#'  standard error of the estimated hsp_p1. See the 'castor' package manual for details;\cr
#'  \emph{- prior_state},for each tip in the "hot tree", the column shows the
#'  \emph{a priori} state assumed by the analysis, based on the specified reference
#'  threshold. Possible values are: 1, 2 or NA for hidden states. When multiple
#'  reference thresholds are used, the different assumed states are shown
#'  separated by "|" symbol following the ascending order of the reference thresholds;\cr
#'  \emph{-prior_state1_probs},for each tip in the "hot tree", the column shows the
#'  \emph{a priori} probability of 'state 1' assumed by the analysis, based on
#'  the specified reference threshold(s). These values express the assumed prior
#'  probability of measuring the 'state 1' conditional upon its "true" state is
#'  indeed 'state 1'and conditional upon its state is "non-hidden" (see the
#'  'castor' package manual). When a single reference threshold is used (i.e.,
#'  when 'method' is set to "single_thr", "kmeans" or "pam") probabilities are
#'  shown in the form of a vector of length 2. Error-free measurement (i.e.,
#'  c(1,0)) is assumed for the known medicinal species and for the outgroup
#'  taxon (see Details). For convenience, the same is printed for taxa with
#'  prior hidden state (i.e., "NA") even if the information is not used in
#'  computation by the function. When multiple reference thresholds are used
#'  (i.e., when 'method' is set to "multi_thr") a modified format is used where
#'  the different assumed probabilities for 'state 1' are shown separated by "|"
#'  symbol following the ascending order of the reference thresholds. Since the
#'  assumed value for 'state 2' never changes, only one value is shown. For
#'  example, the modified format c(0.88792|0.88792|1,0) shows the assumed
#'  probabilities for a tip, based on three different reference thresholds).\cr
#'
#'  A second pdf file is created when the 'kmeans' or 'pam' methods are
#'  selected, showing the optimal number of clusters obtained through the
#'  "average silhouette" method. The file is named according to the
#'  format:t\emph{i}_method_k_disease.pdf, where "t\emph{i}" is the "hot tree"
#'  analysed (with \emph{i} corresponding to the \emph{i}-th tree in the
#'  "IndepTrees" file), "method" the method applied, "k" is the optimal number
#'  of clusters identified and "disease" is the name of the disease of interest
#'  taken directly from the name of the sample_file (i.e.,it is equal to the
#'  'disout' parameter of the 'sample4nodesigl()' function).
#'
#'@author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#'
#'@references
#' - Maechler M., Rousseeuw P., Struyf A., Hubert M., Hornik K. (2023).
#'   cluster: Cluster Analysis Basics and Extensions.
#'   https://CRAN.R-project.org/package=cluster.
#' - Louca S, Doebeli M (2017). “Efficient comparative phylogenetics on large trees.”
#'   Bioinformatics. doi:10.1093/bioinformatics/btx701.
#' - Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics
#'   and evolutionary analyses in R.” Bioinformatics, 35, 526-528.
#'   doi:10.1093/bioinformatics/bty633.
#' - Kassambara A, Mundt F. (2022)   factoextra: Extract and Visualize the Results
#'   of Multivariate Data Analyses. https://CRAN.R-project.org/package=factoextra
#' - Revell L.J.(2024). phytools 2.0: an updated R ecosystem for phylogenetic
#'   comparative methods (and other things). PeerJ 12:e16505
#'   https://doi.org/10.7717/peerj.16505
#' - Webb, C. O.; Ackerly, D. D. & Kembel, S. W. (2008) Phylocom: software for the
#'   analysis of phylogenetic community structure and trait evolution.
#'   Bioinformatics 24: 2098-2100.
#' - Zaneveld K.R., and Thurber  R. L. V. (2014). Hidden state prediction: A
#'   modification of classic ancestral state reconstruction algorithms helps unravel
#'   complex symbioses. Frontiers in Microbiology. 5:431.
#' - Zecca, G., Toini, E., Labra, M, Grassi, F. (2025) Accelerating the
#'   identification and the prioritisation of new plants with medicinal
#'   potential: the  pm4mp R package.\cr
#'
#'@seealso [clean_and_match()], [sample4nodesigl()], [tree4nodesigl()]
#'
#'@examplesIf interactive()
#' # WARNING: to run the example a copy of the Phylocom (Webb et al. 2008) binaries MUST be in the
#' # current temp directory. # DO NOT RUN otherwise!
#' # Use tempdir() to find the current temp directory.
#'
#'  library(ape)
#'  library(castor)
#'  library(factoextra)
#'  library(cluster)
#'  library(stats)
#'  library(phytools)
#'  library(ggplot2)
#'
#'  # Generate a random tree with 450 tips and save it in the current temp folder
#'  set.seed(1) # do not change it
#'  tr<-rtree(450)
#'
#'  # Add dummy genera
#'  tr$tip.label<-paste0(rep(paste0("Genus", LETTERS[1:18]),each=25),"_",tr$tip.label)
#'
#'  # Plot the tree
#'  tipcol<-c(rep("black",450))
#'  tipcol[c(10,11,12,13,14,25,26,27,28,29,30,32,35,41,42,43,46,47,48,55,56,57,
#'  58,59,60,101,102,103,105,106,108,109,110,309, 313, 318,319,320,325,350,351,
#'  352,353,354,355,356,357, 358,362,364,368,374,380,399,400,403,407,408,411)]<-"red"
#'  #known medicinal taxa in red
#'  plot(tr,tip.color =tipcol, cex=0.7, underscore=TRUE)
#'
#'  # Show medicinal taxa labels
#'  tr$tip.label[c(10,11,12,13,14,25,26,27,28,29,30,32,35,41,42,43,46,47,48,55,56,
#'  57,58,59,60,101,102,103,105,106,108,109,110,309, 313, 318,319,320,325,350,351,
#'  352,353,354,355,356,357, 358,362,364,368,374,380,399,400,403,407,408,411)]
#'
#'  # Write the tree in the current temp folder
#'  tree<-write.tree(tr)
#'  treefile <-tempfile("tree", fileext = ".tree")
#'  cat(tree, file = treefile, sep = "\n")
#'
#'  # Create a sample file and write it in the current temp folder
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
#'  # Get the temporary directory path and set it as current working directory
#'  path<-paste0(normalizePath(tempdir(), winslash = "/"),"/")
#'  setwd(path)
#'
#'  # Output files are in the current temp folder (type 'path')
#'  tree4nodesigl(treefile, "example")
#'  # it may take a few minutes to complete
#'  nodesiglR(start=1,stop=500, s=samplefile, f="phylo_nobl_example")
#'  nodesigl_harvesteR(sample_file=basename(samplefile),tree="phylo_nobl_example",fract=5)
#'  files<-list.files(path, pattern ="IndependentHot")
#'
#'  #Test. Output files are in the current temp folder
#'  # Single threshold set at: 0.70
#'  hmpp(reftree="phylo_nobl_example",sample_file=basename(samplefile),
#'                  IndepTrees=files[2],method="single_thr", thr_level=0.7,
#'                  cut_off=100,min_revealed = 2,max_STE=0.25,seed=123)
#'
#'  # Multiple (5) thresholds set at: 0.75, 0.80, 0.85, 0.90, 0.95
#'  hmpp(reftree="phylo_nobl_example",sample_file=basename(samplefile),
#'                 IndepTrees=files[2],method="multi_thr", thr_level=c(0.75,0.95),
#'                 by=0.05, cut_off=100,min_revealed= 2,max_STE =0.25,seed=123, cf=0)
#'
#'  # kmeans method testing a maximum of 5 clusters
#'  hmpp(reftree="phylo_nobl_example",sample_file=basename(samplefile),
#'                  IndepTrees=files[2], method="kmeans", cut_off=100,kmax=5,
#'                  min_revealed = 2,max_STE =0.25,seed=123)
#'
#'  # PAM method testing a maximum of 20 clusters
#'  hmpp(reftree="phylo_nobl_example",sample_file=basename(samplefile),
#'                  IndepTrees=files[2], method="pam", cut_off=100,kmax =20,
#'                  min_revealed = 2,max_STE =0.25,seed=123)
#'
#'  # Tidy up. Remove the output files from the current temporary directory.
#'  toremove<-list.files(path, pattern = ".tree|phylo_|disease|Nodesigl_|hmpp_")
#'  file.remove(paste0(path,toremove))
#'
#'  # You can delete the Phylocom copy yourself if you want
#'
#'@export
#'

hmpp<-function(reftree,sample_file,IndepTrees,method, thr_level,by=NULL,
               cut_off=100,kmax =10,min_revealed = 2,max_STE =0.25,seed=123, cf=0.1){
  #Input files & cut_off...
  set.seed(seed)
  Btree<-ape::read.tree(reftree)
  Medtips<-unlist(utils::read.delim(sample_file,header=FALSE,colClasses=c("NULL","NULL",NA)),
                  use.names = FALSE)
  HotTrees<-ape::read.tree(file=IndepTrees)
  if(inherits(HotTrees,"phylo")){
    HotTrees<-phytools::as.multiPhylo(HotTrees)}
  names(HotTrees)<-paste0("t",seq_along(HotTrees))
  cutoff<-sapply(HotTrees, \(x) ape::Ntip(x)>cut_off)
  HotTrees<-HotTrees[cutoff]
  #Checks...
  if(length(HotTrees)==0){
    stop("No tree has passed the set 'cut_off' threshold.")
  }
  if(!is.null(by) && method!="multi_thr" || is.null(by)&& method=="multi_thr"){
    stop("To apply the multi-threshold method
         set both the 'method' and 'by'parameters appropriately")
  }else if(!is.null(by)&& method=="multi_thr"){
    if(round((thr_level[2]-thr_level[1])/ by, digits=5)%%1 ==0){
      thr_level<-seq(thr_level[1],thr_level[2], by=by)
    }else{
      stop("To apply the multi-threshold method,the difference between the maximum
           and minimum threshold must be a multiple of the 'by' parameter")
    }
  }
  #Getting new trees with a random non-med outgroup...
  to_root_out<-to_root(HotTrees,Btree,Medtips)
  RHotTrees<-to_root_out[["RHotTrees"]]
  outgrnames<-to_root_out[["outgrnames"]]
  # Identify a reference threshold (4 methods) from minimal distances background,
  # find the tips above this threshold and calculete their PP...
  min_dist_out<-min_dist(HotTrees,Btree,Medtips)
  TreesHotTipsList<-min_dist_out[["TreesHotTipsList"]]
  TreesColdTipsList<-min_dist_out[["TreesColdTipsList"]]
  MinPhyDistList<-min_dist_out[["MinPhyDistList"]]
  #PP
  EmCDFs<-lapply(MinPhyDistList, \(x) stats::ecdf(x))
  #Find tips above threshold(s)
  mthds_out<-mthds(method,thr_level,kmax,MinPhyDistList,TreesColdTipsList)
  TipsAboveT<-mthds_out[["TipsAboveT"]]
  if (length(mthds_out)>1){
    bestk<- mthds_out[["bestk"]]
    ks<- mthds_out[["ks"]]}
  #Prepare fixed tips states...
  tip_states<-lapply(RHotTrees, \(x) rep(NA, ape::Ntip(x)))
  # ...and prepare 'state 1' PP
  st1_probs<-lapply(RHotTrees, \(x) matrix(rep(c(1,0), ape::Ntip(x)),ncol=2,
                                           nrow= ape::Ntip(x), byrow=TRUE))
  for(i in seq_along(RHotTrees)){
    names(tip_states[[i]]) <- RHotTrees[[i]]$tip.label
    tip_states[[i]][TreesHotTipsList[[i]]]<-1
    tip_states[[i]][outgrnames[[i]]]<- 2
    rownames(st1_probs[[i]])<- RHotTrees[[i]]$tip.label
  }
  if (method != "multi_thr"){
    one4all_out<-one4all(MinPhyDistList, TipsAboveT, EmCDFs, RHotTrees, tip_states,
                         st1_probs, min_revealed, max_STE,cf)
    p1 <-one4all_out[["p1"]]
    ste<-one4all_out[["ste"]]
    tip_states<- one4all_out[["tip_states"]]
    st1_p<-one4all_out[["st1_p"]]
  }else{
    all4all_out<-all4all(TipsAboveT, MinPhyDistList, EmCDFs, RHotTrees, thr_level,
                         tip_states, st1_probs, min_revealed, max_STE,cf)
    p1 <-all4all_out[["p1"]]
    ste<-all4all_out[["ste"]]
    tip_states<-all4all_out[["tip_states"]]
    st1_p<-all4all_out[["st1_p"]]
  }
  #Output
  tips<-lapply(RHotTrees, \(x) x$tip.label)
  outNumb<-mapply(\(x,y) which(y %in% x), x=outgrnames, y=tips, SIMPLIFY=FALSE)
  for(i in seq_along(tips)){
    tips[[i]][outNumb[[i]]]<-paste0("OUTGROUP_",outgrnames[[i]])
  }
  p2<-lapply(p1, \(x) 1-x)
  kmp<-lapply(tips, \(x) x %in% Medtips)
  p1max<-lapply(seq_along(p1), \(x) max(p1[[x]][!kmp[[x]]]))
  maxest<-mapply(\(x,y) y==x, y=p1, x=p1max, SIMPLIFY=FALSE)
  p1_r<-mapply(\(x,y) y/x, y=p1, x=p1max, SIMPLIFY=FALSE)
  p1_r<-lapply(p1_r, \(x) replace(x, x>1, 1))
  df_p<- mapply(cbind, species=tips,known_med_plants=kmp, hsp_p1=p1,
                max_est_p1=maxest,hsp_p1_ratio=p1_r, hsp_p2=p2, hsp_STE_p1=ste,
                prior_state=tip_states,prior_state1_probs=st1_p,SIMPLIFY=FALSE)
  disease<- sub("sample_", "", basename(sample_file))
  if(method == "kmeans" || method=="pam"){
    invisible(lapply(seq_along(df_p),
                     \(x) utils::write.csv(df_p[[x]],file=paste0("hmpp_",names(HotTrees)[x],"_",
                                                                 method,"_k",bestk[x],"_", disease,".csv"), row.names = FALSE)))
    names(ks)<-names(HotTrees)
    for(i in seq_along(ks)){
      print(ks[i])
      ggplot2::ggsave(filename=paste0(names(HotTrees)[i],"_", method,"_k",bestk[i],
                                      "_", disease,".pdf") , device ="pdf", dpi="retina")
    }
    grDevices::dev.off()
  }else{
    if (min(thr_level)!= max(thr_level)){thr_level<-paste0(min(thr_level),"_",max(thr_level),"_by",by)}
    invisible(lapply(seq_along(df_p),
                     \(x) utils::write.csv(df_p[[x]], file=paste0("hmpp_",names(HotTrees)[x],
                                                                  "_", method, "_", thr_level, "_", disease, ".csv"),row.names = FALSE)))
  }
}

