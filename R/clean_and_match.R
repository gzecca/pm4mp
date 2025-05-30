#'
#' @title match_n_qmatch()
#' @usage NULL
#' @description Internal function to find exact matches and quasi-matches (up to
#'   a user-set acceptance threshold).
#' @details This function compares the names of the reference species with the
#'   names of the medicinal plant identifying their exact correspondences. Fuzzy
#'   matches can be optionally detected to account for typos and misnamed taxa.
#' @param disSpecies,refSpecies,match_thr See the 'Clean_and_match' function
#'   documentation for details.
#' @return A list of R objects.
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#' @examplesIf interactive()
#' #NOT RUN
#'   library(stringdist)
#'   out_tabs<-match_n_qmatch(disSpecies=disSpecies, refSpecies=refSpecies,
#'   match_thr=match_thr)
#'   rm(out_tabs)
#' @keywords internal
#' @noRd
#'


match_n_qmatch<-function(disSpecies, refSpecies, match_thr){
  matching_tab<-as.data.frame(intersect(disSpecies, refSpecies))
  names(matching_tab)<-"Plant.Name"
  nomatching<-setdiff(disSpecies, refSpecies)

  if(isTRUE(length(nomatching)>0)){
    refDiff<-setdiff(refSpecies, matching_tab$Plant.Name) #avoid double matching...
    distlist<-lapply(nomatching, \(x) stringdist::stringsim(x,refDiff, method = "dl"))
    maxsim<-sapply(distlist, \(x) max(x))
    maxpos<-sapply(distlist, \(x) which.max(x))
    above_thr<- maxsim >= match_thr
    smatch_tab<-as.data.frame(cbind(Reference = refDiff[maxpos],
                                  DataBase = nomatching,
                                  Similarity = maxsim ))
    qmatch_tab<-smatch_tab[above_thr,]
    if(nrow(qmatch_tab)==0){
      qmatch_tab<-c("\n\nNone of the similarity scores calculated between potentially mismatched name pairs
                    \nwere equal to or higher than the chosen acceptance threshold (i.e., the value set
                    \nfor the match_thr parameter). Should you wish to save these potential fuzzy matches
                    \nto a file, please consider lowering the value of the match_thr parameter.\n\n\n")
    }
  }else{
    qmatch_tab<-c("All the medicinal plants provided by the user are listed in the reference list.\n\n\n")
  }
  exact100 <-length(matching_tab$Plant.Name)/length(disSpecies)*100
  cat(paste0("\nMedicinal plants exactly matching to the reference:\n n matches = ",
      length(matching_tab$Plant.Name), "\n % matches = ",round(exact100,3), "\n\n"))
  out_tabs<-list(exact_match = matching_tab, quasi_match=qmatch_tab)
}

#'
#' @title short_form()
#' @usage NULL
#' @description Internal function to change hybrid species names to the short
#'   form "Genus_x".
#' @details This function modifies the names of the hybrid species of medicinal
#'   plants and possibly of the reference species in the abbreviated form
#'   "Genus_x".
#' @param refSpecies,clean_ref See the 'Clean_and_match' function documentation
#'   for details.
#' @param out_tabs A list of R objects produced by the 'match_n_qmatch'
#'   function.
#' @return A list of R objects (a modified version of out_tabs list).
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#' @examplesIf interactive()
#' #NOT RUN
#' out_tabs<-short_form(out_tabs=out_tabs,refSpecies=refSpecies,
#'   clean_ref=clean_ref)
#'   rm(out_tabs)
#' @keywords internal
#' @noRd
#'

short_form<-function(out_tabs,refSpecies, clean_ref){
  out_tabs$exact_match$Plant.Name<-sub("(.*?_x)_.*", "\\1",out_tabs$exact_match$Plant.Name)
  if (clean_ref){
    refSpecies<-sub("(.*?_x)_.*", "\\1",refSpecies)
    Reference.Species<-list(Reference.Species=refSpecies)
    out_tabs<- append(out_tabs, Reference.Species)
  }
  out_tabs
}

#'
#' @title check_duplicates()
#' @usage NULL
#' @description Internal function to check for species names duplicated.
#' @details This function checks for duplicate names in the reference list and
#'   medicinal plant list, possibly produced by the 'clean_ref' and 'shortHybr'
#'   options.
#' @param refSpecies, clean_ref, shortHybr See the 'Clean_and_match' function
#'   documentation for details.
#' @param out_tabs A list of R objects produced by the 'match_n_qmatch' and
#'   'short_form' functions.
#' @return Duplicate names, if any, are printed on the screen with a warning. If
#'   no duplicate is found, this is stated on the screen.
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#' @examplesIf interactive()
#' #NOT RUN
#'   check_duplicates(out_tabs=out_tabs,refSpecies=refSpecies,
#'   clean_ref=clean_ref, shortHybr=shortHybr)
#' @keywords internal
#' @noRd
#'

check_duplicates<-function(out_tabs,refSpecies, clean_ref, shortHybr){
  dupMed<-duplicated(out_tabs$exact_match$Plant.Name)
  if(sum(dupMed)>0){
    cat("\nWARNING:\n the shortHybr option and/or the cleaning steps resulted in duplicate names in the medicinal taxa list.\n\n",
        "Duplicates:\n", paste(unique(out_tabs$exact_match$Plant.Name[dupMed]),
                         collapse="\n "), "\n\n")
  }else{
    cat("\nNo duplicate names found in in the medicinal taxa list.\n\n")
  }
  if(clean_ref && shortHybr){
    dupRef<-duplicated(out_tabs$Reference.Species)
  }else if(clean_ref && !shortHybr|!clean_ref){
    dupRef<-duplicated(refSpecies)
  }
  if(sum(dupRef)>0) {
    cat("\nWARNING:\n the shortHybr and/or clean_ref options resulted in duplicate names in the reference list.\n\n",
        "Duplicates:\n", paste(unique(refSpecies[dupRef]), collapse="\n "), "\n\n")
  }else{
    cat("\nNo duplicate names found in in the reference list.\n\n")
  }
}

#'
#' @title Clean and match species names
#'
#' @description This function performs several steps of cleaning species names
#'   linked to a specific disease to standardize their format before attempting
#'   to match them to the user-provided species list for reference. It can also
#'   optionally standardize the names in the reference species list if requested
#'   by the user. Fuzzy matches can be optionally detected to account for typos
#'   and misnamed taxa. Unless otherwise specified, input files are assumed to
#'   be in the current working directory.\cr Though this function allows for
#'   some flexibility, it is impossible to predict all possible eventualities.
#'   It is therefore advisable to carefully examine the results before
#'   proceeding further with subsequent analyses.
#'
#' @details A preliminary cleaning is applied to species names downloaded from
#'   the plants-diseases databases following Zecca et all. 2025. Several
#'   cleaning steps are performed. In particular:
#' \itemize{
#' \item  leading and/or trailing whitespace are removed from species names;
#' \item  Internal whitespaces are replaced by underscores;
#' \item  any multiple underscores are replaced by a single underscore;
#' \item  information relating subspecies (including hybrid subspecies),
#'         variants, cultivars etc. are discarded;
#' \item  all species names are written in the case sensitive format 'Genus_species',
#'         except for hybrids species;
#' \item  hybrids names are written in the format 'Genus_x_species' or,optionally,
#'        'Genus_x';
#' \item  ambiguous names, such as those named "Genus_sp/_spp/_st/_sp./_spp./_st.",
#'        are discarded from the list of plants related to the disease of interest.
#' }
#'   If requested by the user, similar cleaning steps, except deletion of
#'   accessions, can be applied also to the reference species list. \cr
#'   After the cleaning phases the function searches for an exact match between the
#'   names in the reference species list and the names found in the list of
#'   plants related to the disease of interest. Fuzzy matches are also
#'   identified against a threshold set by the user (see parameters). \cr
#'   Due to software constraints or readability needs, the reference phylogeny
#'   may not include full species names. For this reason the list of reference
#'   species can be provided in different ways by the user (see parameters).
#'
#' @param reference The list of reference species. It can be provided via one of
#'   the following options:
#' \itemize{
#' \item a "path_tree_file_name" relating to a tree in Newick format whose tips
#'      are the reference species;
#' \item a character vector containing the reference species names;
#' \item a "path_alignment_file_name" relating to an alignment file in sequential
#'        PHYLIP format from which the reference phylogeny was constructed
#'       (provided the species are the same and full names are included);
#' \item a character vector of length 2 with named elements specifying the alignment
#'       and the phylogeny of reference. The alignment and the tree used must be
#'       in in sequential PHYLIP and a Newick format, respectively.
#'       This option was expressly designed to work with the data of Zanne et al. (2014),
#'       in which case the 'reference' parameter must be:
#'       c(alignment="Vascular_Plants_32223_taxon.phy", tree="Vascular_Plants_rooted.dated.tre").
#'       In this case the alignment is used exclusively to reconstruct the complete
#'       names of the hybrid species while the tree is used for matching
#'       taxa after hybrids names correction.
#' }
#'
#' @param reftype Character, specifying the format of the 'reference' parameter.
#'   It can be one of:
#' \itemize{
#' \item "align": for an alignment in sequential PHYLIP format;
#' \item "tr": for a tree in Newick format;
#' \item "vect": for vector specified by the user;
#' \item "z": for a vector with the following form:
#'     c(alignment="path_alignment_file_name", tree="path_tree_file_name"),
#'     where alignment and tree are in sequential PHYLIP and Newick format, respectively.
#' }
#'   No other forms are allowed. Please note that the  "z" option will
#'   automatically set clean_ref = FALSE and shortHybr = TRUE, overriding any
#'   different settings chosen by the user.
#'
#' @param clean_ref Logical, specifying whether the reference species list
#'   should be cleaned or not before matching step (default=FALSE). If TRUE
#'   cleaned names of reference species are save in a separate file.
#'
#' @param shortHybr Logical, specifying whether hybrid species names should be
#'   shortened in the format 'Genus_x' as in the
#'   "Vascular_Plants_rooted.dated.tre" by Zanne et al. (2014) or not (deafult =
#'   TRUE).
#'
#' @param check_dup Logical, specifying whether the function should check for
#'   duplicate names in the reference list and in the list of medicinal plants,
#'   possibly produced by the 'clean_ref' and 'shortHybr' options.
#'
#' @param match_thr Numeric between 0 and 1 (inclusive). This parameter defines
#'   the acceptance threshold for string similarity score in fuzzy text
#'   searches, being 1 = exact match and 0 = no match at all. Similarity score
#'   are calculated via the 'stringsim' function ('stringdist' package) using
#'   the Full Damerau-Levenshtein distance (see the 'stringdist' package for
#'   details). Default value = 0.9.
#'
#' @param disease_file A csv file including the name of taxa linked to a
#'   specific disease of interest obtained from a medicinal plant database such
#'   as those downloaded using the functions pr2d_CMAUPv1 and pr2d_CMAUPv2.
#'
#' @param outpath Path to the output files. When 'outpath' is set to NULL
#'   (default) the current working directory is used.
#'
#' @return The number and percentage of medicinal plants that produced an exact
#'   match with the reference taxa are printed on the screen. Duplicate names,
#'   if any, are also printed on the screen when check_dup=TRUE. Two csv files
#'   named "diseasename_cleanedup_exmatch.csv" and
#'   "diseasename_cleanedup_qumatch_threshold.csv" are saved. The first one is
#'   single column file including all medicinal species that produced an exact
#'   match with the reference taxa. The latter is a three-columns file
#'   containing the name found in the reference list (Reference column), the
#'   name found in the list of medicinal plants (Database column) and the
#'   similarity score for the pair (Similarity column). Only pairs with
#'   similarity score equal or above the acceptance threshold are reported. The
#'   final choice about these matches is left to the user. If clean_ref = TRUE a
#'   third file named "Reference_species_cleanedup.csv" is saved which includes
#'   the reference names as they resulted at the end of the preliminary cleaning
#'   and formatting steps. This file can be used to update reference names
#'   before going further with subsequent analyses.
#'
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}, with the
#'   contribution of Francesco Artusa,
#'   \email{francesco.artusa@studenti.unimi.it} and Elisa Toini,
#'   \email{e.toini@@campus.unimib.it}
#'
#' @references
#' - Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics
#'   and evolutionary analyses in R.”  Bioinformatics, 35, 526-528.
#'   doi:10.1093/bioinformatics/bty633. \cr
#' - van der Loo M.P.J.,The stringdist Package for Approximate String Matching.
#'   The R Journal (2014) 6:1, pages 111-122.\cr
#' - Zanne, A., Tank, D., Cornwell, W. et al. Three keys to the radiation of
#'   angiosperms into freezing environments.  Nature 506, 89–92 (2014).
#'   https://doi.org/10.1038/nature12872
#' - Zecca, G., Toini, E., Labra, M, Grassi, F. (2025) Accelerating the
#'   identification and the prioritisation of new plants with medicinal
#'   potential: the  pm4mp R package.\cr
#'
#' @seealso [pr2d_CMAUPv1()], [pr2d_CMAUPv2()]
#'
#' @examplesIf interactive()
#' # NOT RUN
#' # The following example show you the setting used by Zecca et al. (2025)
#' # to prepare the input files for subsequent analysis.
#' # Plants linked to Malaria included in the file
#' # "Malaria_ unspecified_ICD-11-01_ALL_718_CMAUP2.csv" were analysed using
#' # reference alignment and tree fron Zanne et al. (2014).
#' # Fuzzy matches are detected imposing an acceptance threshold of 0.9.
#' # 'clean_ref' and 'shortHybr' parameters are automatically set to FALSE and TRUE,
#' # respectively. The current working directory is used for output files.
#' # library(ape)
#' # library(stringdist)
#' # Clean_and_match(reference=c(alignment="Vascular_Plants_32223_taxon.phy",
#' #    tree="Vascular_Plants_rooted.dated.tre"),reftype="z",
#' #		disease_file="Malaria_ unspecified_ICD-11-01_ALL_718_CMAUP2.csv", match_thr=0.9)
#'
#'
#' #'  # Toy example. Use tempdir() to find the current temp directory.
#'  library(stringdist)
#'  library(ape)
#'
#'  # Create a temporary reference alignment file.
#'  align<-c("18 10
#'  Genus1__species1\t	AGTAAGGTTC
#'  Genus1_species2\t	AGTAAGGTTC
#'  Genus2_species2\t	AGTAAGCTCC
#'  Genus2_species3\t	AGTAAGGTTC
#'  Genus2_Species4\t	AGTTACGGTC
#'  Genus3_species1\t	AGTAAGGTTC
#'  Genus3_species2\t	AGTACGGTTC
#'  Genus3_specIes3\t	AGTAAGGATG
#'  Genus3_species4\t	AGTAAGGTTC
#'  Genus3_species5\t	ACTTACGTTC
#'  Genus4_species1\t	ATTACCGTTC
#'  Genus4_species2\t	AGTCAGGTTC
#'  Genus5_species1\t	AGTAAGCCTC
#'  Genus6_x_species1\t	AGTTTGGTTC
#'  Genus6__x_species2\t	AATAAGGCTC
#'  Genus7_species1\t	AGTAAAGTAC
#'  Genus7_species2\t	AGAAAGGTCC
#'  Genus7_species3\t	AGTCACGTTC")
#'  alignfile <- tempfile("align", fileext = ".phy")
#'  cat(align, file = alignfile, sep = "\n")
#'  ref<-unlist(read.table(alignfile, colClasses=c(NA,"NULL"), skip=1))
#'  names(ref)<-NULL
#'
#'  # Create a temporary reference tree file.
#'  t<-rtree(18)
#'  t$tip.label<-ref
#'  tree<-write.tree(t)
#'  treefile <- tempfile("tree", fileext = ".tree")
#'  cat(tree, file = treefile, sep = "\n")
#'  #plotreference tree
#'  plot(t, underscore=TRUE)
#'
#'  # Create a temporary medicinal plants csv file.
#'  med<-c("Plant.Name,Plant Genus,Plant Family
#'Genus1\t  x species1 \t var pulcherrima, 'Genus1,ABB
#'Genus1\t\tx  species3, Genus1,ABB
#'genus2  Species1 sp. rubra, Genus2,CDD
#'genus2\u0020  Species2, Genus2,CDD
#'genus2  SPEcies3, Genus2,CDD
#'gENus2 species4, Genus2,CDD
#'Genus3  species1, Genus3,DCE
#'Genus3\u0020\u0020\u0020 species2\u0020 ssp ilex, Genus3,DCE
#'Genus3  species31, Genus3,DCE
#'Genus3 sp, Genus3,DCE
#'Genus3 \u0020\u0020 spp, Genus3,DCE
#'Genus4  species1 cv glauca, Genus4,ABB
#'Genus4  species2, Genus4,ABB
#'genus5  x\t species1 st variosa, Genus5,DCE
#'Genus6\t x \u0020\u0020 species1, Genus6,FGF
#'Genus6 species2, Genus6,FGF
#'Genus7\t  species1 x speciesT, Genus7,CDD
#'Genus7 x species21, Genus7,CDD
#'Genus7 x \u0020\u0020\u0020 species321, Genus7,CDD")
#'  medfile <- tempfile("med", fileext = ".csv")
#'  cat(med, file = medfile, sep ="\n")
#'
#'  #Get the temporary directory path and set it as current working directory.
#'  path<-paste0(normalizePath(tempdir(), winslash = "/"),"/")
#'  setwd(path)
#'
#'  #Test.
#'  clean_and_match(reference=alignfile,reftype="align",	disease_file=medfile,
#'  clean_ref=FALSE,shortHybr=TRUE,check_dup=TRUE, match_thr=0.9)
#'  cat("Results have been saved in the temporary folder:\n",path, "\n\n")
#'
#'  # WARNING: this chunk of code partially overwrite previous results!
#'  clean_and_match(reference=treefile,reftype="tr",	disease_file=medfile,
#'  clean_ref=TRUE,shortHybr=TRUE,check_dup=TRUE,match_thr=0.8)
#'  cat("Results have been saved in the temporary folder:\n", path, "\n\n")
#'
#'  #Tidy up. Remove the output files from the current temporary directory.
#'  file.remove(alignfile, treefile, medfile)
#'  toremove<-list.files(path, pattern =
#'  "_cleanedup_exmatch.csv$|_cleanedup_qumatch_|Reference_species_cleanedup.csv")
#'  file.remove(paste0(path,toremove))
#'
#' @export
#'

clean_and_match<-function(reference,reftype=c("align","tr", "vect", "z"),
                          clean_ref=FALSE,shortHybr=TRUE,check_dup=FALSE,
                          disease_file, match_thr=0.9, outpath=NULL){
  if (match_thr>1||match_thr<0){
    stop("'match_thr' argument must be a number between 0 and 1 (inclusive).")}
  if(reftype=="z"){clean_ref<-FALSE;shortHybr<-TRUE}
  write_ref<-ifelse(clean_ref, TRUE,FALSE)
  switch(reftype,
    align = {refSpecies<-unlist(utils::read.table(reference, colClasses=c(NA,"NULL"),
                                                  skip=1))
          names(refSpecies)<-NULL},
    tr = {mytree<-ape::read.tree(reference)
          refSpecies<- mytree$tip.label},
    vect = {refSpecies <- reference},
    z = {mytips<-ape::read.tree(reference["tree"])$tip.label
          myalign<-unlist(utils::read.table(reference["alignment"],
                                            colClasses=c(NA,"NULL"), skip=1))
          shortnames<-sort(grep("_x$", mytips, value=TRUE))
          longnames<-sort(setdiff(grep("_x_", myalign, value=TRUE),
                                  grep("subsp",myalign, value=TRUE)))
          mytips[match(shortnames, mytips)]<-longnames
          refSpecies<- mytips},
    stop("Invalid 'reftype' value!
         'reftype' argument must be one of: 'align','tr','vect','z'."))
  disease<-utils::read.csv(disease_file)
  names(disease)[which(names(disease) == "Plants")]<-"Plant.Name"
  disSpecies<-disease$Plant.Name
  disSpecies<-trimws(disSpecies)
  disSpecies<-gsub("\\h+","_",disSpecies, perl=TRUE)
  disSpecies<-gsub("_{2,}","_",disSpecies)
  disSpecies<-gsub("(^.)","\\U\\1",tolower(disSpecies), perl=TRUE)
  todel<-grep("_sp$|_spp$|_st$|_sp\\.$|_spp\\.$|_st\\.$", disSpecies)
  if (length(todel) >0){
    disSpecies<-disSpecies[-todel]}
  disHybr<- grep("_x_", disSpecies, fixed=TRUE, value=TRUE)
  disSpecies[- which(disSpecies %in% disHybr)]<- sub("(.*?_.*?)_.*", "\\1",
                                                 disSpecies[- which(disSpecies %in% disHybr)])
  if(clean_ref){
    refSpecies<-trimws(refSpecies)
    refSpecies<-gsub("\\h+","_",refSpecies, perl=TRUE)
    refSpecies<-gsub("_{2,}","_",refSpecies)
    refSpecies<-gsub("(^.)","\\U\\1",tolower(refSpecies), perl=TRUE)
    refHybr<-setdiff(grep("(.*?_x_.*?)",refSpecies,value=TRUE),
                     grep(".*?_.*?_.*?_x_.*?",refSpecies,value=TRUE))
    refSpecies[- which(refSpecies %in% refHybr)]<- sub("(.*?_.*?)_.*", "\\1",
                                                   refSpecies[- which(refSpecies %in% refHybr)])
  }
  out_tabs<-match_n_qmatch(disSpecies=disSpecies, refSpecies=refSpecies, match_thr=match_thr)
  if (shortHybr){
    out_tabs<-short_form(out_tabs=out_tabs,refSpecies=refSpecies, clean_ref=clean_ref)
  }
  if(check_dup){
    check_duplicates(out_tabs=out_tabs,refSpecies=refSpecies, clean_ref=clean_ref,
                     shortHybr=shortHybr)
  }
  prefix<-gsub(".csv","", disease_file, fixed=TRUE)
  utils::write.csv(out_tabs$exact_match,
                   file= paste0(outpath,prefix,"_cleanedup_exmatch.csv"),
                   quote = FALSE,row.names =FALSE)
  if (is.data.frame(out_tabs$quasi_match)){
    utils::write.csv(out_tabs$quasi_match,
                   file= paste0(outpath,prefix,"_cleanedup_qumatch_",match_thr,".csv"),
                   quote = FALSE,row.names =FALSE)
  }else{
    cat(out_tabs$quasi_match)
  }
  if (write_ref){
    Ref.Species<-ifelse (shortHybr, as.data.frame(out_tabs$Reference.Species),
                         as.data.frame(refSpecies))
    names(Ref.Species)<-"Reference.Species"
    utils::write.csv(Ref.Species, file= paste0(outpath,"Reference_species_cleanedup.csv"),
                     quote = FALSE,row.names =FALSE)
  }
}

