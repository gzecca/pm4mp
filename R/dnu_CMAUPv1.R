#'
#' @title Get Diseases, Numbers and URLs from CMAUP v1.0
#'
#' @description Scraping disease names, disease numbers and disease URLs from
#'   CMAUP v1.0.
#'
#' @details This function retrieve names and URLs of all the diseases included
#'   on the CMAUP v1.0 website. Since diseases are archived following an
#'   internal numbering (please note that disease names and  numbering are
#'   different between CMAUP v1.0 and CMAUP v2.0),the number corresponding to
#'   each disease is also provided. These numbers are useful to set the "dn"
#'   parameter in the 'pr2d_CMAUPv1' function.
#'
#' @param outpath path for the output files. When outpath=NULL (default) the
#'   current working directory is used.
#'
#' @return A text file ("Diseases_Numbers_URLs_v1") with three columns named
#'   "Disease name", "Disease_number" and "URL", respectively. Values are
#'   separate by comma and space.
#'
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#'
#' @references
#' - Perepolkin D. (2023) polite: Be Nice on the Web.
#'   https://github.com/dmi3kno/polite \cr
#' - Wickham H (2023). rvest: Easily Harvest (Scrape) Web Pages.
#'   https://rvest.tidyverse.org/, https://github.com/tidyverse/rvest. \cr
#' - Xian Zeng, Peng Zhang, Yali Wang, et al. CMAUP: a database of collective
#'   molecular activities of useful plants. Nucleic Acids Research 2019;
#'   47(D1): D1118-D1127  DOI:doi.org/10.1093/nar/gky965.
#' - Zecca, G., Toini, E., Labra, M, Grassi, F. (2025) Accelerating the
#'   identification and the prioritisation of new plants with medicinal
#'   potential: the  pm4mp R package.\cr
#'
#' @seealso [pr2d_CMAUPv1()], [dnu_CMAUPv2()]
#'
#' @examplesIf interactive()
#' # Use tempdir() to find the current temp directory.
#' library(rvest)
#' library(polite)
#' # Use a temporary directory.
#' path<-paste0(normalizePath(tempdir(), winslash = "/"),"/")
#' # Scraping data from CMAUP v1.0
#' dnu_CMAUPv1(outpath=path)
#' #Remove the output file from the current temporary directory.
#' toremove<-"Diseases_Numbers_URLs_v1"
#' file.remove(paste0(path,toremove))
#'
#'
#' @import rvest
#'
#' @export
#'

dnu_CMAUPv1<-function(outpath=NULL){
  if (is.null(outpath)){
    outpath<-paste0(getwd(),"/")
  }else{
    outpath<-outpath
  }
  bow<-polite::bow("https://bidd.group/CMAUP-2019/browseplantbydis.php")
  html<-polite::scrape(bow)
  htab<- html %>% rvest::html_node("table")
  tab<-htab %>% rvest::html_table()
  links<- htab %>% rvest::html_elements("a") %>% rvest::html_attr ("href")
  httplinks<-paste0("https://bidd.group/CMAUP-2019/",links)
  cbind(Disease_name=tab$Disease,
        Disease_number=gsub("searchresults.php?disease=Disease","",links, fixed=TRUE),
        URL=httplinks)%>%
  utils::write.table(file =paste0(outpath,"Diseases_Numbers_URLs_v1"),
                     sep=", ", quote=FALSE, row.names=FALSE)
}
