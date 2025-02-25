#'
#' @title get_stuff1()
#' @usage NULL
#' @description Internal function for collecting information used by the
#'   'pr2d_CMAUPv1' function.
#' @details This function retrieves names, ICD-10 codes and URLs of all diseases
#'   archived on the CMAUP v1.0 website.
#' @return A list of R objects.
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#' @examplesIf interactive()
#' #NOT RUN
#' library(rvest)
#' library(polite)
#' # Scraping data from CMAUP v1.0
#' get_stuff1()
#' # Remove stuff
#' rm(stuff)
#' @import rvest
#' @keywords internal
#' @noRd
#'

get_stuff1<-function(){
  bow<-polite::bow("https://bidd.group/CMAUP-2019/browseplantbydis.php")
  html<-polite::scrape(bow)
  mytab<- html %>% rvest::html_node("table")
  links<- mytab %>% rvest::html_elements("a") %>% rvest::html_attr ("href")
  dnames<- mytab %>% rvest::html_elements("a")%>%rvest::html_text2()
  dnames<-gsub("/", "_",dnames)
  dcod<- mytab %>%rvest::html_elements("tr")%>%rvest::html_element("td")%>%
    rvest::html_text2()
  dcodes<-sapply(strsplit(dcod[-1], ";"), \(x) paste0("ICD-10-", x[1]))
  stuff<-list(bow=bow,links=links, dnames=dnames, dcodes=dcodes)
}



#' @title Scraping Plants Related to Diseases (CMAUP v1.0)
#'
#' @description The lists of plants associated with the requested diseases are
#'   retrieved from CMAUP v1.0 site.
#'
#' @details This function creates several files storing lists of plants
#'   associated with the required diseases. Please note that at least one of the
#'   parameters ALL and MED must be set as TRUE.
#'
#' @param ALL Logical (default=TRUE), indicating whether files including all
#'   plants related to the required diseases must be created.
#'
#' @param MED Logical (default=TRUE), indicating whether files including only
#'   the "Medicinal Plants" (as defined by the CMAUP's Plant Usage Classes) must
#'   be created for the requested diseases.
#'
#' @param dn Integer vector, specifying the diseases of interest. Diseases are
#'   selected using the internal numbering used by CMAUP v1.0. The number
#'   corresponding to each disease can be easily obtained using the
#'   'dnu_CMAUPv1' function. Please note that disease names and numbers differ
#'   between CMAUP v1.0 and CMAUP v2.0. When dn= NULL (default), information
#'   relating to all the diseases present on the CMAUP v.10 website are included
#'   in downloads, according to the values of the ALL and MED parameters. Keep
#'   in mind that in this case the function may take a long time to complete its
#'   work (polite webscraping is applied and hundreds of file are created). See
#'   examples below.
#'
#' @param outpath path for the output files. When outpath=NULL (default) the
#'   current working directory is used.
#'
#' @return For each disease a three-columns csv file is created containing
#'   information relating to the name, genus and family of plants linked to that
#'   disease. Values are separate by comma. Files are named according to the
#'   following scheme: "Disease name_ICD-10 code_ALL/MED_number of plants
#'   included in the file_CMAUP1.csv"
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
#'   47(D1): D1118-D1127  DOI: doi.org/10.1093/nar/gky965.
#' - Zecca, G., Toini, E., Labra, M, Grassi, F. (2025) Accelerating the
#'   identification and the prioritisation of new plants with medicinal
#'   potential: the  pm4mp R package.\cr
#'
#' @seealso [dnu_CMAUPv1]
#'
#' @examplesIf interactive()
#' library(rvest)
#' library(polite)
#'
#'  # NOT RUN!
#'  # All diseases present on CMAUP v1.0 are selected. If n is the number of
#'  # diseases, 2n file are created (i.e., n_ALL + n_MED). It may take a long time
#'  # to complete!
#'
#'  # pr2d_CMAUPv1()
#'
#'  # All diseases (n) present on CMAUP v1.0 are selected, but only the n_MED
#'  # files are created. It may take a long time to complete!
#'
#'  # pr2d_CMAUPv1(ALL=FALSE)
#'
#'  # All diseases (n) present on CMAUP v1.0 are selected, but only the n_ALL
#'  # files are created. It may take a long time to complete!
#'
#'  # pr2d_CMAUPv1(MED=FALSE)
#'
#'  # RUN
#'  # Use tempdir() to find the current temp directory.
#'
#'  path<-paste0(normalizePath(tempdir(), winslash = "/"),"/")
#'
#'  # Only disease#1 is selected and only one file is created in the temporary folder.
#'
#'  pr2d_CMAUPv1(dn=c(1), MED=FALSE, outpath=path)
#'
#'  # Only disease# 2,7 and 25 are selected and only 3_ALL files are
#'  # created in the temporary folder.
#'
#'  pr2d_CMAUPv1(dn=c(2,7,25), MED=FALSE, outpath=path)
#'
#'  #Remove the output file from the current temporary directory.
#'  toremove<-list.files(path, pattern = "*_CMAUP1.csv$")
#'  file.remove(paste0(path,toremove))
#'
#' @import rvest
#' @export
#'

pr2d_CMAUPv1<-function(ALL=TRUE, MED=TRUE,dn=NULL, outpath=NULL){
  if (is.null(outpath)){
    outpath<-paste0(getwd(),"/")
  }else{
    outpath<-outpath
  }
  stuff<-get_stuff1()
  if (is.null(dn)){
    query<-stuff$links
  }else{
    query<-paste0("searchresults.php?disease=Disease",dn)
  }
  for (i in which(stuff$links %in% query)){
    name<-stuff$dnames[i]
    code<-stuff$dcodes[i]
    dpage<- polite::nod(bow=stuff$bow, path=paste0("CMAUP-2019/",stuff$link[i]))%>%
      polite::scrape()
    tab_all<-dpage%>% rvest::html_node("table")%>% rvest::html_table()
    if(isTRUE(ALL)){
      n_all<-as.character(length(tab_all$"Plant Name"))
      tab_all[c("Plant Name","Plant Genus","Plant Family")] %>%
        utils::write.table(file= paste0(outpath,name,"_",code,"_ALL_",n_all,"_CMAUP1.csv"),
                           sep=",", row.names=FALSE)
    }
    if(isTRUE(MED)){
      img<-dpage%>% rvest::html_elements('tr') %>%
        rvest::html_element('img[src="img/icon_medHerb.png"]')
      med<-!is.na(img[-1])
      tab_med<-tab_all[med,]
      n_med<-as.character(length(tab_med$"Plant Name"))
      tab_med[c("Plant Name","Plant Genus","Plant Family")] %>%
        utils::write.table(file= paste0(outpath,name,"_",code,"_MED_",n_med,"_CMAUP1.csv"),
                           sep=",", row.names=FALSE)
    }
  }
  print("Job done.", quote=FALSE)
}


