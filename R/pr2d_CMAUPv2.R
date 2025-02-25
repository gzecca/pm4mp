#'
#' @title get_stuff2()
#' @usage NULL
#' @description Internal function for collecting information used by the
#'   'pr2d_CMAUPv2' function.
#' @details This function retrieves names, ICD-11 codes and URLs of all diseases
#'   archived on the CMAUP v2.0 website.
#' @return A list of R objects.
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#' @examplesIf interactive()
#' library(rvest)
#' library(polite)
#' # Scraping data from CMAUP v2.0
#' get_stuff2()
#' #Remove stuff
#' rm(stuff)
#' @import rvest
#' @keywords internal
#' @noRd
#'

get_stuff2<-function(){
  bow<-polite::bow("https://bidd.group/CMAUP/browseplantbydis.php")
  html<-polite::scrape(bow)
  mytab<- html %>% rvest::html_node("table")
  links<- mytab %>% rvest::html_elements("a") %>% rvest::html_attr ("href")
  dnames<- mytab %>% rvest::html_elements("a")%>% rvest::html_text2()
  dnames<-gsub("/|,", "_",dnames)
  dcod<- mytab %>%rvest::html_elements("tr")%>% rvest::html_element("td")%>%
    rvest::html_text2()
  dcodes<-sapply(strsplit(dcod[-1], " "), \(x) paste0("ICD-11-",x[1]))
  stuff<-list(bow=bow,links=links, dnames=dnames, dcodes=dcodes)
}


#' @title Scraping Plants Related to Diseases (CMAUP v2.0)
#'
#' @description The lists of plants associated with the requested diseases are
#'   retrieved from CMAUP v2.0 site.
#'
#' @details This function creates several files storing lists of plants
#'   associated with the required diseases within the current working directory.
#'   Please note that at least one among the parameters ALL, MED and DRUG must
#'   be set as TRUE.
#'
#' @param ALL Logical (default=TRUE), indicating whether files including all
#'   plants related to the required diseases must be created. At least one of
#'   the ALL and MED parameters must be set as TRUE.
#'
#' @param MED Logical (default=TRUE), indicating whether files including only
#'   the "Medicinal Plants" (as defined by the CMAUP's Plant Usage Classes) must
#'   be created for the requested diseases.
#'
#' @param DRUG Logical (default=TRUE), indicating whether files including only
#'   the "Drug-Producing Plants" (as defined by the CMAUP's Plant Usage Classes)
#'   must be created for the requested diseases.
#'
#' @param dn Integer vector, specifying the diseases of interest. Diseases are
#'   selected using the internal numbering used by CMAUP v2.0. The number
#'   corresponding to each disease can be easily obtained using the
#'   'dnu_CMAUPv2' function. Please note that disease names and numbers differ
#'   between CMAUP v1.0 and CMAUP v2.0. When dn= NULL (default), information
#'   relating to all the diseases present on the CMAUP v2.0 website are included
#'   in downloads, according to the values of the ALL, MED and DRUG parameters.
#'   Keep in mind that in this case the function may take a long time to
#'   complete its work (polite webscraping is applied and hundreds of file are
#'   created in the current working directory). See examples below.
#'
#' @param outpath path for the output files. When outpath=NULL (default) the
#'   current working directory is used.
#'
#' @return For each disease a three-columns csv file is created containing
#'   information relating to the name, genus and family of plants linked to that
#'   disease. Values are separate by comma. Files are named according to the
#'   following scheme: "Disease name_ICD-11 code_ALL/MED/DRUG_number of plants
#'   included in the file_CMAUP2.csv"
#'
#' @author Giovanni Zecca, \email{giovanni.zecca@@unimib.it}
#'
#' @references
#' - Dongyue Hou, Hanbo Lin, Yuhan Feng, et al. CMAUP database update 2024: extended
#'   functional and association information of of useful plants for biomedical research.
#'   Nucleic Acids Research 2024; DOI: doi.org/10.1093/nar/gky965.\cr
#' - Perepolkin D. (2023) polite: Be Nice on the Web. https://github.com/dmi3kno/polite \cr
#' - Wickham H (2023). rvest: Easily Harvest (Scrape) Web Pages.
#'   https://rvest.tidyverse.org/, https://github.com/tidyverse/rvest.
#' - Zecca, G., Toini, E., Labra, M, Grassi, F. (2025) Accelerating the
#'   identification and the prioritisation of new plants with medicinal
#'   potential: the  pm4mp R package.\cr
#'
#' @seealso [dnu_CMAUPv2]
#'
#' @examplesIf interactive()
#' library(rvest)
#' library(polite)
#'
#'  # NOT RUN
#'  # All diseases present on CMAUP v2.0 are selected. If n is the number of
#'  # diseases, 3n file are  created in the current working directory (i.e.,
#'  # n_ALL + n_MED + n_DRUG).
#'  # It may take a long time to complete!
#'
#'  # pr2d_CMAUPv2()
#'
#'  # All diseases (n) present on CMAUP v2.0 are selected, but only the n_MED +
#'  # n_DRUG files # are created. It may take a long time to complete!
#'
#'  # pr2d_CMAUPv2(ALL=FALSE)
#'
#'  # All diseases (n) present on CMAUP v2.0 are selected, but only the n_ALL
#'  # files are created. It may take a long time to complete!
#'
#'  # pr2d_CMAUPv2(MED=FALSE, DRUG=FALSE)
#'
#'  #RUN
#'  # Use a temporary directory.
#'  path<-paste0(normalizePath(tempdir(), winslash = "/"),"/")
#'
#'  # Only disease#1 is selected and only one file is created in the temporary folder.
#'
#'  pr2d_CMAUPv2(dn=c(1), MED=FALSE, DRUG=FALSE, outpath=path)
#'
#'  # Only disease 2,7 and 25 are selected and only six files are created (ie.,
#'  # 3_ALL + 3_DRUG) in the temporary folder.
#'
#'  pr2d_CMAUPv2(dn=c(2,7,25), MED=FALSE, outpath=path)
#'
#'  #Remove the output file from the current temporary directory.
#'  toremove<-list.files(path, pattern = "*_CMAUP2.csv$")
#'  file.remove(paste0(path,toremove))
#'
#' @import rvest
#' @export
#'

pr2d_CMAUPv2<-function(ALL=TRUE, MED=TRUE, DRUG=TRUE,dn=NULL, outpath=NULL){
  if (is.null(outpath)){
    outpath<-paste0(getwd(),"/")
  }else{
    outpath<-outpath
  }
  stuff<-get_stuff2()
  if (is.null(dn)){
    query<-stuff$links
  }else{
    query<-paste0("searchresults.php?disease=disease",dn)
  }
  for (i in which(stuff$links %in% query)){
    name<-stuff$dnames[i]
    code<-stuff$dcodes[i]
    dpage<- polite::nod(bow=stuff$bow, path=paste0("CMAUP/",stuff$links[i]))%>%
      polite::scrape()
    tab_all<-dpage%>% rvest::html_node("table")%>% rvest::html_table()
    if(isTRUE(ALL)){
      n_all<-as.character(length(tab_all$"Plant Name"))
      tab_all %>% utils::write.table(file= paste0(outpath,name,"_",code,"_ALL_",n_all,"_CMAUP2.csv"),
                                     sep=",", row.names=FALSE)
    }
    if(isTRUE(MED)){
      imgm<-dpage%>% rvest::html_elements('tr') %>%
        rvest::html_element('img[src="img/icon_medHerb.png"]')
      med<-!is.na(imgm[-1])
      tab_med<-tab_all[med,]
      n_med<-as.character(length(tab_med$"Plant Name"))
      tab_med %>% utils::write.table(file= paste0(outpath, name,"_",code,"_MED_",n_med,"_CMAUP2.csv"),
                                     sep=",", row.names=FALSE)
    }
    if(isTRUE(DRUG)){
      imgd<-dpage%>% rvest::html_elements('tr') %>%
        rvest::html_element('img[src="img/icon_drug_producing_plant.png"]')
      drug<-!is.na(imgd[-1])
      tab_drug<-tab_all[drug,]
      n_drug<-as.character(length(tab_drug$"Plant Name"))
      tab_drug %>% utils::write.table(file= paste0(outpath, name,"_",code,"_DRUG_",n_drug,"_CMAUP2.csv"),
                                      sep=",", row.names=FALSE)
    }

  }
  print("Job done.", quote=FALSE)
}

