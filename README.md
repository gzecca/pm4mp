
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pm4mp

<!-- badges: start -->
<!-- badges: end -->

The aim of the `pm4mp` package is to support phylogenetic-based
identification and prioritisation of species with high medicinal
potential. The new methods implemented in the package propose solutions
to some of the inherent limitations of using ‘hot nodes’ in
bioprospecting, extend the ‘hot node’ approach to the concept of ‘hot
trees’ and offer new ways of prioritising potential new medicinal
species. These new functionalities are presented and discussed in Zecca
et al. (2025).

The functions included allow to:

- download lists of medicinal plants from [CMAUP v1.0 &
  v2.0](https://bidd.group/CMAUP/) public database (Xian et al. 2019;
  Dongyue et al. 2024);
- check and match taxa names against a reference phylogeny provided by
  the user;
- prepare input files for subsequent analyses;
- call an external executable of the software
  [Phylocom](https://github.com/phylocom/phylocom) (Webb et al.2008) to
  perform multiple replicates of the ‘nodsigl’ analysis;
- import, process and summarise the outputs obtained from multiple
  replicates of the ‘nodsigl’ analysis according to predefined criteria
  and identify the ‘hot trees’ for the studied disease;
- compute the Hot Ancestry Score for all species present in the
  identified ‘hot trees’;
- compute the (relative) probability of being a medicinal species for
  all the taxa present in the provided (hot) tree(s);
- show the results graphically.

**NOTE:The [Phylocom](https://github.com/phylocom/phylocom) application,
which is required to perform the ‘nodesigl’ analysis, is not included in
this package and must be downloaded separately by the user from the
author’s site.**

However, not all functions of the package require the presence of
Phylocom to be used. In particular the following functions can be used
without Phylocom being present:

- `clean_and_match()`
- `dnu_CMAUPv1()`
- `dnu_CMAUPv2()`
- `hmpp()`
- `pr2d_CMAUPv1()`
- `pr2d_CMAUPv2()`
- `sample4nodesigl()`
- `tree4nodesigl()`
- `visual_comp()`

Conversely, the following function depends directly on Phylocom being
present for it to work:

- `nodesiglR()`

Finally, the following functions, although not directly dependent on the
presence of Phylocom, directly or indirectly use the results of its
‘nodesigl’ analysis:

- `has()`(indirect use)
- `hot_tree_painteR()`(indirect use)
- `nodesigle_harvesteR()`(direct use)

### Installation

You can install the last version of `pm4mp` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gzecca/pm4mp")
```

### Load & attach

You can load and attach the installed version of `pm4mp` with:

``` r
library(pm4mp)
```

### Getting help

You can use the `help()` function and the `?` help operator to access
the documentation pages for `pm4mp` functions.

``` r
# load and attach
library(pm4mp)

# access the documentation page of the function has()
help("has")

# access the documentation page of the function hmpp()
?hmpp
```

### About examples in the documetation pages

- Currently, to illustrate how the package works, the examples in the
  documentation pages are built from scratch, without relying on
  pre-existing data. This should make the format of the files used and
  the workflow clearer. However, we recognise that this is not optimal
  for several reasons. We therefore plan to provide a vignette and
  example datasets in the near future.

- Since nearly all of the functions in the package generate output
  files, all of the examples presented in the documentation pages write
  files to the user’s tempdir() during their execution. At the end of
  each example, the code needed to clean up the user’s temp folder is
  provided.

### References

- Dongyue Hou, Hanbo Lin, Yuhan Feng, et al. CMAUP database update 2024:
  extended functional and association information of useful plants for
  biomedical research. Nucleic Acids Research 2024;
  <DOI:doi.org/10.1093/nar/gkad921>

- Webb, C. O.; Ackerly, D. D. & Kembel, S. W. (2008) Phylocom: software
  for the analysis of community structure and trait evolution.
  Bioinformatics 24: 2098-2100. (PDF)

- Xian Zeng, Peng Zhang, Yali Wang, et al. CMAUP: a database of
  collective molecular activities of useful plants. Nucleic Acids
  Research 2019; 47(D1): D1118-D1127; <DOI:doi.org/10.1093/nar/gky965>

- Zecca, G., Toini, E., Labra, M, Grassi, F. (2025) Accelerating the
  identification and the prioritisation of new plants with medicinal
  potential: the pm4mp R package.
