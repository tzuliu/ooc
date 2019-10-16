# OOC: Ordered Optimal Classification 
[![Build Status](https://travis-ci.org/tzuliu/ooc.svg?branch=master)](https://travis-ci.org/tzuliu/ooc)
![license](https://img.shields.io/github/license/mashape/apistatus.svg)
[![Github All Releases](https://img.shields.io/github/downloads/tzuliu/ooc/total.svg)]()
[![GitHub Release](https://github-basic-badges.herokuapp.com/release/tzuliu/ooc.svg)]()
[codecov](https://codecov.io/github/tzuliu/ooc/branch/master/graphs/badge.svg)](https://codecov.io/gh/tzuliu/ooc) 
<!--[![Github all releases](https://img.shields.io/github/downloads/tzuliu/ooc/total.svg)](https://GitHub.com/tzuliu/ooc/releases/)-->
<!--[![GitHub Download Count](https://github-basic-badges.herokuapp.com/downloads/tzuliu/ooc/total.svg)]()-->



This is a repository for the project "OOC: Ordered Optimal Classification"

---

## To install the "ooc" Package in R, first download and load "devtools" as below:
````
install.packages("devtools", dependencies=TRUE)
library(devtools)
````

## Then install and load the development version from GitHub:
````
devtools::install_github('tzuliu/ooc')
library(ooc)
````

---
## Prerequisites:

### Mac Users

* Make sure you have already installed ***gfortran*** (reference can be reached [here](https://cran.r-project.org/bin/macosx/tools/))
   - Trouble Shooting for Warning Message in R regarding ***gfortran*** not found:
      - Try **Homebrew** to install ***gfortran*** through ````brew install gcc```` and ````brew link gcc````
      - Try to open the file ***Makeconf*** in the folder ````/Library/Frameworks/R.framework/Resources/etc````, and then change **FLIBS** to ````FLIBS = -L/usr/local/Cellar/gcc/"YOUR GCC VERSION"/lib/gcc/"6, 7, 8, or 9 (depends on the version)"```` (reference can be reached [here](https://octaviancorlade.github.io/compile-rcpparmadillo-glibfortran-high-sierra/))
* Make sure you have already installed ***Xcode Developer Tools***

### Windows Users

* Make sure you have already installed ***RTools***
---
## Paper:

For the final version of the paper "What Ordered Optimal Classification Reveals about Ideological Structure Cleavages and Polarization" and the related/replication codes, please refer to [this repository](https://github.com/tzuliu/What-Ordered-Optimal-Classification-Reveals-about-Ideological-Structure-Cleavages-and-Polarization)
