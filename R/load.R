# do this once
# devtools::install_bitbucket("remkoduursma/plotby")
library(plotBy)

r <- require(HIEv)
if(!r)stop("Install the HIEv R package from bitbucket.org/remkoduursma/hiev")

Library <- function(pkg, ...){
  
  PACK <- .packages(all.available=TRUE)
  pkgc <- deparse(substitute(pkg))
  
  if(pkgc %in% PACK){
    library(pkgc, character.only=TRUE)
  } else {
    install.packages(pkgc, ...)
    library(pkgc, character.only=TRUE)
  }
  
}


Library(stringr)
Library(doBy)
Library(plyr)
Library(RColorBrewer)
Library(reshape2)

source("R/functions.R")

if(!dir.exists("Output"))dir.create("Output")
if(!dir.exists("cache"))dir.create("cache")

setToken(tokenfile="HIEv_token.txt")
