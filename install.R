#!/usr/bin/Rscript
source("https://bioconductor.org/biocLite.R")
args <- commandArgs()

inst.spec <- args[6]
inst.spec <- "yeast"

options(warn = -1)
cat("This script will automatically install R packages!", sep = "\n")
cat("Avallible species: yeast, mouse, human, rat, fly, worm, all (install all species from list)", sep = "\n")
cat("Warning! Stains will work only with ENSEMBL gtf files", sep = "\n")

spec.df <- data.frame(db = c("org.Sc.sgd.db", "org.Mm.eg.db" ,"org.Hs.eg.db", "org.Rn.eg.db", "org.Dm.eg.db", "org.Ce.eg.db"),
                      spec = c("yeast", "mouse", "human", "rat", "fly", "worm"),
                      stringsAsFactors = F)


instdb <- spec.df[spec.df$spec == inst.spec,]$db
if(identical(instdb, character(0))){
  cat("Stain is not set correctly! Avallible species: yeast, mouse, human, rat, fly, worm", sep = "\n")
} else {
  cat(paste("Installing ", instdb), sep = "\n")
  biocLite(instdb)
}

if(tolower(inst.db) == "all"){
  cat("Installing all annotation packages!", sep = "")
  for(i in spec.df$db){
    cat(paste("Installing", i), sep = "\n")
    biocLite(i)
  }
}

local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org" 
       options(repos=r)
})

need.pack <- c("rtracklayer","stringr", "parallel", "plyr", "dplyr", "Biostrings", "ggplot2", "gridExtra", "xtable", "AnnotationDbi")
inst <- installed.packages()

for (f in need.pack){
  g <- grep(f, inst[,1])
  if(length(g)<1){
    cat(paste(f, "is not found in you system! Installing..."), sep = "\n")
    install.packages(f)
  } else {
    cat(paste(f, "is found!"), sep = "\n")
  }
  
}


cat("All R packages installed!", sep = "\n")



