#!/usr/bin/Rscript
options(warn = -1)
print("This script will automatically install R packages!")

local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org" 
       options(repos=r)
})

need.pack <- c("rtracklayer","stringr", "parallel", "plyr", "dplyr", "Biostrings", "ggplot2", "gridExtra", "xtable")
inst <- installed.packages()

for (f in need.pack){
  g <- grep(f, inst[,1])
  if(length(g)<1){
    print(paste(f, "is not found in you system! Installing..."))
    install.packages(f)
  } else {
    print(paste(f, "is found!"))
  }
  
}


print("All R packages installed!")
