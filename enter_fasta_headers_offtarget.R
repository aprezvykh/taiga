#!/usr/bin/Rscript
args <- commandArgs()
library("Biostrings", quietly = TRUE,warn.conflicts = FALSE)
print("Inserting spacer indexes!")
dir <- args[6]
file <- "off_target.ngg.space.fasta"
df <- read.delim(file, stringsAsFactors = F, header = F)
numbers <- seq(1:length(df[grep(">", df$V1),]))
pams <- paste(">offtarget_gRNA_", numbers, sep = "")
sequences <- df[grep(">", df$V1, invert = T),]
sequences <- unique(unlist(strsplit(sequences, " ")))
system("touch offtarget.ngg.headers.fasta")
cat(paste(pams, '\n', sequences, sep = ''),sep = '\n',file = "offtarget.ngg.headers.fasta")
