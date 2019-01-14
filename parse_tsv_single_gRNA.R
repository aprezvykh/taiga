#!/usr/bin/Rscript
args <- commandArgs()
library(rtracklayer)
library(stringr)
library(parallel)
library(plyr)
library(dplyr)

exit <- function() {
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}


dir <- args[6]
gtf.path <- args[7]
prefix <- args[8]
threads <- args[9]
seed.mismatch <- args[10]
non.seed.mismatch <- args[11]
protein.coding <- args[12]

#dir <- "~/ropsir.TESTING/"
#gtf.path <- "~/git/ropsir/data/genome.gtf"
#prefix <- "test.4"
#threads <- 32
#seed.mismatch <- 2
#non.seed.mismatch <- 4
#protein.coding <- "F"

seed.mismatch <- as.numeric(seed.mismatch)
non.seed.mismatch <- as.numeric(non.seed.mismatch)
###
cl <- makeCluster(threads,type = "FORK")
setwd(dir)

print(paste("Using", threads, "threads!"))

tab <- read.delim("blast.outfmt6",header = F,stringsAsFactors = F)
pam <- read.delim("blast.tsv",header = F, stringsAsFactors = F)
df <- bind_cols(tab,pam)
nrow(df) < 1

if(nrow(df) < 2){
  cat("Check your gRNA! No hits found!", sep = "\n")
  exit()
}

letter.freq <- function(x){
  ss <- summary(as.factor(unlist(strsplit(x, NULL))))
  gc <- 100*((ss[3] + ss[2])/23)
  return(gc)
}

get.mm.pos <- function(x){
  paste(as.character(gregexpr(pattern ='X',x)[[1]]), collapse = ",")
}

parse.mismatch.string <- function(x){
  cc.cp <- as.numeric(strsplit(x, ",")[[1]])
  if(length(which(cc.cp < 10)) > non.seed.mismatch){
    return("novalid")
  } else if(length(which(cc.cp > 10)) > seed.mismatch){
    return("novalid")
  } else {
    return("valid")
  }
  
}

reconstruct.cigar <- function(x){
  st <- x[["sticks"]]
  ll <- seq(x[["qstart"]], x[["qend"]])
  ll.min <- as.numeric(min(ll))
  ll.max <- as.numeric(max(ll))
  right.ll <- length(seq(ll.max,23))
  left.ll <- length(seq(1,ll.min))
  pp <- paste(paste(rep(" ", left.ll), sep = "", collapse = ""),
              st,
              paste(rep(" ", right.ll),sep = "", collapse = ""),
              collapse = "", sep = "")
  pp <- str_sub(pp, start=2)
  pp <- gsub('.$', "", pp)
  return(pp)
}


count.total.mismatches <- function(x){
  str_count(x,pattern = "X")
}

get.loci <- function(x){
  sgtf <- gtf[gtf[["seqnames"]] == x[["sseqid"]],]
  sub.sgtf <- sgtf[sgtf[["start"]] <= x[["sstart"]],]
  sub.sgtf <- sub.sgtf[sub.sgtf[["end"]] >= x[["send"]],]
  ret.gene <- unique(sub.sgtf[["gene_id"]])
  ret.regions <- paste(as.character(unique(sub.sgtf[["type"]])),collapse = ",")
  if(identical(ret.gene, character(0))){
    ret.gene <- c("intergenic")
  } else {
    ret.gene <- unique(sub.sgtf[["gene_id"]])
  }
  return(paste(ret.gene[1]))
}


print("importing GTF...")
gtf <- as.data.frame(import(gtf.path))
fasta <- read.delim("testgrna.fasta", header = F)
energies <- read.table("energies.txt", header = F)
headers <- c("test-gRNA")
seqs <- fasta$V1
spacer.seqs <- data.frame(headers,seqs)
spacer.seqs$headers <- gsub(">", "", spacer.seqs$headers)
energies$name <- spacer.seqs$headers
names(energies) <- c("val", "name")

names(df) <- c("qseqid",
               "sseqid",
               "pident",
               "length",
               "mismatch",
               "gapopen", 
               "qstart", 
               "qend", 
               "sstart", 
               "send", 
               "evalue", 
               "bitscore",
               "aa",
               "bb",
               "cc",
               "dd",
               "ee",
               "seq",
               "sticks")
df$ee <- NULL

mm.sum <- as.numeric(seed.mismatch + non.seed.mismatch)
clusterExport(cl, "gtf")
print("reconstructing cigar...")
df$recon.cigar <- parApply(cl = cl, X = df, MARGIN = 1, FUN = reconstruct.cigar)
df$recon.cigar <- gsub(" ", "X", df$recon.cigar)
print("Counting mismatches...")

df$total.mm <- unlist(lapply(df$recon.cigar, count.total.mismatches))
#df <- df[df$total.mm < mm.sum,]
print("getting mismatch position...")
df$mm.pos <- unlist(parLapply(cl = cl, X = df$recon.cigar,fun = get.mm.pos))
print("parsing mismatch string...")
df$val <- unlist(parLapply(cl = cl, X = df$mm.pos, fun = parse.mismatch.string))
print("parsing annotation file...")
df$loc <- parApply(cl = cl, X = df, MARGIN = 1, FUN = get.loci)
print("Constructing final data frame!")
final.df <- data.frame()


df$qseqid <- c("test-gRNA")
f <- unique(df$qseqid)
sa <- df[df$qseqid == f,]
en <- energies[energies$name == f,]
sa$energy <- en$val
sa$pam.fasta <- spacer.seqs[spacer.seqs$headers == f,]$seqs
sa <- data.frame(sa, stringsAsFactors = F)
sa$gc.content <- letter.freq(unique(as.character(sa$pam.fasta)))
final.df <- rbind(sa, final.df)


final.df <- final.df[grep("XXX$|XX$", final.df$recon.cigar, invert = T),]
final.df$mm.pos <- gsub("-1", "0", final.df$mm.pos)

final.df$aa <- NULL
final.df$bb <- NULL
final.df$cc <- NULL
final.df$dd <- NULL

final.df$coord <- paste(final.df$sseqid, ":",final.df$sstart, "-", final.df$send, sep = "")
final.df$pident <- NULL
final.df$length <- NULL
final.df$mismatch <- NULL
final.df$qstart <- NULL
final.df$qend <- NULL
final.df$sstart <- NULL
final.df$send <- NULL
final.df$sticks <- NULL
final.df$gapopen <- NULL

names(final.df) <- c("gRNA.id", "chr", "evalue", "bitscore",
                     "Aligned.sequence", "cigar.string", "total.mismatch.N", "mismatch.position",
                     "validation", "Locus", "gRNA.energy", "PAM.sequence", "GC.content", "Genomic.coordinate")

write.csv(final.df, paste(prefix, "-single-gRNA-results.csv", sep = ""))
stopCluster(cl = cl)
