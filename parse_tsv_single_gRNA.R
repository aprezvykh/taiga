#!/usr/bin/Rscript
shhh <- suppressPackageStartupMessages 
shhh(library(rtracklayer, warn.conflicts = FALSE, quietly = TRUE))
shhh(library(stringr, warn.conflicts = FALSE, quietly = TRUE))
shhh(library(parallel, warn.conflicts = FALSE, quietly = TRUE))
shhh(library(plyr, warn.conflicts = FALSE, quietly = TRUE))
shhh(library(dplyr, warn.conflicts = FALSE, quietly = TRUE))

args <- commandArgs()

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
stain <- args[13]

debug = T

if(debug == T){
  dir <- c("~/taiga.TESTING/")
  gtf.path <- c("~/git/taiga/data/genome.gtf")
  prefix <- c("test.3")
  threads <- 32
  seed.mismatch <- 2
  non.seed.mismatch <- 4
  protein.coding <- "F"
  test.gene <- "nogene"
  paralogs <- "F"
  stain <- "yeast"
}

spec.df <- data.frame(db = c("org.Sc.sgd.db", "org.Mm.eg.db" ,"org.Hs.eg.db", "org.Rn.eg.db", "org.Dm.eg.db", "org.Ce.eg.db"),
                      spec = c("yeast", "mouse", "human", "rat", "fly", "worm"),
                      stringsAsFactors = F)

if(stain != "unknown"){
  cat(paste("Annotation stain is set to", stain), sep = "\n")
  cat(paste("Loading", stain), sep = "\n")
  instdb <- spec.df[spec.df$spec == stain,]$db
  shhh(library(paste0(instdb), warn.conflicts = FALSE,quietly = TRUE,character.only = T))
} else {
  cat("Stain is not set!", sep = "\n")
}


seed.mismatch <- as.numeric(seed.mismatch)
non.seed.mismatch <- as.numeric(non.seed.mismatch)
###
cat("Starting cluster", sep = "\n")
cl <- makeCluster(threads,type = "FORK")
cat("Cluster started!", sep = "\n")

setwd(dir)

cat(paste("Using", threads, "threads!"), sep = "\n")

tab <- read.delim("blast.outfmt6",header = F,stringsAsFactors = F)
pam <- read.delim("blast.tsv",header = F, stringsAsFactors = F)
df <- bind_cols(tab,pam)

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


cat("importing GTF...", sep = "\n")
gtf <- as.data.frame(import(gtf.path))
fasta <- read.delim("testgrna.fasta", header = F)
energies <- read.table("energies.txt", header = F)
headers <- c("test-gRNA")
seqs <- fasta$V1
spacer.seqs <- data.frame(headers,seqs)
spacer.seqs$headers <- gsub(">", "", spacer.seqs$headers)
energies$name <- spacer.seqs$headers
names(energies) <- c("val", "name")

names(df) <- c("qseqid", "sseqid", "pident", "length", 
               "mismatch", "gapopen", "qstart", "qend", 
               "sstart", "send", "evalue", "bitscore",
               "aa", "bb", "cc", "dd", "ee", "seq",
               "sticks")
df$ee <- NULL

mm.sum <- as.numeric(seed.mismatch + non.seed.mismatch)
clusterExport(cl, "gtf")
cat("reconstructing cigar...", sep = "\n")
df$recon.cigar <- parApply(cl = cl, X = df, MARGIN = 1, FUN = reconstruct.cigar)
df$recon.cigar <- gsub(" ", "X", df$recon.cigar)
cat("Counting mismatches...")

df$total.mm <- unlist(lapply(df$recon.cigar, count.total.mismatches))
#df <- df[df$total.mm < mm.sum,]
cat("getting mismatch position...", sep = "\n")
df$mm.pos <- unlist(parLapply(cl = cl, X = df$recon.cigar,fun = get.mm.pos))
cat("parsing mismatch string...", sep = "\n")
df$val <- unlist(parLapply(cl = cl, X = df$mm.pos, fun = parse.mismatch.string))
cat("parsing annotation file...", sep = "\n")
df$loc <- parApply(cl = cl, X = df, MARGIN = 1, FUN = get.loci)
cat("Constructing final data frame!", sep = "\n")
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

if(nrow(final.df < 2)){
  cat("Cannot found any valid gRNA! exitting!", sep = "\n")
  exit()
}

for(i in c("aa", "bb", "cc", "dd")){
  final.df[[i]] <- NULL
}

final.df$coord <- paste(final.df$sseqid, ":",final.df$sstart, "-", final.df$send, sep = "")

for(i in c("pident", "length", "mismatch", "qstart",
           "qend", "sstart", "send", "sticks", 
           "gapopen")){
  final.df[[i]] <- NULL
}

names(final.df) <- c("gRNA.id", "chr", "evalue", "bitscore",
                     "Aligned.sequence", "cigar.string", "total.mismatch.N", "mismatch.position",
                     "validation", "Locus", "gRNA.energy", "PAM.sequence", "GC.content", "Genomic.coordinate")

final.df <- final.df[order(final.df$total.mismatch.N,decreasing = F),]
final.df$Number.of.genes.with.full.match <- length(unique(final.df[final.df$cigar.string == "|||||||||||||||||||||||"]$Locus))
final.df$Number.of.different.mutation.locations <- length(unique(final.df$cigar.string))

if(stain != "unknown"){
    if(tolower(protein.coding) == "f"){
      big.final$gene.name <- mapIds(org.Sc.sgd.db,
                                    keys = as.character(big.final$Locus),
                                    column = "DESCRIPTION",
                                    keytype = "ENSEMBL",
                                    multiVals = "first")
      
    } else if(tolower(protein.coding) == "t"){
      big.final$gene.name <- mapIds(org.Sc.sgd.db,
                                    keys = as.character(big.final$gene.id),
                                    column = "DESCRIPTION",
                                    keytype = "ENSEMBL",
                                    multiVals = "first")
    }
}



write.csv(final.df, paste(prefix, "-single-gRNA-results.csv", sep = ""))
system(paste("ssconvert ", prefix, "-single-gRNA-results.csv ", prefix, "-single-gRNA-results.xls 2> /dev/null", sep = ""))
stopCluster(cl = cl)
cat("Done!", sep = "\n")

write.csv(final.df, paste(prefix, "-single-gRNA-results.csv", sep = ""))
stopCluster(cl = cl)
