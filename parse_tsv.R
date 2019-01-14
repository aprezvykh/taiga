#!/usr/bin/Rscript
library(Rcpp, warn.conflicts = FALSE, quietly = TRUE)
library(rtracklayer,warn.conflicts = FALSE,quietly = TRUE)
library(stringr,warn.conflicts = FALSE,quietly = TRUE)
library(parallel,warn.conflicts = FALSE,quietly = TRUE)
library(plyr,warn.conflicts = FALSE,quietly = TRUE)
library(dplyr,warn.conflicts = FALSE,quietly = TRUE)
library(xtable,warn.conflicts = FALSE,quietly = TRUE)

args <- commandArgs()

dir <- args[6]
gtf.path <- args[7]
prefix <- args[8]
threads <- args[9]
seed.mismatch <- args[10]
non.seed.mismatch <- args[11]
protein.coding <- args[12]
test.gene <- args[13]
files.dir <- args[14]
ropsir.dir <- args[15]
paralogs <- args[16]

seed.mismatch <- as.numeric(seed.mismatch)
non.seed.mismatch <- as.numeric(non.seed.mismatch)

cat(paste("Current data dir is", dir, sep = " "),sep = "\n")
cat(paste("GTF file is", gtf.path, sep = " "),sep = "\n")
cat(paste("Threads number is", threads, sep = " "),sep = "\n")
cat(paste("Seed mismatch number is", seed.mismatch, sep = " "),sep = "\n")
cat(paste("Non-seed mismatch number is", non.seed.mismatch, sep = " "),sep = "\n")
cat(paste("Protein coding is", protein.coding, sep = " "),sep = "\n")
cat(paste("Tested gene name is", test.gene, sep = " "),sep = "\n")

#dir <- c("~/ropsir.TESTING/")
#gtf.path <- c("~/git/ropsir/data/genome.gtf")
#prefix <- c("test.3")
#threads <- 32
#seed.mismatch <- 2
#non.seed.mismatch <- 4
#protein.coding <- "T"
#test.gene <- "YAL005C"
#paralogs <- "T"


cl <- makeCluster(threads,type = "FORK")
setwd(dir)

if(identical(tolower(paralogs), "f")){
    tab <- read.delim("blast.outfmt6",header = F,stringsAsFactors = F)
    if(nrow(tab) < 1){
      cat("No hits in blast! Exitting", sep = "\n")
      exit()
    } else {
      cat(paste("We got", nrow(tab), "blast hits!"), sep = "\n")
    }
    pam <- read.delim("blast.tsv",header = F, stringsAsFactors = F)
    df <- bind_cols(tab,pam)
}

if(identical(tolower(paralogs), "t")){
  tab <- read.delim("blast.outfmt6",header = F,stringsAsFactors = F)
  pam <- read.delim("blast.tsv",header = F, stringsAsFactors = F)
  tab.offtarget <- read.delim("blast.offtarget.outfmt6",header = F,stringsAsFactors = F)
  pam.offtarget <- read.delim("blast.offtarget.tsv",header = F, stringsAsFactors = F)
  
  target.df <- bind_cols(tab,pam)
  target.df$target <- c("target")
  
  offtarget.df <- bind_cols(tab.offtarget,pam.offtarget)
  offtarget.df$target <- c("offtarget")
  df <- bind_rows(target.df, offtarget.df)
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
  if(length(which(cc.cp < 10))> non.seed.mismatch){
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
construct.final.df <- function(x){
  sa <- df[df[["qseqid"]] == x,]
  en <- energies[energies[["name"]] == x,]
  sa[["energy"]] <- en[["val"]]
  sa[["pam.fasta"]]<- spacer.seqs[spacer.seqs[["headers"]] == x,][["seqs"]]
  sa <- data.frame(sa, stringsAsFactors = F)
  sa[["gc.content"]] <- letter.freq(unique(as.character(sa[["pam"]])))
  return(sa)
}
get.gene.coord <- function(x){
  sub <- gtf[grep(x, gtf$gene_id,fixed = T),]
  if(nrow(sub)<1){
    return("intergenic")
  }
  sub <- sub[sub$type == "gene",]
  return(paste(sub$seqnames, ":", sub$start, "-", sub$end, sep = ""))
}
exit <- function() {
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}


if(identical(test.gene, "nogene")){
  gtf <- as.data.frame(import(gtf.path))
} else if(nchar(test.gene) > 3){
  gtf <- as.data.frame(import(gtf.path))
  gtf.sub <- gtf[gtf$gene_id == test.gene,]
  cat(paste("Selected gene contains", nrow(gtf.sub), "features!", sep = " "), sep = "\n")
  if(nrow(gtf.sub) < 1){
    warning("Cannot find gene in GTF file!")
    cat(paste("Check your gene identifyer format! It should be like:", head(unique(gtf$gene_id),5), sep = ""),sep = "\n")
    exit()
  }
}

fasta <- read.delim("ngg.headers.fasta", header = F)
energies <- read.table("energies.txt", header = F)
headers <- as.character(fasta[grep(">", fasta$V1),])
seqs <- as.character(fasta[grep(">", fasta$V1, invert = T),])
spacer.seqs <- data.frame(headers,seqs)
spacer.seqs$headers <- gsub(">", "", spacer.seqs$headers)
energies$name <- spacer.seqs$headers
names(energies) <- c("val", "name")

if(identical(tolower(paralogs), "t")){
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
                   "sticks",
                   "target")
} else if(identical(tolower(paralogs), "f")){

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

}

if(identical(tolower(paralogs), "t")){
    df.tar <- df[df$target == "target",]
    df.offtar <- df[df$target == "offtarget",]
    tar.gene.coords <- data.frame(unlist(lapply(df.tar$sseqid, get.gene.coord)),stringsAsFactors = F)
    names(tar.gene.coords) <- c("coord")
    tar.gene.coords$chr <- unlist(lapply(strsplit(tar.gene.coords$coord,":"), function(x)x[[1]]))
    tmp <- unlist(lapply(strsplit(tar.gene.coords$coord,"-"), function(x)x[[1]]))
    tar.gene.coords$start <- unlist(lapply(strsplit(tmp, ":"), function(x)x[[2]]))
    tar.gene.coords$stop <- unlist(lapply(strsplit(tar.gene.coords$coord,"-"), function(x)x[[2]]))
    tar.gene.coords$true.start <- as.numeric(tar.gene.coords$start) + df.tar$sstart
    tar.gene.coords$true.end <- as.numeric(tar.gene.coords$stop) + df.tar$send
    df.tar$sseqid <- tar.gene.coords$chr
    df.tar$sstart <- tar.gene.coords$true.start
    df.tar$send <- tar.gene.coords$true.end
    df <- data.frame()
    df <- bind_rows(df.tar, df.offtar)
}


mm.sum <- as.numeric(seed.mismatch + non.seed.mismatch)
clusterExport(cl, "gtf")
cat("reconstructing cigar...", sep = "\n")
df$recon.cigar <- parApply(cl = cl, X = df, MARGIN = 1, FUN = reconstruct.cigar)
df$recon.cigar <- gsub(" ", "X", df$recon.cigar)
cat("Counting mismatches...", sep = "\n")
df$total.mm <- unlist(parLapply(cl = cl, X = df$recon.cigar, fun = count.total.mismatches))
df <- df[df$total.mm < mm.sum,]
cat("getting mismatch position...", sep = "\n")
df$mm.pos <- unlist(parLapply(cl = cl, X = df$recon.cigar,fun = get.mm.pos))
cat("parsing mismatch string...", sep = "\n")
df$val <- unlist(parLapply(cl = cl, X = df$mm.pos, fun = parse.mismatch.string))
cat("parsing annotation file. This can take a while...", sep = "\n")
df$loc <- parApply(cl = cl, X = df, MARGIN = 1, FUN = get.loci)
cat("Constructing final data frame!", sep = "\n")


grna.ids <- unique(df[["qseqid"]])
clusterExport(cl,varlist = list("df", "grna.ids", "energies", "spacer.seqs", "letter.freq"))
ll <- parLapply(cl = cl, X = grna.ids, fun = construct.final.df)
final.df <- data.frame(do.call("rbind", ll))
final.df <- final.df[grep("XXX$|XX$", final.df$recon.cigar, invert = T),]
final.df$mm.pos <- gsub("-1", "0", final.df$mm.pos)

if(identical(test.gene, "nogene")){
  big.final <- final.df
} else if (nchar(test.gene) > 0){
  big.final <- final.df[final.df$sseqid == test.gene,]
  if(nrow(final.df) < 1){
    warning("0 gRNAs found for your gene!")
  }
  cat(paste(nrow(final.df), "gRNAs found for your gene!", sep = " "),sep = "\n")
} 


final.df$aa <- NULL
final.df$bb <- NULL
final.df$cc <- NULL
final.df$dd <- NULL


final.df$gc.content <- as.numeric(as.character(unlist(lapply(as.character(final.df$pam.fasta), letter.freq))))
 
if(identical(tolower(paralogs), "t")){
  final.df$mismatch <- NULL
  final.df$gapopen <- NULL
  final.df$evalue <- NULL
  final.df$qstart <- NULL
  final.df$qend <- NULL
  final.df$sticks <- NULL
  final.df$pident <- NULL
  final.df$loc <- NULL
  big.final <- final.df
  names(big.final) <- c("gRNA.id","gene.id","gRNA.alignment.length", "gRNA.alignment.start",
                        "gRNA.alignment.end", "bitscore", "evalue", "aligned.sequence", "target", "cigar.string", "total.mismatch.N",
                        "mismatch.position", "validation", "gRNA.energy", "PAM.sequence", "GC.content")
  write.csv(big.final, paste(prefix, "-results.csv", sep = ""))
  big.final.for.html <- big.final
  big.final.for.html[,1] <- NULL
  write.table(big.final.for.html, paste(prefix, "-results.tsv", sep = ""),quote = F,sep = "\t")
  stopCluster(cl = cl)
  exit()
}

###parsing final data frame
if(identical(test.gene, "nogene")){
    if(tolower(protein.coding) == "t"){
      final.df$mismatch <- NULL
      final.df$gapopen <- NULL
      final.df$evalue <- NULL
      final.df$qstart <- NULL
      final.df$qend <- NULL
      final.df$sticks <- NULL
      final.df$pident <- NULL
      final.df$loc <- NULL
      cat("Ordering data frame...", sep = "\n")
      order.highest <- as.character(data.frame(table(t(final.df$qseqid)))[order(data.frame(table(t(final.df$qseqid)))$Freq, decreasing = T),]$Var1)
      big.final <- data.frame()
      for(f in order.highest){
        cat(f, sep = "\n")
        df <- final.df[final.df$qseqid == f,]
        big.final <- rbind(df, big.final)
      }
      
      big.final <- big.final[seq(dim(big.final)[1],1),]
      names(big.final) <- c("gRNA.id","gene.id","gRNA.alignment.length", "gRNA.alignment.start",
                            "gRNA.alignment.end", "bitscore", "evalue", "aligned.sequence", "cigar.string", "total.mismatch.N",
                            "mismatch.position", "validation", "gRNA.energy", "PAM.sequence", "GC.content")
      #big.final$PAM.sequence <- as.character(big.final$PAM.sequence)
      #big.final$GC.content <- as.character(unlist(lapply(big.final$PAM.sequence, letter.freq)))
    } else if (tolower(protein.coding) == "f") {
      final.df$mismatch <- NULL
      final.df$gapopen <- NULL
      final.df$evalue <- NULL
      final.df$qstart <- NULL
      final.df$qend <- NULL
      final.df$sticks <- NULL
      final.df$pident <- NULL
      final.df$coord <- paste(final.df$sseqid, ":", final.df$sstart, "-", final.df$send, sep = "")
      #final.df$coord <- NULL
      final.df$sseqid <- NULL
      final.df$sstart <- NULL
      final.df$send <- NULL
      final.df$sticks <- NULL
      cat("Ordering data frame...", sep = "\n")
      order.highest <- data.frame(table(t(final.df$qseqid)))[order(data.frame(table(t(final.df$qseqid)))$Freq, decreasing = T),]$Var1
      big.final <- data.frame()
      for(f in order.highest){
        cat(f, sep = "\n")
        df <- final.df[final.df$qseqid == f,]
        big.final <- rbind(df, big.final)
      }
      big.final <- big.final[seq(dim(big.final)[1],1),]
      names(big.final) <- c("gRNA.id","gRNA.alignment.length", "bitscore", "evalue",
                            "aligned.sequence", "cigar.string", "total.mismatch.N",
                            "mismatch.position", "validation", "Locus", "gRNA.energy", "PAM.sequence", "GC.content", "Genome.cordinate")
      big.final$PAM.sequence <- as.character(big.final$PAM.sequence)
      big.final$GC.content <- as.character(unlist(lapply(big.final$PAM.sequence, letter.freq)))
    }
}


write.csv(big.final, paste(prefix, "-results.csv", sep = ""))
html.filename <- paste(files.dir, "/", prefix, "-output.html", sep = "")
print(html.filename)
cat("Saving output to HTML table!")
sink(html.filename)
print(xtable(big.final), type = "html")
sink()

###report
if(identical(tolower(paralogs), "f")){
  top.grnas <- as.character(unique(big.final$gRNA.id)[1:20])  
  big.final.cutted <- big.final[big.final$gRNA.id %in% top.grnas,]
  write.csv(big.final.cutted, paste(prefix, "-cutted-results.csv", sep = ""))
  cat("Converting to XLS! (ssconvert warning about X11 display is non-crucial, just skip it :) )", sep = "\n")
  system(paste("ssconvert ", prefix, "-cutted-results.csv ", prefix, "-cutted-results.xls", sep = ""))
}

stopCluster(cl = cl)



