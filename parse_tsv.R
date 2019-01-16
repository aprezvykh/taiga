#!/usr/bin/Rscript
shhh <- suppressPackageStartupMessages 
shhh(library(Rcpp, warn.conflicts = FALSE, quietly = TRUE))
shhh(library(rtracklayer,warn.conflicts = FALSE,quietly = TRUE))
shhh(library(stringr,warn.conflicts = FALSE,quietly = TRUE))
shhh(library(parallel,warn.conflicts = FALSE,quietly = TRUE))
shhh(library(plyr,warn.conflicts = FALSE,quietly = TRUE))
shhh(library(dplyr,warn.conflicts = FALSE,quietly = TRUE))
shhh(library(xtable,warn.conflicts = FALSE,quietly = TRUE))
shhh(library(ggplot2,warn.conflicts = FALSE,quietly = TRUE))
shhh(library(AnnotationDbi,warn.conflicts = FALSE,quietly = TRUE))

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
stain <- args[17]

seed.mismatch <- as.numeric(seed.mismatch)
non.seed.mismatch <- as.numeric(non.seed.mismatch)

cat(paste("Current data dir is", dir, sep = " "),sep = "\n")
cat(paste("GTF file is", gtf.path, sep = " "),sep = "\n")
cat(paste("Threads number is", threads, sep = " "),sep = "\n")
cat(paste("Seed mismatch number is", seed.mismatch, sep = " "),sep = "\n")
cat(paste("Non-seed mismatch number is", non.seed.mismatch, sep = " "),sep = "\n")
cat(paste("Protein coding is", protein.coding, sep = " "),sep = "\n")
cat(paste("Tested gene name is", test.gene, sep = " "),sep = "\n")

debug = F

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



cat("Starting cluster", sep = "\n")
cl <- makeCluster(threads,type = "FORK")
cat("Cluster started!", sep = "\n")

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
    names(df) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                   "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                   "aa","bb", "cc", "dd", "ee", "seq", "sticks", "target")
} else if(identical(tolower(paralogs), "f")){
  names(df) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                 "aa","bb", "cc", "dd", "ee", "seq", "sticks")
  
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

for(i in c("aa", "bb", "cc", "dd")){
  final.df[[i]] <- NULL
}

final.df$gc.content <- as.numeric(as.character(unlist(lapply(as.character(final.df$pam.fasta), letter.freq))))

###parsing final data frame
if(identical(test.gene, "nogene")){
    if(tolower(protein.coding) == "t"){
      for(i in c("mismatch", "gapopen", "evalue", "qstart",
                 "qend", "sticks", "pident", "loc")){
        final.df[[i]] <- NULL
      }
      
      cat("Ordering data frame...", sep = "\n")
      order.highest <- as.character(data.frame(table(t(final.df$qseqid)))[order(data.frame(table(t(final.df$qseqid)))$Freq, decreasing = T),]$Var1)
      big.final <- data.frame()
      for(f in order.highest){
        df <- final.df[final.df$qseqid == f,]
        df <- df[order(df$total.mm, decreasing = F),]
        df$num.full.genes <- length(unique(df[df$recon.cigar == "|||||||||||||||||||||||",]$sseqid))
        df$diff.mut <- length(unique(df$recon.cigar)) -1
        big.final <- rbind(df, big.final)
        
      }
      
      big.final <- big.final[seq(dim(big.final)[1],1),]
      names(big.final) <- c("gRNA.id","gene.id","gRNA.alignment.length", "gRNA.alignment.start",
                            "gRNA.alignment.end", "bitscore", "evalue", "aligned.sequence", "cigar.string", 
                            "total.mismatch.N", "mismatch.position", "validation", "gRNA.energy", "PAM.sequence", 
                            "GC.content", "Number.of.genes.with.full.match", "Number.of.different.mutation.locations")
    } else if (tolower(protein.coding) == "f") {
      final.df$coord <- paste(final.df$sseqid, ":", final.df$sstart, "-", final.df$send, sep = "")
      for(i in c("mismatch", "gapopen", "evalue", "qstart",
                 "qend", "sticks", "pident", "sseqid", 
                 "sstart", "send")){
        final.df[[i]] <- NULL
      }
      cat("Ordering data frame...", sep = "\n")
      order.highest <- as.character(data.frame(table(t(final.df$qseqid)))[order(data.frame(table(t(final.df$qseqid)))$Freq, decreasing = T),]$Var1)
      big.final <- data.frame()
      for(f in order.highest){
        df <- final.df[final.df$qseqid == f,]
        df <- df[order(df$total.mm, decreasing = F),]
        df$num.full.genes <- length(unique(df[df$recon.cigar == "|||||||||||||||||||||||",]$loc))
        df$diff.mut <- length(unique(df$recon.cigar)) -1
        big.final <- rbind(df, big.final)
        
      }
      
      big.final <- big.final[seq(dim(big.final)[1],1),]
      names(big.final) <- c("gRNA.id","gRNA.alignment.length", "bitscore", "evalue",
                            "aligned.sequence", "cigar.string", "total.mismatch.N",
                            "mismatch.position", "validation", "Locus", "gRNA.energy", 
                            "PAM.sequence", "GC.content", "Genome.cordinate", 
                            "Number.of.genes.with.full.match", "Number.of.different.mutation.locations")
      big.final$PAM.sequence <- as.character(big.final$PAM.sequence)
      big.final$GC.content <- as.character(unlist(lapply(big.final$PAM.sequence, letter.freq)))
    }
}

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




write.csv(big.final, paste(prefix, "-results.csv", sep = ""))
html.filename <- paste(files.dir, "/", prefix, "-output.html", sep = "")
cat("Saving output to HTML table!", sep = "\n")
sink(html.filename)
print(xtable(big.final), type = "html")
sink()

###report
if(identical(tolower(paralogs), "f")){
  top.20.grnas <- as.character(unique(big.final$gRNA.id)[1:20])  
  big.final.cutted.20 <- big.final[big.final$gRNA.id %in% top.20.grnas,]
  write.csv(big.final.cutted.20, paste(prefix, "top20-cutted-results.csv", sep = ""))
  
  top.50.grnas <- as.character(unique(big.final$gRNA.id)[1:50])  
  big.final.cutted.50 <- big.final[big.final$gRNA.id %in% top.50.grnas,]
  write.csv(big.final.cutted.50, paste(prefix, "top50-cutted-results.csv", sep = ""))
  
  top.200.grnas <- as.character(unique(big.final$gRNA.id)[1:200])  
  big.final.cutted.200 <- big.final[big.final$gRNA.id %in% top.200.grnas,]
  write.csv(big.final.cutted.200, paste(prefix, "top200-cutted-results.csv", sep = ""))
  
  cat("Converting to XLS", sep = "\n")
  system(paste("ssconvert ", prefix, "top20-cutted-results.csv ", prefix, "top20-cutted-results.xls 2> /dev/null", sep = ""))
  system(paste("ssconvert ", prefix, "top50-cutted-results.csv ", prefix, "top50-cutted-results.xls 2> /dev/null", sep = ""))
  system(paste("ssconvert ", prefix, "top200-cutted-results.csv ", prefix, "top200-cutted-results.xls 2> /dev/null", sep = ""))
  
  system("rm *-cutted-results.csv")
}


cat("Picking best gRNAs!", sep = "\n")
v.sorted <- sort(unique(big.final$Number.of.different.mutation.locations))
best.grna.df <- big.final[which(big.final$Number.of.different.mutation.locations == max(big.final$Number.of.different.mutation.locations) & 
                  big.final$Number.of.genes.with.full.match == 1),]
write.csv(best.grna.df, paste(prefix, "best_gRNA-results.csv", sep = ""))
system(paste("ssconvert ", prefix, "best_gRNA-results.csv ", prefix, "top200-cutted-results.xls 2> /dev/null", sep = ""))
system(paste("rm ", prefix, "best_gRNA-results.csv", sep = ""))

cat("Making graphic report!", sep = "\n")
if(identical(tolower(protein.coding), "f")){
    pdf(paste(prefix, "-graphic-report.pdf", sep = ""))
    top.genes <- as.data.frame(table(t(big.final$Locus)))
    top.genes <- top.genes[order(top.genes$Freq,decreasing = T),][1:10,]
    g <- ggplot(data = top.genes) + geom_bar(aes(x = reorder(Var1, -Freq), y = Freq), stat = "identity") + 
      theme_bw() + ggtitle("Top-10 genes with gRNAs hit") + 
      scale_x_discrete("Gene") + scale_y_continuous("Number of hits")
    print(g)
    plot(density(as.numeric(big.final$GC.content)), main = "GC content density")
    plot(density(as.numeric(big.final$gRNA.energy)), main = "gRNA free energy density")
    dev.off()
} else if(identical(tolower(protein.coding), "t")){
  pdf(paste(prefix, "-graphic-report.pdf", sep = ""))
  top.genes <- as.data.frame(table(t(big.final$gene.id)))
  top.genes <- top.genes[order(top.genes$Freq,decreasing = T),][1:10,]
  g <- ggplot(data = top.genes) + geom_bar(aes(x = reorder(Var1, -Freq), y = Freq), stat = "identity") + 
    theme_bw() + ggtitle("Top-10 genes with gRNAs hit") + 
    scale_x_discrete("Gene") + scale_y_continuous("Number of hits")
  print(g)
  plot(density(as.numeric(big.final$GC.content)), main = "GC content density")
  plot(density(as.numeric(big.final$gRNA.energy)), main = "gRNA free energy density")
  dev.off()
}


stopCluster(cl = cl)
cat("Done!", sep = "\n")






