**ROPSIR**
This software is made to find degenerate CRISPR-CAS9 gRNA targets in genome. It was developed a software to test 
different modificated Cas9 types. Main aim of ropsir is ... <br/>
**Software requirements:** <br/>
*Linux system* (tested in Ubuntu 16.04, kernel version 4.15.0-30-generic, 64 cores, 1TB RAM) <br/>
*Multicore* (8+) <br/>
*Ncbi-blast+* (install - sudo apt-get install ncbi-blast+ on Ubuntu/Debian systems) <br/>
*R language* (install - sudo apt-get install r-base-core Rscript) <br/>
*blastxmlparser* (install -  sudo apt-get install ruby ruby-dev; sudu apt-get install gem; sudo gem install blastxmlparser) <br/>
*RNAfold* (install - sudo apt-add-repository ppa:j-4/vienna-rna; sudo apt-get update; sudo apt-get install vienna-rna) <br/>
*ssconvert* (optional, converts csv file to xls, install - sudo apt-get install gnumeric) <br/>
*gffread* (need to convert genome file to CDS sequences) - see https://github.com/gpertea/gffread <br/>


**R packages: (will be installed automatically, by run install.R)** <br/>
*Biostrings* <br/>
*rtracklayer* <br/>
*stringr* <br/>
*parallel* <br/>
*plyr* <br/>
*dplyr* <br/>


**Run options**
-g - genome file in FASTA file, file <br/>
-a - annotation file in GTF format, file <br/>
-s - length of a spacer sequence exclude PAM site (20, by default), integer <br/>
-u - sequence that not allowed in spacer sequence (TTT, by default), string <br/>
-p - PAM sequence, that unallowed (AA, by default), string <br/>
-t - threads number, integer <br/>
-w - size of the word in blast, integer (10 by default, for debugging purposes only) <br/>
-d - uses only first 10 PAM sequences (T/F), bool <br/>
-pr - prefix for output csv file (string, date by default) <br/>
-sm - number of mismatches in seed region, integer (2 by default) <br/>
-nm - number of mismatches in non-seed region, integer (4 by default) <br/>
-tg - test specific gRNA sequence to find targets (GATTATAATATTCCTTGTGTTAG, for example), string <br/>
-pc - use only protein-coding sequences, ignore intergenic region (T/F), bool <br/>
-ts - test specific gene to find gRNA for its sequence (gene identifyer format should be <br/>
from "gene_id" column from gtf), character <br/>
-o - search for paralogs for gene (T/F) - works with -ts option, script will try to find paralogs for gene specifyed in -ts option <br/>
-c - paralogs cutoff (evalue scoring cutoff for paralogs blast) - float, 0.05 by default (should be set greater then 0.05, if you cannot find paralogs in paralog search mode) <br/>


**Minimum required arguments**
-g - genome file in FASTA file, file <br/>
-t - number of threads, file <br/>
-a - annotation file in GTF format, file <br/>
-pc - use only protein-coding sequences, ignore intergenic region (T/F), bool <br/>
-o - search for paralogs for gene (T/F)
-d - debug mode (T/F)

**Sample runs:** <br/>
*git clone https://github.com/aprezvykh/ropsir* <br/>
*sudo install.R* <br/>

1) Running ropsir to find gRNAs in genome-wide mode (parsing genome file by regular expression, with flags -sm;-nm;-u;-p;-s and finding **all** gRNAs and their targets)  <br/>
This mode can be run in protein-coding mode (-pc T), and all-genome mode (-pc F). This mode is very computational resource-demanding, and should be used in parallel (at least, in 32 threads), <br/>
otherwise, it can take a lot of time <br/>
*~/ropsir_location/./ropsir.sh -g ~/ropsir_location/data/genome.fasta -a ~/ropsir_location/data/genome.gtf -t 32 -d T -t 32 -pc T -o F -pr test.1* <br/>

2) Running ropsir to find target genes for specific gRNA, genome-wide (this mode takes nucleotide sequence from -tg flag, and find targets for ) <br/>
This mode also can be run in protein-coding mode (-pc T), and all-genome mode (-pc F). <br/>
*~/ropsir_location/./ropsir.sh -g ~/ropsir_location/data/genome.fasta -a ~/ropsir_location/data/genome.gtf -t 32 -d T -t 32 -pc F -tg GGTAAGATGAAGGAAACTGCCGA -o F -pr test.5 <br/>

3) To run all tests, just run:  <br/>
*./self_tests.sh* <br/>

Output of this script is table, that presented in csv and html format; it should look like: <br/>
![alt text](https://github.com/aprezvykh/ropsir/blob/master/sample_images/ropsir_image.PNG) <br/>
Column names explained: <br/>
1) gRNA.id - id of guide RNA found <br/>
2) gRNA.alignment.length - length of aligned guide RNA sequence <br/>
3) bitscore - Bitscore parametr from blastn software <br/>
4) evalue - Evalue parametr from blastn software <br/>
5) Aligned.sequence - aligned part of guide RNA <br/>
6) cigar.string - string; showing matches and mismatches of currend guide RNA, length 23 by default (| means match, X means mismatch) <br/>
7) total.mismatch.N - total number of mismatches in alignment <br/>
8) mismatch.position - position of mismatch (from left to right) <br/>
9) validation - shows alignment validation (from mismatch number in seed and non-seed region) <br/>
10) Locus - shows position of guide RNA alignment (gene_id or intergenic) <br/>
11) gRNA.Energy - free fold energy of guide RNA <br/>
12) GC.content - GC content :) <br/>
13) Genomic.coordinate - genomic coordinate in alignment <br/>


If you found an issue, please, report it in current repository or email me: <br/>
*aprezvykh@yandex. ru*
