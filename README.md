**TAIGA** <br/>
![alt text](https://upload.wikimedia.org/wikipedia/commons/2/29/Archangelsk_taiga.JPG) <br/>
This software is made to find degenerate CRISPR-CAS9 gRNA targets in genome. It was developed a software to test 
different modificated Cas9 types. <br/>
**Software requirements:** <br/>
1) *Linux system* (tested in Ubuntu 16.04, kernel version 4.15.0-30-generic, 64 cores, 1TB RAM) <br/>
2) *Multicore* (8+) <br/>
3) *Ncbi-blast+* (install - sudo apt-get install ncbi-blast+ on Ubuntu/Debian systems) <br/>
4) *R language* (install - sudo apt-get install r-base-core Rscript) <br/>
5) *blastxmlparser* (install -  sudo apt-get install ruby ruby-dev; sudu apt-get install gem; sudo gem install blastxmlparser) <br/>
6) *RNAfold* (install - sudo apt-add-repository ppa:j-4/vienna-rna; sudo apt-get update; sudo apt-get install vienna-rna) <br/>
7) *ssconvert* (optional, converts csv file to xls, install - sudo apt-get install gnumeric) <br/>
8) *gffread* (need to convert genome file to CDS sequences) - see https://github.com/gpertea/gffread to install and compile. Gffread must be accessible from command prompt () <br/>

**R packages: (will be installed automatically, by run install.R)** <br/>
1) *Biostrings* <br/>
2) *rtracklayer* <br/>
3) *stringr* <br/>
4) *parallel* <br/>
5) *plyr* <br/>
6) *dplyr* <br/>


**Run options**
1) -g - genome file in FASTA file, file <br/>
2) -a - annotation file in GTF format, file <br/>
3) -s - length of a spacer sequence exclude PAM site (20, by default), integer <br/>
4) -u - sequence that not allowed in spacer sequence (TTT, by default), string <br/>
5) -p - PAM sequence, that unallowed (AA, by default), string <br/>
6) -t - threads number, integer <br/>
7) -w - size of the word in blast, integer (10 by default, **for debugging purposes only**) <br/>
8) -d - uses only first 10 PAM sequences (T/F), bool <br/>
9) -pr - prefix for output csv file (string, date by default) <br/>
10) -sm - number of mismatches in seed region, integer (2 by default) <br/>
11) -nm - number of mismatches in non-seed region, integer (4 by default) <br/>
12) -tg - test specific gRNA sequence to find targets (GATTATAATATTCCTTGTGTTAG, for example), string <br/>
13) -pc - use only protein-coding sequences, ignore intergenic region (T/F), bool <br/>
14) -ts - test specific gene to find gRNA for its sequence (gene identifyer format should be <br/>
from "gene_id" column from gtf) (**NOW UNDER DEVELOPMENT**), character <br/>
15) -o - search for paralogs for gene (T/F) - works with -ts option, script will try to find paralogs for gene specifyed in -ts option (**NOW UNDER DEVELOPMENT, SHOULD BE ALWAYS FALSE**) <br/>
16) -c - paralogs cutoff (evalue scoring cutoff for paralogs blast) - float, 0.05 by default (should be set greater then 0.05, if you cannot find paralogs in paralog search mode) (**NOW UNDER DEVELOPMENT**) <br/>


**Minimum required arguments**
1) -g - genome file in FASTA file, file <br/>
2) -t - number of threads, file <br/>
3) -a - annotation file in GTF format, file <br/>
4) -pc - use only protein-coding sequences, ignore intergenic region (T/F), bool <br/>
5) -o - search for paralogs for gene (T/F)
6) -d - debug mode (T/F)

___

**Sample runs:** <br/>
*git clone https://github.com/aprezvykh/taiga* <br/>
*install.R* (don't run it by sudo!) <br/>

1) Running taiga to find gRNAs in genome-wide mode (parsing genome file by regular expression, with flags -sm;-nm;-u;-p;-s and finding **all** gRNAs and their targets)  <br/>
This mode can be run in protein-coding mode (-pc T), and all-genome mode (-pc F). This mode is very computational resource-demanding, and should be used in parallel (at least, in 32 threads), <br/>
otherwise, it can take a lot of time <br/>
*~/taiga_location/./taiga.sh -g ~/taiga_location/data/genome.fasta -a ~/taiga_location/data/genome.gtf -t 32 -d T -t 32 -pc T -o F -pr test.1* <br/>

2) Running taiga to find target genes for specific gRNA, genome-wide (this mode takes nucleotide sequence from -tg flag, and find targets for ) <br/>
This mode also can be run in protein-coding mode (-pc T), and all-genome mode (-pc F). <br/>
*~/taiga_location/./taiga.sh -g ~/taiga_location/data/genome.fasta -a ~/taiga_location/data/genome.gtf -t 32 -d T -t 32 -pc F -tg GGTAAGATGAAGGAAACTGCCGA -o F -pr test.5 <br/>

3) To run all tests, just run:  <br/>
*./self_tests.sh* <br/>

**Output files:** <br/>
Output of this script in genome-wide mode is presented with several files: <br/>
1) *prefix-graphic-report.pdf* <br/>
2) *prefix-results.csv (full table* in CSV format, useful for analysis in R/Python)  <br/>
3) *prefix-output.html (full table* in HTML format) <br/>
4) *prefixtop20-cutted-results.xls* (top-20 gRNAs) <br/>
5) *prefixtop50-cutted-results.xls* (top-50 gRNAs)<br/>
6) *prefixtop200-cutted-results.xls* (top-200 gRNAs) <br/>

Output of this script in single gRNA-testing mode: <br/>
*prefix-single-gRNA.xlsx*  <br/>
___

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
___

If you found an issue, please, report it in current repository or email me: <br/>
*aprezvykh@yandex. ru*
