#!/bin/bash
ropsir_dir=/mnt/raid/illumina/AlexR/git/ropsir/
exec_dir=$(readlink -e $(pwd))

test1="$ropsir_dir/./ropsir.sh -g $ropsir_dir/data/genome.fasta -a $ropsir_dir/data/genome.gtf -t 32 -d T -pc T -o F -pr test.1"
echo "Evaluating test 1! Command is $test1"
eval "$test1"

test2="$ropsir_dir/./ropsir.sh -g $ropsir_dir/data/genome.fasta -a $ropsir_dir/data/genome.gtf -t 32 -d T -pc F -o F -pr test.2"
echo "Evaluating test 2! Command is $test2"
eval "$test2"

test3="$ropsir_dir/./ropsir.sh -g $ropsir_dir/data/genome.fasta -a $ropsir_dir/data/genome.gtf -t 32 -o F -pc T -tg ATGTCAAAAGCTGTCGGTATTGA -pr test.3"
echo "Evaluating test 3! Command is $test3"
eval "$test3"

test4="$ropsir_dir/./ropsir.sh -g $ropsir_dir/data/genome.fasta -a $ropsir_dir/data/genome.gtf -t 32 -o F -pc T -tg ATGTCAAAAGCTGTCGGTATTGA -pr test.4"
echo "Evaluating test 5! Command is $test4"
eval "$test4"
