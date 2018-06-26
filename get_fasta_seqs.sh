#!/bin/bash

#download fasta seqs given file of uniprot ids

file=$1
name=$2

ids=($(cat ${file} | awk '{print $1}'))
pos=($(cat ${file} | awk '{print $2}'))
wild=($(cat ${file} | awk '{print $3}'))
sub=($(cat ${file} | awk '{print $4}'))


# get ref fasta for each line in file, with custom header attached
for i in "${!ids[@]}" ; do
      echo "#${ids[i]}_${pos[i]}_${wild[i]}_${sub[i]}"; curl -sS "http://www.uniprot.org/uniprot/"${ids[i]}".fasta";
done | sed '/^>/ d' | sed -r 's/[#]+/>/g' | perl -npe 'chomp if ($.!=1 && !s/^>/\n>/)' > $name.snp.fasta
