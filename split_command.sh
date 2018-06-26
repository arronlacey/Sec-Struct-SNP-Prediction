#!/bin/bash
split -d -l 300 -a 3 --additional-suffix=.fasta $1.snp.fasta $_1
