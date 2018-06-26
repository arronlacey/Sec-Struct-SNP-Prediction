#!/bin/bash
sed -e "s/myfile/$1/g" -e "s/23:00/$3/g" -e "s/myfolder/${1%.fasta}/g" < psibatch.sh | bsub -N -u $2
