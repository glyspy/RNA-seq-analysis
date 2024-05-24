#!/bin/bash

# list of samples
samples=("SRR1982462" "SRR1982463" "SRR1982464" "SRR1982465" "SRR1982466" "SRR1982467" "SRR1982468")
#export PATH=/home/user04/apps/samtools:$PATH

# loop through the samples and run the script of the analysis
for sample in "${samples[@]}"
do
        echo "Analyzing sample: $sample ..."
        sh ./single_sample.sh $sample
done
