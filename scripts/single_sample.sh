#!/bin/bash

# Usage: ./analyze_sample.sh <sample> 

# assign input arguments to variables
sample=$1

export PATH=/home/user04/apps/samtools:$PATH

# prepare genome.fa
bowtie-build genome.fa genome

# align reads to the reference genome using Tophat
tophat -p 10 -I 1000 -i 20 --bowtie1 --library-type fr-firststrand -o tophat.$sample.dir genome ${sample}_1_trimmed.fastq.gz ${sample}_2_trimmed.fastq.gz

# rename the alignment (bam) output file
mv tophat.$sample.dir/accepted_hits.bam tophat.$sample.dir/$sample.bam

# index the BAM file using samtools
samtools index tophat.$sample.dir/$sample.bam

# reconstruct transcripts for this sample using cufflinks
cufflinks --overlap-radius 1 --library-type fr-firststrand -o cufflinks.$sample.dir tophat.$sample.dir/$sample.bam

# renaming the cufflinks transcript structure output file
mv cufflinks.$sample.dir/transcripts.gtf cufflinks.$sample.dir/$sample.transcripts.gtf

echo cufflinks.$sample.dir/$sample.transcripts.gtf >> assemblies.txt

# message of analysis completion
echo "Analysis for $sample completed."