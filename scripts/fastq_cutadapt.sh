#!/bin/bash

# create an ssociative array mapping sample names to adapter sequences

declare -A adapter_map

adapter_map=([SRR1982462_1]="GGGGTGAAATCCTTGCCCAGGTGGTGGCCCAGCACAATCACGATCATATTGCCCAGGAGCCTGAAGTTCTCAGGATCCACATGCAGCTTGTCACAGTGGAGCTCACTGAGGCTGGCAAAGGTGCCCTTGAGGCTGTCCAAGTGATTCAGGCCATCGTTAAAGGCAGTTATCACCTTCTTGCCATGGGCCTTCACTTTGGCATTACCCATGATAGCAGAGGCAGAGGATAGGTCTCCAAAGCTATCAAAGTACCGCTGGGTCCAAGGGTAGACAACCAGCAGCCTGCCCAGGGCCTCACCACCAACTTCATCGGAGTTCACCTTTCCCCACAGGCAAGA" [SRR1982462_2]="GGGGTGAAATCCTTGCCCAGGTGGTGGCCCAGCACAATCACGATCATATTGCCCAGGAGCCTGAAGTTCTCAGGATCCACATGCAGCTTGTCACAGTGGAGCTCACTGAGGCTGGCAAAGGTGCCCTTGAGGCTGTCCAAGTGATTCAGGCCATCGTTAAAGGCAGTTATCACCTTCTTGCCATGGGCCTTCACTTTGGCATTACCCATGATAGCAGAGGCAGAGGATAGGTCTCCAAAGCTATCAAAGTACCGCTGGGTCCAAGGGTAGACAACCAGCAGCCTGCCCAGGGCCTCACCACCAACTTCATCGGAGTTCACCTTTCCCCACAGGCAAGA" [SRR1982463_1]="CCGACCCCGAGAGACGCGCACCCGCGGCCGCTCCCGTCCCGTTCCGTTCCGACTGCCGGCGACGGCCGGGTATGGGCCCGACGCTCCAGCGCCATCCATTTTCAGGGCTAGTTGATTCGGCAGGTGAGTTGTTACACACTCCTTAGCGGATTCCGACTTCCATGGCCACCGTCCTGCTGTCTATATCAACCAACACCTTTTCTGGGGTCTGATGAGCGTCGGCATCGGGCGCCTTAACCCGGCGTTCGGTTCATCCCGCAGCGCCAGTTCTGCTTACCAAAAGTGGCCCACTAGGCACTCGCATTCCACGCCCGGCTCCACGCCAGCGAGCCGGGCTTCTTACCCATTTAAAGTTTGAGAATAGGTTGAGATCGTTTCGGCCCCAAGACCTCTAATCATTCGCTTTACCGGATAAAACTGCGTACGTCGGGAGCGAGAGCGCCAGCTATCCTGAGGGAAACTTCGGAGGGAACCAGCTACTAGATGGTTCGATTAGTCTTTCGCCCCTATA" [SRR1982463_2]="CCGACCCCGAGAGACGCGCACCCGCGGCCGCTCCCGTCCCGTTCCGTTCCGACTGCCGGCGACGGCCGGGTATGGGCCCGACGCTCCAGCGCCATCCATTTTCAGGGCTAGTTGATTCGGCAGGTGAGTTGTTACACACTCCTTAGCGGATTCCGACTTCCATGGCCACCGTCCTGCTGTCTATATCAACCAACACCTTTTCTGGGGTCTGATGAGCGTCGGCATCGGGCGCCTTAACCCGGCGTTCGGTTCATCCCGCAGCGCCAGTTCTGCTTACCAAAAGTGGCCCACTAGGCACTCGCATTCCACGCCCGGCTCCACGCCAGCGAGCCGGGCTTCTTACCCATTTAAAGTTTGAGAATAGGTTGAGATCGTTTCGGCCCCAAGACCTCTAATCATTCGCTTTACCGGATAAAACTGCGTACGTCGGGAGCGAGAGCGCCAGCTATCCTGAGGGAAACTTCGGAGGGAACCAGCTACTAGATGGTTCGATTAGTCTTTCGCCCCTATA" [SRR1982464_1]="GGGGTGAAATCCTTGCCCAGGTGGTGGCCCAGCACAATCACGATCATATTGCCCAGGAGCCTGAAGTTCTCAGGATCCACATGCAGCTTGTCACAGTGGAGCTCACTGAGGCTGGCAAAGGTGCCCTTGAGGCTGTCCAAGTGATTCAGGCCATCGTTAAAGGCAGTTATCACCTTCTTGCCATGGGCCTTCACTTTGGCATTACCCATGATAGCAGAGGCAGAGGATAGGTCTCCAAAGCTATCAAAGTACCGCTGGGTCCAAGGGTAGACAACCAGCAGCCTGCCCAGGGCCTCACCACCAACTTCATCGGAGTTCACCTTTCCCCACAGGCAAGA" [SRR1982464_2]="GGGGTGAAATCCTTGCCCAGGTGGTGGCCCAGCACAATCACGATCATATTGCCCAGGAGCCTGAAGTTCTCAGGATCCACATGCAGCTTGTCACAGTGGAGCTCACTGAGGCTGGCAAAGGTGCCCTTGAGGCTGTCCAAGTGATTCAGGCCATCGTTAAAGGCAGTTATCACCTTCTTGCCATGGGCCTTCACTTTGGCATTACCCATGATAGCAGAGGCAGAGGATAGGTCTCCAAAGCTATCAAAGTACCGCTGGGTCCAAGGGTAGACAACCAGCAGCCTGCCCAGGGCCTCACCACCAACTTCATCGGAGTTCACCTTTCCCCACAGGCAAGA" [SRR1982465_1]="GCGGGGGTGAAATC" [SRR1982465_2]="GCCTGTGGGGAAAGGTGAAC" [SRR1982466_1]="GCCTGTGGGGAAAGGTGAAC" [SRR1982466_2]="GCCTGTGGGGAAAGGTGAAC", [SRR1982467_1]="GCGGGGGTGAAATCCTTGCCCAGGTGGTGGCCCAGCACAATCACGATCATATTGCCCAGGAGCCTGAAGTTCTCAGGATCCACATGCAGCTTGTCACAGTGGAGCTCACTGAGGCTGGCAAAGGTGCCCTTGAGGCTGTCCAAGTGATTCAGGCCATCGTTAAAGGC" [SRR1982467_2]="TCGTTAAAGGCAGTTATCACCTT" [SRR1982468_1]="GGGGTGAAATCCTTGCCCAGGTGGTGGCCCAGCACAATCACGATCATATTGCCCAGGAGCCTGAAGTTCTCAGGATCCACATGCAGCTTGTCACAGTGGAGCTCACTGAGGCTGGCAAAGGTGCCCTTGAGGCTGTCCAAGTGATTCAGGCCATCGTTAAAGGCAGTTATCACCTTCTTGCCATGGGCCTTCACTTTGGCATTACCCATGATAGCAGAGGCAGAGGATAGGTCTCCAAAGCTATCAAAGTACCGCTGGGTCCAAGGGTAGACAACCAGCAGCCTGCCCAGGGCCTCACCACCAACTTCATCGGAGTTCACCTTTCCCCACAGGCAAGAGACAGCAGCCTTCTCAGCATCAGTCAGGTGCACCATGATGTCTGTTTCTGGGGTTGTGAGTCA" [SRR1982468_2]="GGGGTGAAATCCTTGCCCAGGTGGTGGCCCAGCACAATCACGATCATATTGCCCAGGAGCCTGAAGTTCTCAGGATCCACATGCAGCTTGTCACAGTGGAGCTCACTGAGGCTGGCAAAGGTGCCCTTGAGGCTGTCCAAGTGATTCAGGCCATCGTTAAAGGCAGTTATCACCTTCTTGCCATGGGCCTTCACTTTGGCATTACCCATGATAGCAGAGGCAGAGGATAGGTCTCCAAAGCTATCAAAGTACCGCTGGGTCCAAGGGTAGACAACCAGCAGCCTGCCCAGGGCCTCACCACCAACTTCATCGGAGTTCACCTTTCCCCACAGGCAAGAGACAGCAGCCTTCTCAGCATCAGTCAGGTGCACCATGATGTCTGTTTCTGGGGTTGTGAGTCA")

fastq_dir="/data/ncbi_data/project5"
output_dir="/home/user04/final_project/fastq_trimmed"

# loop through FastQ files
for fastq_file in "$fastq_dir"/*.fastq.gz; do
    # Extract sample name from the file name (adjust the pattern as needed)
    sample_name=$(basename "$fastq_file" .fastq.gz)

    # Check if the sample name is in the adapter_map
    if [ -n "${adapter_map[$sample_name]}" ]; then
        adapter_sequence="${adapter_map[$sample_name]}"
        output_file="$output_dir/${sample_name}_trimmed.fastq.gz"

        # Use cutadapt tool with the adapter_sequence
        echo "Processing $fastq_file with adapter sequence: $adapter_sequence"

        cutadapt -a "$adapter_sequence" -o "$output_file" "$fastq_file"
    else
        echo "No adapter sequence found for $sample_name"
    fi
done