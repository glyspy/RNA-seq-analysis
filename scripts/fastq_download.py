import subprocess

# samples correspond to Het_1, Het_2, Imm_1, Imm_2
sra_numbers = [
    "SRR1982462", "SRR1982463", "SRR1982464", "SRR1982465", "SRR1982466", "SRR1982467", "SRR1982468"
    ]

# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
#for sra_id in sra_numbers:
#   print ("Currently downloading: " + sra_id)
#    prefetch = "~/Downloads/sratoolkit.3.0.10-mac-arm64/bin/prefetch " + sra_id
#    print ("The command used was: " + prefetch)
#    subprocess.call(prefetch, shell=True)

# this will extract the .sra files from above into a folder named 'fastq'
for sra_id in sra_numbers:
    print ("Generating fastq for: " + sra_id)
    fastq_dump = "/Users/glykeriasp/Downloads/sratoolkit.3.0.10-mac-arm64/bin//fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ~/ncbi/public/sra/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)