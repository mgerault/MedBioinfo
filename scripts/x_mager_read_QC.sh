#!/bin/bash
echo "script start: download and initial sequencing read quality control"
date

echo "export run accession number from assigned patients"
sqlite3 -batch /proj/applied_bioinformatics/common_data/sample_collab.db "SELECT run_accession FROM sample_annot spl LEFT JOIN sample2bioinformatician s2b using(patien>

echo "create subdirectory folder to store raw FASTQ files and download fastq files from assigned accessions"
mkdir ./data/sra_fastq/
cat ./analyses/x_mager_run_accessions.txt | srun --cpus-per-task=1 --time=00:30:00 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif xargs fastq-dump --split-files --gzip --readids --outdir ./data/sra_fastq/ --disable-multithreading 

# count number of reads per FASTQ file
# ls -d $PWD/data/sra_fastq/* | xargs wc -l
# base call quality score are encoded in ASCII character: in increasing order of quality: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
# print statistics from each FASTQ files with seqkit
srun --cpus-per-task=1 --time=00:30:00 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif seqkit stats ./data/sra_fastq/*.fastq.gz --threads 1

echo "check if FASTQ files have been de-replicated and remove duplicates from all files at once"
ls -d $PWD/data/sra_fastq/* | xargs -I % sh -c 'zcat % | srun --cpus-per-task=1 --time=00:30:00 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif seqkit rmdup -s -o "${1/.fastq.gz/}_clean.${1#*.}" --threads 1' -- %


echo "check if FASTQ files already trimmed"
#  (on de-replicated FASTQ obtained before): full adapter + reversed one (33 base)
ls -d $PWD/data/sra_fastq/* | grep clean | xargs -I % sh -c 'echo %; zcat % | srun --cpus-per-task=1 --time=00:30:00 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif seqkit locate -p "^AGATCGGAAGAGCACACGTCTGAACTCCAGTCA|^AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|ACTGACCTCAAGTCTGCACACGAGAAGGCTAGA$|TGTGAGAAAGGGATGTGCTGCGAGAAGGCTAGA$" -r --threads 1'
# on shoter adapter: 20, 10, 5
ls -d $PWD/data/sra_fastq/* | grep clean | xargs -I % sh -c 'echo %; zcat % | srun --cpus-per-task=1 --time=00:30:00 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif seqkit locate -p "^AGATCGGAAGAGCACACGTC|^AGATCGGAAGAGCGTCGTGT|ACTGACCTCAAGTCTGCACA$|TGTGAGAAAGGGATGTGCTG$" -r --threads 1'
ls -d $PWD/data/sra_fastq/* | grep clean | xargs -I % sh -c 'echo %; zcat % | srun --cpus-per-task=1 --time=00:30:00 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif seqkit locate -p "^AGATCGGAAG|^AGATCGGAAG|ACTGACCTCA$|TGTGAGAAAG$" -r --threads 1'
ls -d $PWD/data/sra_fastq/* | grep clean | xargs -I % sh -c 'echo %; zcat % | srun --cpus-per-task=1 --time=00:30:00 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif seqkit locate -p "^AGATC|^AGATC|ACTGA$|TGTGA$" -r --threads 1'


echo "Quality control + merging paire end reads"
### Quality control
mkdir ./analyses/fastqc
srun --cpus-per-task=2 --time=00:30:00 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif xargs -I{} -a ./analyses/x_mager_run_accessions.txt fastqc ./data/sra_fastq/{}_1.fastq.gz  ./data/sra_fastq/{}_2.fastq.gz -o ./analyses/fastqc --noextract -t 2

### Merging paired end reads
mkdir ./data/merged_pairs
srun --cpus-per-task=2 --time=00:30:00 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif xargs -n 1 -a ./analyses/x_mager_run_accessions.txt -I{} flash -t 2 -d ./data/merged_pairs -z -o {}.flash -M 150 ./data/sra_fastq/{}_1.fastq.gz ./data/sra_fastq/{}_2.fastq.gz 2>&1 | tee -a ./analyses/x_mager_flash.log


echo "Searching for alignment with potential contaminants"
### Use read mapping to check for PhiX contamination 
## doxnload Phix genome
mkdir ./data/reference_seqs
singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif efetch -db nuccore -id NC_001422 -format fasta > ./data/reference_seqs/PhiX_NC_001422.fna

## Create bowtie2 indexed database
mkdir ./data/bowtie2_DBs
srun --cpus-per-task=1 --time=00:10:00 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif bowtie2-build -f ./data/reference_seqs/PhiX_NC_001422.fna ./data/bowtie2_DBs/PhiX_bowtie2_DB

## Run bowtie2
mkdir ./analyses/bowtie
srun --cpus-per-task=8 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif bowtie2 -x ./data/bowtie2_DBs/PhiX_bowtie2_DB -U ./data/merged_pairs/ERR*.extendedFrags.fastq.gz -S ./analyses/bowtie/x_mager_merged2PhiX.sam --threads 8 --no-unal 2>&1 | tee ./analyses/bowtie/x_mager_bowtie_merged2PhiX.log

### same procedure for SARS-CoV-2
singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif efetch -db nuccore -id NC_045512 -format fasta > ./data/reference_seqs/PhiX_NC_045512.fna
srun --cpus-per-task=1 --time=00:10:00 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif bowtie2-build -f ./data/reference_seqs/PhiX_NC_045512.fna ./data/bowtie2_DBs/SC2_bowtie2_DB
srun --cpus-per-task=8 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif bowtie2 -x ./data/bowtie2_DBs/SC2_bowtie2_DB -U ./data/merged_pairs/ERR*.extendedFrags.fastq.gz -S ./analyses/bowtie/x_mager_merged2SC2.sam --threads 8 --no-unal 2>&1 | tee ./analyses/bowtie/x_mager_bowtie_merged2SC2.log

### same oral bacterias: Streptococcus oralis
singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif efetch -db nuccore -id NZ_KQ970792 -format fasta > ./data/reference_seqs/PhiX_NZ_KQ970792.fna
srun --cpus-per-task=1 --time=00:10:00 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif bowtie2-build -f ./data/reference_seqs/PhiX_NZ_KQ970792.fna ./data/bowtie2_DBs/StrO_bowtie2_DB
srun --cpus-per-task=8 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif bowtie2 -x ./data/bowtie2_DBs/StrO_bowtie2_DB -U ./data/merged_pairs/ERR*.extendedFrags.fastq.gz -S ./analyses/bowtie/x_mager_merged2StrO.sam --threads 8 --no-unal 2>&1 | tee ./analyses/bowtie/x_mager_bowtie_merged2StrO.log

### same oral bacterias: Streptococcus salivarius  --> 0.44% were aligned with this strain !
singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif efetch -db nuccore -id NZ_CP040804 -format fasta > ./data/reference_seqs/PhiX_NZ_CP040804.fna
srun --cpus-per-task=1 --time=00:10:00 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif bowtie2-build -f ./data/reference_seqs/PhiX_NZ_CP040804.fna ./data/bowtie2_DBs/StrS_bowtie2_DB
srun --cpus-per-task=8 singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif bowtie2 -x ./data/bowtie2_DBs/StrS_bowtie2_DB -U ./data/merged_pairs/ERR*.extendedFrags.fastq.gz -S ./analyses/bowtie/x_mager_merged2StrS.sam --threads 8 --no-unal 2>&1 | tee ./analyses/bowtie/x_mager_bowtie_merged2StrS.log


echo "Using multiqc to report results"
### Combine quality control results into one unique report for all samples analysed
srun singularity exec /proj/applied_bioinformatics/users/x_mager/FASTQ_assignment/FASTQ_assign_image.sif multiqc --force --title "x_mager sample sub-set" ./data/merged_pairs/ ./analyses/fastqc/ ./analyses/x_mager_flash.log ./analyses/bowtie/

date
echo "script end."



















