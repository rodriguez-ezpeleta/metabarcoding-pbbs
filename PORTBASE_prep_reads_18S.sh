#####################################################
#     PREPARE READS FOR 18S BARCODE                 #
#####################################################

# Primers used are 1389F and 1510R from Amaral-Zettler, McCliment, Ducklow, & Huse, 2009 and target V9 hypervariable region of the 18S rRNA 

# This script is for:
# 1- Check quality of demultiplex paired-end reads with FASTQC (Andrew, 2010) and create a file summarizing the number of reads kept after each step
# 2- Remove primers and adapters with Cutadapt (Martin, 2011)
# 3- Merge forward and reverse reads with PEAR (Zhang, Kobert, Flouri, & Stamatakis, 2014)
# 4- Remove reads with average Phred quality score lower than 20 with Trimmomatic (Bolger, Lohse, & Usadel, 2014)
# 5- Create a fasta file of each sample, merge all .fasta files and create a group file to assign sequences to a specific sample with Mothur (Schloss et al., 2009)


### 1- Check quality of demultiplex paired-end reads with FASTQC (Andrew, 2010) and create a file summarizing the number of reads kept after each step
# File where I the pathway to raw data sequences for 18S
names=(`cut -f 1 samples_18S.txt`)
pathtoR1=(`cut -f 2 samples_18S.txt`)
pathtoR2=(`cut -f 3 samples_18S.txt`)

# Create a folder where we will put all the fastq files to check for quality of demultiplex paired-end reads
mkdir check_quality_fastq

# Create a file where we summarize the numbers of reads kept after each step for each sample with:
## sample = sampleID
## raw = # of reads of raw sequences
## rm_primers = # of reads after primers removal
## merged = # of reads after merging of forward and reverse reads
## quality filtering = # of reads after the removal of low quality reads
echo "sample raw rm_primer merged quality_filtering" > reads_stats.txt  

for i in `seq 0 79`; do
    echo -n ${names[i]} " " >> reads_stats.txt	
    gunzip -c ${pathtoR1[i]} > ${names[i]}_R1.fastq
    gunzip -c ${pathtoR2[i]} > ${names[i]}_R2.fastq
# Check sequences quality with FASTQC 
    fastqc -o check_quality_fastq/ ${names[i]}_R1.fastq
    fastqc -o check_quality_fastq/ ${names[i]}_R2.fastq
# Count # of reads remaining after this step
    raw="$(grep -c '^@M0' ${names[i]}_R2.fastq)"
    echo -n $raw " " >> reads_stats.txt

### 2- Remove primers and adapters with Cutadapt (Martin, 2011)
# Options:
# - Linkedadapter (to remove primer and adapter sequences that flank your sequence): -g for forward (R1) and -G for reverse (R2) 
# - Discard-untrimmed : to throw away all read pairs in which R1 doesn’t start with FWDPRIMER or in which R2 does not start with REVPRIMER.
# - Minimum lenght (-m): discard reads shorter than 50bp 

    cutadapt -g TTGTACACACCGCCC...GTAGGTGAACCTGCRGAAGG -G CCTTCYGCAGGTTCACCTAC...GGGCGGTGTGTACAA --discard-untrimmed -m 50 -o ${names[i]}_R1_cut.fastq -p ${names[i]}_R2_cut.fastq ${names[i]}_R1.fastq ${names[i]}_R2.fastq >> ${names[i]}.log

# Count # of reads remaining after this step
    rm_primer="$(grep -c '^@M0' ${names[i]}_R1_cut.fastq)"
    echo -n $rm_primer " " >> reads_stats.txt

### 3- Merge forward and reverse reads with PEAR (Zhang, Kobert, Flouri, & Stamatakis, 2014)
# Options:
# - minimum overlap lenght = 80bp 
# - maximum possible sequence lenght = 190bp

   pear -f ${names[i]}_R1_cut.fastq -r ${names[i]}_R2_cut.fastq -v 80 -m 190 -o ${names[i]} >> ${names[i]}.log

# Count # of reads remaining after this step

   merged="$(grep -c '^@M0' ${names[i]}.assembled.fastq)"
   echo -n $merged " " >> reads_stats.txt

### 4- Remove reads with average Phred quality score lower than 20 with Trimmomatic (Bolger, Lohse, & Usadel, 2014)
   trimmomatic SE -phred33 -trimlog ${names[i]}_assembled.trimmed.logfile ${names[i]}.assembled.fastq ${names[i]}_assembled.trimmed.fastq AVGQUAL:20 >> ${names[i]}.log

# Count # of reads remaining after this step
   quality="$(grep -c '^@M0' ${names[i]}_assembled.trimmed.fastq)"
   echo  $quality " " >> reads_stats.txt

# Check sequences quality  
   fastqc -o check_quality_fastq/ ${names[i]}_assembled.trimmed.fastq

# 5- Create a fasta file of each sample, merge all .fasta files and create a group file to assign sequences to a specific sample with Mothur (Schloss et al., 2009)
# Split fastq file in the two files fasta and qual 
    mothur "#fastq.info(fastq=${names[i]}_assembled.trimmed.fastq)"

# Merge the fasta files generated for each samples into one fasta file (here called PORTBASE_18S.fasta)
    cat ${names[i]}_assembled.trimmed.fasta >> PORTBASE_18S.fasta

# Preparation of files for make.group
    echo ${names[i]} | tr "\n$" "\-">> group.txt
    echo ${names[i]}_assembled.trimmed.fasta | tr "\n$" "\-">> fasta_file.txt
done

sed -i "s/.$//" fasta_file.txt
sed -i "s/.$//" group.txt

file="fasta_file.txt"
fastas=$(cat "$file")
file="group.txt"
groups=$(cat "$file")

# Create a group file to assign sequences to a specific sample
mothur "#make.group(fasta=$fastas, groups=$groups)"

# Rename the file to PORTBASE_18S.groups
mv mergegroups  PORTBASE_amaral.groups

# The two files that we will need for the next steps of the analyses_18S.sh script are ready:
# Fasta file containing all samples : PORTBASE_18S.fasta
# Group file assigning each sequence to its sample : PORTBASE_18S.groups







