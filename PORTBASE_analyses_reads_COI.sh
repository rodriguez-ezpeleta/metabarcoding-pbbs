#####################################################
#     ANALYSE OF READS FOR COI BARCODE              #
#####################################################

# Files needed:
# PORTBASE_COI.fasta (recovered from PORTBASE_prep_reads_COI.sh)
# PORTBASE_COI.groups (recovered from PORTBASE_prep_reads_COI.sh)
# R Script PREP_SWARM.R 
# BOLD database aligned formatted for Mothur Bold_Database_COI.align
# BOLD database fasta formatted for Mothur COI-BOLD-may2018_good.fasta
# BOLD database taxonomy formatted for Mothur COI-BOLD-may2018_better_v2.tax
# samples_data_COI.txt file containing the metadata of each sample

# This script is for:
# 1- Dereplicate sequences to process faster using Mothur(Schloss et al., 2009)
# 2- Align sequences against BOLD database using Mothur
# 3- Removal of chimeras using Mothur
# 4- OTU clustering using SWARM (Mahé, Rognes, Quince, de Vargas, & Dunthorn, 2014)
# 5- Taxonomic assignement of representative sequences of each OTU using Mothur
# 6- Convert Mothur files into a .biom file using Mothur


### 1- Dereplicate sequences (Retain only unique reads) to process faster using Mothur(Schloss et al., 2009)
mothur "#unique.seqs(fasta=PORTBASE_COI.fasta)"

### 2- Align sequences against databases using Mothur
mothur "#align.seqs(fasta=PORTBASE_COI.unique.fasta, reference=Bold_Database_COI.align, processors=12, flip=T)"
mothur "#summary.seqs(fasta=PORTBASE_COI.unique.align, name=PORTBASE_COI.names, processors=12)"

# Keep sequences aligned in the barcode region and sequences with a minimum lenght threshold of 250bp and no ambiguous base 
mothur "#screen.seqs(fasta=PORTBASE_COI.unique.align, name=PORTBASE_COI.names, group=PORTBASE_COI.groups, minlength=250, start=365, end=635, maxambig=0, processors=12)"
mothur "#summary.seqs(fasta=PORTBASE_COI.unique.good.align, name=PORTBASE_COI.good.names, processors=12)"

# Count the number of sequences per sample
mothur "#count.groups(group=PORTBASE_COI.good.groups)"

# Remove the column where only gaps are present 
mothur "#filter.seqs(fasta=PORTBASE_COI.unique.good.align, vertical=T, processors=8)"

# Dereplicate sequences again 
mothur "#unique.seqs(fasta=PORTBASE_COI.unique.good.filter.fasta, name=PORTBASE_COI.good.names)"
mothur "#summary.seqs(fasta=PORTBASE_COI.unique.good.filter.unique.fasta, name=PORTBASE_COI.unique.good.filter.names, processors=12)"

# Denoise sequences 
mothur "#pre.cluster(fasta=PORTBASE_COI.unique.good.filter.unique.fasta, name=PORTBASE_COI.unique.good.filter.names, group=PORTBASE_COI.good.groups, diffs=2, processors=12)"
mothur "#summary.seqs(fasta=PORTBASE_COI.unique.good.filter.unique.precluster.fasta, name=PORTBASE_COI.unique.good.filter.unique.precluster.names, processors=12)"

### 3- Removal of chimeras using Mothur
# Identify chimera with uchime (Edgar et al., 2011) denovo mode
mothur "#chimera.uchime(fasta=PORTBASE_COI.unique.good.filter.unique.precluster.fasta, name=PORTBASE_COI.unique.good.filter.unique.precluster.names, group=PORTBASE_COI.good.groups, dereplicate=t, processors=12)"

# Removal of chimera
mothur "#remove.seqs(accnos=PORTBASE_COI.unique.good.filter.unique.precluster.denovo.uchime.accnos, fasta=PORTBASE_COI.unique.good.filter.unique.precluster.fasta, name=PORTBASE_COI.unique.good.filter.unique.precluster.names, group=PORTBASE_COI.good.groups, dups=T)"
mothur "#summary.seqs(fasta=PORTBASE_COI.unique.good.filter.unique.precluster.pick.fasta, name=PORTBASE_COI.unique.good.filter.unique.precluster.pick.names, processors=12)"
mothur "#count.groups(group=PORTBASE_COI.good.pick.groups)"

# Change names for easier names
cp PORTBASE_COI.unique.good.filter.unique.precluster.pick.fasta PORTBASE_COI_all.fasta
cp PORTBASE_COI.unique.good.filter.unique.precluster.pick.names PORTBASE_COI_all.names
cp PORTBASE_COI.good.pick.groups PORTBASE_COI_all.groups

### 4- OTU clustering using SWARM ((Mahé, Rognes, Quince, de Vargas, & Dunthorn, 2014)

# Prepare mothur files to the SWARM format
Rscript PREP_SWARM.R PORTBASE_COI_all.names PORTBASE_COI_all.fasta

# Clustering with SWARM
swarm PORTBASE_COI_all_swarm.fasta -d 1 -t 10 -f -l PORTBASE_COI_all_swarm_d1.log -r -o PORTBASE_COI_all_swarm_d1.out -s PORTBASE_COI_all_swarm_d1.stats -w PORTBASE_COI_all_swarm_d1_rep.fasta 

# Prepare SWARM file to Mothur format
awk '/^>/{printf ">Otu%04d\n",++i; next}{print}' PORTBASE_COI_all_swarm_d1_rep.fasta > PORTBASE_COI_all_swarm_d1_rep_otuName.fasta
sed -E  s/\(_[0-9]*,\)/","/g PORTBASE_COI_all_swarm_d1.out   | sed -E s/\(_[0-9]*\\t\)/"\\t"/g | sed -E  s/\(_[0-9]*$\)//g | sed s/swarm/swarm_1/g >  PORTBASE_COI_all_swarm_d1.list

### 5- Taxonomic assignement of representative sequences of each OTU using Mothur
mothur "#classify.seqs(fasta=PORTBASE_COI_all_swarm_d1_rep_otuName.fasta, template=COI-BOLD-may2018_good.fasta, taxonomy=COI-BOLD-may2018_better_v2.tax, method=wang, cutoff=50, processors=12)"

### 6- Convert Mothur files into a .biom file using Mothur
mothur "#make.table(name=PORTBASE_COI_all.names, group=PORTBASE_COI_all.groups, processors=12)"
mothur "#make.shared(list=PORTBASE_COI_all_swarm_d1.list, count=PORTBASE_COI_all.count_table, label=swarm_1)"
awk '{print $1 "\t1\t" $2}' PORTBASE_COI_all_swarm_d1_rep_otuName.COI_BOLD_may2018_better_v2.wang.taxonomy >> PORTBASE_COI_temp1.taxonomy
sed '1 i\OTU\tSize\tTaxonomy' PORTBASE_COI_temp1.taxonomy >> PORTBASE_COI_all_swarm_d1_rep_otuName.COI_BOLD_may2018_better_v2.wang.cons.taxonomy
mothur "#make.biom(shared=PORTBASE_COI_all_swarm_d1.shared, constaxonomy=PORTBASE_COI_all_swarm_d1_rep_otuName.COI_BOLD_may2018_better_v2.wang.cons.taxonomy, metadata=samples_data_COI.txt)"
mv PORTBASE_COI_all_swarm_d1.swarm_1.biom PORTBASE_COI_swarm_rep_otu.biom


# This file is ready to be imported in Phyloseq R package (McMurdie and Holmes (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE. 8(4):e61217) with the function import_biom to get the otu_table
# A phyloseq object has 3 components:
# OTU table
# sample data (metadata for each sample)
# taxonomy data