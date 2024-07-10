#!/bin/bash
#SBATCH --job-name=paleo

####################################################
#	This is the first part of the metage pipeline 
# 	which was used to analyse the sedimental DNA
#	from the PADME project. From the adapter 
# 	removal suig cutadapt until the geneous 
#	adapter removal.
####################################################

date
# creating the folders for the temporary files
mkdir -p OUT 1_Cutadapt 2_Rmadapt_Fastp 

DirProjet="/path/to/the/your/project/"
DirRaw="/path/to/the/raw/files/" #Change here the path to the raw files
DirCut=$DirProjet/1_Cutadapt/
DirRmad=$DirProjet/2_Rmadapt_Fastp/
DirOut=$DirProjet/OUT/

######################## ANALYSIS
for file in $DirRaw*_R1.fastq.gz ; do
	name=${file##*/}  # retain the part after the last slash
	number=${name%_*} # keep the sample number before the "_"

	echo ""
	echo "_________________________ Sample "$number" _________________________"

	gunzip -c $DirRaw${number}_R1.fastq.gz > $DirRaw${number}_R1.fq
	gunzip -c $DirRaw${number}_R2.fastq.gz > $DirRaw${number}_R2.fq

# ETAPE 1 - Removal of the adapters on R1 and R2 with CUTADAPT

	echo "___________ 1 - RM ADAPT - CUTADAPT ___________" 
	source /local/env/envcutadapt-1.15.sh
 
	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -O 3 -e 0.5 -m 1 -o $DirCut${sample_tag}_R1_cutadapt.fastq.gz --untrimmed-output=$DirCut${sample_tag}_R1_untrim.fastq.gz $DirRaw${sample_tag}*_R1_001.fastq.gz > $DirOut${sample_tag}_R1_cutadapt_report.txt 2>&1
	cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -O 3 -e 0.5 -m 1 -o $DirCut${sample_tag}_R2_cutadapt.fastq.gz --untrimmed-output=$DirCut${sample_tag}_R2_untrim.fastq.gz $DirRaw${sample_tag}*_R2_001.fastq.gz > $DirOut${sample_tag}_R2_cutadapt_report.txt 2>&1

# ETAPE 2 - Removal of the adapters on R1 and R2 with Fastp

	echo "___________ 2 - RM ADAPT - fastp ___________" 
	fastp_env=/my_env/Fastp/fastp_0_23_1

	$fastp_env -i $DirCut${sample_tag}_R1_cutadapt.fastq.gz -o $DirRmad${sample_tag}_R1_rmadapt.fastq.gz --adapter_fasta "/groups/Paleogenomics/DOG/Adapters_&_SNIP/"adapters.fasta --disable_quality_filtering --disable_length_filtering --disable_trim_poly_g --dont_eval_duplication -h $DirRmad${sample_tag}_R1_rmadapt.html>$DirOut${sample_tag}_R1_rmadapt_report.txt 2>&1
	$fastp_env -i $DirCut${sample_tag}_R2_cutadapt.fastq.gz -o $DirRmad${sample_tag}_R2_rmadapt.fastq.gz --adapter_fasta "/groups/Paleogenomics/DOG/Adapters_&_SNIP/"adapters.fasta --disable_quality_filtering --disable_length_filtering --disable_trim_poly_g --dont_eval_duplication -h $DirRmad${sample_tag}_R2_rmadapt.html>$DirOut${sample_tag}_R2_rmadapt_report.txt 2>&1

done 
date 

# ETAPE 3 - Geneious

# Manually removing residual adapters then mapping against canFan3.