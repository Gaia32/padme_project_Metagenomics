#!/bin/bash
#SBATCH --job-name=paleo
# Thx Deborah Diquelou for some parts of her script.

####################################################
#       This is metage pipeline part 1             #
####################################################

date
#mkdir OUT 0_Raw 1_Cutadapt 2_Filtration_Fastp 2_bis_Geneious 3_Ddup_FastP 4_Merging_FastP 5_Align_BWA

DirProjet=`pwd`
DirRaw=$DirProjet/ #mettre les fichiers raw dans le répertoire où on exécute le script
DirCut=$DirProjet/1_Cutadapt/
DirFastAdapt=$DirProjet/2_Filtration_Fastp/
DirGen=$DirProjet/2_bis_Geneious/
DirDupl=$DirProjet/3_Ddup_FastP/
DirMerg=$DirProjet/4_Merging_FastP/ 
DirAlign=$DirProjet/5_Align_BWA/
DirOut=$DirProjet/OUT/

mkdir -p $DirRaw $DirCut $DirFastAdapt $DirGen $DirDupl $DirMerg $DirAlign $DirOut

path_fastp=/groups/Paleogenomics/DOG/ENV/Fastp_0_23/bin
export PATH=$PATH:$path_fastp

#entête du fichier csv
echo "file;tool;report;amount;percent"> ${DirOut}report.csv

#début de l'analyse
for file in $DirRaw*_R1.fastq.gz ; do
	name=${file##*/}  # retain the part after the last slash
	number=${name%_*} # keep the sample number before the "_"

	echo ""
	echo "_________________________ Sample "$number" _________________________"

	gunzip -c $DirRaw${number}_R1.fastq.gz > $DirRaw${number}_R1.fq
	gunzip -c $DirRaw${number}_R2.fastq.gz > $DirRaw${number}_R2.fq

# ETAPE 1 - Elimination des adaptateurs sur R1 et R2
## ETAPE 1.1 - CUTADAPT

	echo ""
	echo "___________ 1 - RMADAPT - CUTADAPT ___________"
############### Environnement CUTADAPT
	source /local/env/envcutadapt-1.15.sh

############### Script

(	cutadapt $DirRaw${number}_R1.fq -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -O 3 -e 0.5 -o $DirCut${number}_R1_cutadapt.fq --untrimmed-output=$DirCut${number}_R1_cutadapt_untrimm.fq > $DirOut${number}_R1_report_cutadapt.txt 2>&1
	nf=${number}_R1_report_cutadapt.txt
	cat $DirOut/$nf | sed -n -e "
/Total reads processed:/p
/Reads with adapters:/p
/Total basepairs processed:/p
/Total written (filtered):/p " | sed -e "
s/^\(.*\):[ \t]*\([0-9,]*\)[ \t]*\(.*\)$/$nf;cutadapt;\1;\2;\3/" | sed -e "
s/,//g
s/; *bp */;/
s/;(/;/
s/)$//
s/%$//
"
	cutadapt $DirRaw${number}_R2.fq -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -O 3 -e 0.5 -o $DirCut${number}_R2_cutadapt.fq --untrimmed-output=$DirCut${number}_R2_cutadapt_untrimm.fq > $DirOut${number}_R2_report_cutadapt.txt 2>&1

	nf2=${number}_R2_report_cutadapt.txt
	cat $DirOut/$nf2 | sed -n -e "
/Total reads processed:/p
/Reads with adapters:/p
/Total basepairs processed:/p
/Total written (filtered):/p " | sed -e "
s/^\(.*\):[ \t]*\([0-9,]*\)[ \t]*\(.*\)$/$nf2;cutadapt;\1;\2;\3/" | sed -e "
s/,//g
s/; *bp */;/
s/;(/;/
s/)$//
s/%$//
")>>${DirOut}report.csv

## ETAPE 1.2 - FAST
	echo ""
	echo "___________ 2 - RMADAPT - FASTP ___________"

############### Script
(	fastp -i $DirCut${number}_R1_cutadapt.fq --adapter_fasta "/groups/Paleogenomics/DOG/Adapters_&_SNIP/"adapters.fasta --disable_quality_filtering --disable_length_filtering -o $DirFastAdapt${number}_R1_wtout_qul_len.fq -h $DirOut${number}_R1_wtout_qul_len.html > $DirOut${number}_R1_report_fastp_adapt_wtout_qul_len.txt 2>&1

	nf=${number}_R1_report_fastp_adapt_wtout_qul_len.txt
	cat $DirOut$nf | sed -n -e "

/reads passed filter:/p
/reads failed due to too many N: /p
/bases trimmed due to adapters:/p
/reads with adapter trimmed: /p
 " | sed -e "
s/^\(.*\):[ \t]*\([0-9,]*\)[ \t]*$/$nf;fastp_adapters;\1;\2/"

	fastp -i $DirCut${number}_R2_cutadapt.fq --adapter_fasta "/groups/Paleogenomics/DOG/Adapters_&_SNIP/"adapters.fasta --disable_quality_filtering --disable_length_filtering -o $DirFastAdapt${number}_R2_wtout_qul_len.fq -h $DirOut${number}_R2_wtout_qul_len.html > $DirOut${number}_R2_report_fastp_adapt_wtout_qul_len.txt 2>&1
#enlever les filtres l et q, pdt l'étape de merging, voir le % de séquences enlevées en fonction de juste adapt remove
nf2=${number}_R2_report_fastp_adapt_wtout_qul_len.txt
	cat $DirOut$nf2 | sed -n -e "

/reads passed filter:/p
/reads failed due to too many N: /p
/bases trimmed due to adapters:/p
/reads with adapter trimmed: /p
" | sed -e "
s/^\(.*\):[ \t]*\([0-9,]*\)[ \t]*$/$nf2;fastp_adapters;\1;\2/"
)>>${DirOut}report.csv
done 
date 
