#!/bin/bash
#SBATCH --job-name=paleo

####################################################
#       This is metage pipeline part 2             #
####################################################
date
#DirProjet=(/home/genouest/cnrs_umr6553/mahautgr/DOGS/New_Pipeline) 

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

for file in $DirRaw*_R1.fq ; do
	name=${file##*/}  # retain the part after the last slash
	number=${name%_*} # keep the sample number before the "_"

	echo ""
	echo "_________________________ Sample "$number" _________________________"

echo ""
echo "___________ 2_bis -  GENEIOUS _________________"
# On traite ici la sortie de fastp pour enlever les derniers adaptateurs, ainsi que les untrimmed de cutadapt qui contiennet encore des adaptateurs
# ces fichiers sont placés dans le dir: DirCut_3, nommés sous la forme:
		# 1804_R1_adrestUnpaired.fastq pour les séquences sortie de fastp sans adaptateurs
		# 1804_R1_untrimadaptUnpaired.fastq: les séquences untrimmed de cutadapt extraites avec geneious qui avaient des adaptateurs
		# 1804_R1_unusedUnpaired.fastq: les séquences en sortie de fastp qui avaient des adaptateurs
# on fusionne les séquences sortie de fastp sans adaptateurs (fichier unused), les séquences untrimmed qui avaient des adaptateurs extraites avec geneious (fichier untrimadapt), les séquences en sortie de fastp qui avaient des adaptateurs (fichier adrest)
# ATTENTION : Geneious ajoute parfois des mots à la fin des noms des séquences, on supprime avec awk le mot "extraction" présent dans les séquences extraites
	for file in ${DirGen}*.fastq.gz ; do
		file_name_d=`echo $file | sed -e "s/fastq.gz/fastq/"`
		if [[ ! -f $file_name_d  ]]; then
			gunzip -c $file > $file_name_d
		fi
	done
cat $DirGen${number}_R1_adrestUnpaired.fastq $DirGen${number}_R1_untrimadaptUnpaired.fastq $DirGen${number}_R1_unusedUnpaired.fastq | awk '
	{
		a=gensub(/extraction */,"","g",$0)
		print a
	}
	'> $DirGen${number}_R1_geneious.fq

	cat $DirGen${number}_R2_adrestUnpaired.fastq $DirGen${number}_R2_untrimadaptUnpaired.fastq $DirGen${number}_R2_unusedUnpaired.fastq | awk '
	{
		a=gensub(/extraction */,"","g",$0)
		print a
	}
	'> $DirGen${number}_R2_geneious.fq


echo ""
	echo "___________ 3 -  repair ___________"
source /local/env/envconda.sh
conda activate /home/genouest/cnrs_umr6553/ddiquelou/bbtools
repair.sh in=$DirGen${number}_R1_geneious.fq in2=$DirGen${number}_R2_geneious.fq --overwrite=true out=$DirGen${number}_R1_geneious_rep.fq out2=$DirGen${number}_R2_geneious_rep.fq outs=outfile_singleton.fq

wc -l ${DirGen}* | awk '{print gensub(/.*\//,"","g",$2)";concat_after_geneious;seq_number;" $1/4}'>>${DirOut}report.csv

# ETAPE 3 - Merging des R1 et R2
	echo ""
	echo "___________ 4 - MERGING - FASTP ___________" #laisser qualité, longueur, histogramme de taille des pairés mergés
## OutilL - FASTP
############### Script
(fastp --in1 $DirGen${number}_R1_geneious_rep.fq --in2 $DirGen${number}_R2_geneious_rep.fq -m --merged_out $DirMerg${number}_merged.fq --out1 $DirMerg${number}_R1_unmerged.fq --out2 $DirMerg${number}_R2_unmerged.fq --failed_out $DirMerg${number}_merged_failed.fq --disable_length_filtering --average_qual 20 --trim_poly_g --trim_poly_x --cut_right 20 --disable_adapter_trimming>$DirOut${number}_report_merge.txt 2>&1

nf=${number}_report_merge.txt
cat $DirOut/$nf | sed -n -e "
/reads passed filter:/p
/Read pairs merged:/p"| sed -e "
s/^\(.*\):[ \t]*\([0-9,]*\)[ \t]*$/$nf;fastp_merge;\1;\2/")>>${DirOut}report.csv

cat $DirMerg${number}_merged.fq $DirMerg${number}_R1_unmerged.fq $DirMerg${number}_R2_unmerged.fq>$DirMerg${number}_all_merged.fq

wc -l ${DirMerg}* | awk '
{print gensub(/.*\//,"","g",$2)";concat_after_merge;seq_number;" $1/4}'>>${DirOut}report.csv

	echo ""
	echo "___________ 5 FILTRE DE LONGUEUR ___________"

(fastp -i $DirMerg${number}_all_merged.fq -o $DirMerg${number}_all_merged_len_filter.fq -l 10 --disable_quality_filtering --disable_adapter_trimming -h $DirDupl${number}_all_merged_len_filtere.html>$DirOut${number}_all_merged_len_filter.txt 2>&1

nf=${number}_all_merged_len_filter.txt
cat $DirOut/$nf | sed -n -e "
/reads passed filter:/p
/reads failed due to too short:/p
"| sed -e "
s/^\(.*\):[ \t]*\([0-9,]*\)[ \t]*$/$nf;len_filter;\1;\2/")>>${DirOut}report.csv

# ETAPE 2 - Suppression des duplicats de PCR
	echo ""
	echo "___________ 6 - DDUP - FASTP ___________"

## Outil - FASTP
############### Script
(	fastp -i $DirMerg${number}_all_merged_len_filter.fq --dedup -o $DirDupl${number}_rmduplicate.fq --disable_length_filtering --disable_quality_filtering --disable_adapter_trimming -h $DirDupl${number}_rmduplicate.html>$DirOut${number}_report_rmduplicates.txt 2>&1
nf=${number}_report_rmduplicates.txt
cat $DirOut/$nf | sed -n -e "
/reads passed filter:/p
/reads failed due to too many N:/p
/Duplication rate (may be overestimated since this is SE data):/p " | sed -e "
s/^\(.*\):[ \t]*\([0-9,]*\)[ \t]*$/$nf;rm_duplicates;\1;\2/
s/^\(.*\):[ \t]*\([0-9,\.]*\)%[ \t]*$/$nf;rm_duplicates;\1;;\2/

")>>${DirOut}report.csv



# ETAPE 4 - Alignement sur le génome de référence


done

date