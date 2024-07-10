#!/bin/bash
#SBATCH --job-name=premap_p2

####################################################
#       Pipeline premapping part 2             #
####################################################
date

DirProjet=/groups/Paleogenomics/Copro/00_Raw
DirRaw=$DirProjet/année3-CROC-Novaseq_DEBRUYNE_041122_metaG/fastq/
#DirCut=$DirProjet/Cutadapt/
DirFastAdapt=$DirProjet/Fastp_adapters/to_do/
DirGen=$DirProjet/Geneious/
DirRep=$DirProjet/Repair/copro/
DirMerg=$DirProjet/Merging_FastP/copro/ 
DirDupl=$DirProjet/Ddup_FastP/
DirLength=$DirProjet/Length_filter
DirLength_merged=$DirLength/Merged/
DirLength_unmerged=$DirLength/Unmerged/
DirLength_singletons=$DirLength/Singletons/
DirAlign=$DirProjet/Align_BWA/
DirOut=$DirProjet/OUT/

mkdir -p $DirRep $DirMerg $DirDupl $DirLength $DirLength_merged $DirLength_unmerged $DirLength_singletons $DirAlign $DirOut

path_fastp=/groups/Paleogenomics/ENV/Fastp_0_23/bin
export PATH=$PATH:$path_fastp

#entete du fichier csv
echo "file;tool;report;amount;percent"> ${DirOut}report_essai_part2.csv


# Décompression éventuelle des fichiers gz

#ls -c1 ${DirRaw}*.fastq.gz 2>/dev/null | while read file ; do 
#echo $file trouvé
	#file_name_d=`echo $file | sed -e "s/fastq.gz/fq/"`
	#if [[ ! -f $file_name_d ]]; then
#echo decompression $file en $file_name_d
	#	gunzip -c $file>$file_name_d
	#fi
#done

for file in $DirRaw*_R1_001.fastq.gz ; do
	name=${file##*/}  # retain the part after the last slash
	number=${name%_R1*} # keep the sample number before the "_"

	echo ""
	echo "_________ Sample "$number" _________"

	echo ""
	#echo "____  GENEIOUS ______"
	# On traite ici la sortie de fastp pour enlever les derniers adaptateurs, ainsi que les untrimmed de cutadapt qui contiennet encore des adaptateurs
	# ces fichiers sont places dans le dir: DirCut_3, nommes sous la forme:
			# 1804_R1_adrestUnpaired.fastq pour les sequences sortie de fastp sans adaptateurs
			# 1804_R1_untrimadaptUnpaired.fastq: les sequences untrimmed de cutadapt extraites avec geneious qui avaient des adaptateurs
			# 1804_R1_unusedUnpaired.fastq: les sÃ©quences en sortie de fastp qui avaient des adaptateurs
	# on fusionne les sequences sortie de fastp sans adaptateurs (fichier unused), les sÃ©quences untrimmed qui avaient des adaptateurs extraites avec geneious (fichier untrimadapt), les sÃ©quences en sortie de fastp qui avaient des adaptateurs (fichier adrest)
	# ATTENTION : Geneious ajoute parfois des mots Ã  la fin des noms des sequences, on supprime avec awk le mot "extraction" prÃ©sent dans les sequences extraites
	
	#cat $DirGen${number}_R1_postgeneious.fastq $DirGen${number}_R1_untrimadaptUnpaired.fastq $DirGen${number}_R1_unusedUnpaired.fastq | awk '
		#{
		#	a=gensub(/extraction */,"","g",$0)
		#	print a
		#}
		#'> $DirGen${number}_R1_geneious.fq

	#cat $DirGen${number}_R2_postgeneious.fastq $DirGen${number}_R2_untrimadaptUnpaired.fastq $DirGen${number}_R2_unusedUnpaired.fastq | awk '
		#{
		#	a=gensub(/extraction */,"","g",$0)
		#	print a
		#}
		#'> $DirGen${number}_R2_geneious.fq

	echo ""
	echo "____ 3 -  REPAIR ____"
	source /local/env/envconda.sh
	#conda activate /home/genouest/cnrs_umr6553/ddiquelou/bbtools
	conda activate /groups/Paleogenomics/ENV/bbtools
	repair.sh in=$DirFastAdapt${number}_R1_wtout_qul_len.fq in2=$DirFastAdapt${number}_R2_wtout_qul_len.fq --overwrite=true out=$DirRep${number}_R1_rep.fq out2=$DirRep${number}_R2_rep.fq outs=$DirRep${number}_outfile_singleton.fq >${DirOut}${number}_report_repair.txt 2>&1

	wc -l ${DirRep}* | awk '{print gensub(/.*\//,"","g",$2)";rep_after_fastP;seq_number;" $1/4}'>>${DirOut}report_essai_part2.csv


	echo ""
	echo "____ 4 - MERGING - FASTP ____"
	## Outil - FASTP
	############### Script
	(fastp --in1 $DirRep${number}_R1_rep.fq --in2 $DirRep${number}_R2_rep.fq -m --merged_out $DirMerg${number}_merged.fq --out1 $DirMerg${number}_R1_unmerged.fq --out2 $DirMerg${number}_R2_unmerged.fq --failed_out $DirMerg${number}_merged_failed.fq --disable_length_filtering --average_qual 20 --trim_poly_g --trim_poly_x --cut_right --cut_right_mean_quality 20 --disable_adapter_trimming --overlap_len_require 5 >$DirOut${number}_report_merge.txt 2>&1
		nf=${number}_report_merge.txt
		cat $DirOut/$nf | sed -n -e "
		/reads passed filter:/p
		/Read pairs merged:/p"| sed -e "
		s/^\(.\):[ \t]\([0-9,]\)[ \t]$/$nf;fastp_merge;\1;\2/"
	)>>${DirOut}report_essai_part2.csv


	#
	#echo ""
	#echo "____ 5 - QUALITE SINGLETONS - FASTP ____"

	(
		fastp --in1 $DirRep${number}_outfile_singleton.fq --out1 $DirRep${number}_outfile_singleton_qual_filter.fq --disable_length_filtering --average_qual 20 --trim_poly_g --trim_poly_x --cut_right 20 --disable_adapter_trimming > ${DirOut}${number}_outfile_singleton.txt 2>&1
		nf=${number}_outfile_singleton.txt
		cat $DirOut/$nf | sed -n -e "
		/reads passed filter:/p
		"| sed -e "
		s/^\(.\):[ \t]\([0-9,]\)[ \t]$/$nf;filtre_qualite_singletons;\1;\2/"
	)>>${DirOut}report_essai_part2.csv

# concat unmerged R1 et R2
	(
		cat $DirMerg${number}_R1_unmerged.fq $DirMerg${number}_R2_unmerged.fq >$DirMerg${number}_all_unmerged.fq

		wc -l ${DirMerg}* | awk '
		{print gensub(/.*\//,"","g",$2)";concat_unM_after_merge;seq_number;" $1/4}'
	)>>${DirOut}report_essai_part2.csv



# ETAPE 2 - Suppression des duplicats de PCR
	echo ""
	echo "____ 6 - DDUP - FASTP ____"

	## Outil - FASTP
	############### Script
	(	
		fastp -i $DirMerg${number}_merged.fq --dedup -o $DirDupl${number}_merged_rmduplicate.fq --disable_length_filtering --disable_quality_filtering --disable_adapter_trimming -h $DirDupl${number}_merged_rmduplicate.html>$DirOut${number}_report_merged_rmduplicates.txt 2>&1
		nf=${number}_report_merged_rmduplicates.txt
		cat $DirOut/$nf | sed -n -e "
		/merged reads passed filter:/p
		/merged reads failed due to too many N:/p
		/merged Duplication rate (may be overestimated since this is SE data):/p " | sed -e "
		s/^\(.\):[ \t]\([0-9,]\)[ \t]$/$nf;rm_duplicates;\1;\2/
		s/^\(.\):[ \t]\([0-9,\.]\)%[ \t]$/$nf;rm_duplicates;\1;;\2/
		"
	)>>${DirOut}report_essai_part2.csv


(	
		fastp -i $DirMerg${number}_all_unmerged.fq --dedup -o $DirDupl${number}_unmerged_rmduplicate.fq --disable_length_filtering --disable_quality_filtering --disable_adapter_trimming -h $DirDupl${number}_unmerged_rmduplicate.html>$DirOut${number}_report_unmerged_rmduplicates.txt 2>&1
		nf=${number}_report_unmerged_rmduplicates.txt
		cat $DirOut/$nf | sed -n -e "
		/unmerged reads passed filter:/p
		/unmerged reads failed due to too many N:/p
		/unmerged Duplication rate (may be overestimated since this is SE data):/p " | sed -e "
		s/^\(.\):[ \t]\([0-9,]\)[ \t]$/$nf;rm_duplicates;\1;\2/
		s/^\(.\):[ \t]\([0-9,\.]\)%[ \t]$/$nf;rm_duplicates;\1;;\2/
		"
	)>>${DirOut}report_essai_part2.csv

(	
		fastp -i $DirRep${number}_outfile_singleton_qual_filter.fq --dedup -o $DirDupl${number}_singleton_rmduplicate.fq --disable_length_filtering --disable_quality_filtering --disable_adapter_trimming -h $DirDupl${number}_singleton_rmduplicate.html>$DirOut${number}_report_singleton_rmduplicates.txt 2>&1
		nf=${number}_report_singleton_rmduplicates.txt
		cat $DirOut/$nf | sed -n -e "
		/singleton reads passed filter:/p
		/singleton reads failed due to too many N:/p
		/singleton Duplication rate (may be overestimated since this is SE data):/p " | sed -e "
		s/^\(.\):[ \t]\([0-9,]\)[ \t]$/$nf;rm_duplicates;\1;\2/
		s/^\(.\):[ \t]\([0-9,\.]\)%[ \t]$/$nf;rm_duplicates;\1;;\2/
		"
	)>>${DirOut}report_essai_part2.csv

	echo ""
	echo "____ 7 - FILTRE DE LONGUEUR ____"

	(fastp -i $DirDupl${number}_merged_rmduplicate.fq -o $DirLength_merged${number}_merged_rmduplicate_sup35bp.fq -l 34 --disable_quality_filtering --disable_adapter_trimming -h $DirLength_merged${number}_merged_rmduplicate_sup35bp.html>$DirOut${number}_merged_rmduplicate_sup35bp.txt 2>&1

		nf=${number}_merged_rmduplicate_sup35bp.txt
		cat $DirOut/$nf | sed -n -e "
		/merged reads passed filter:/p
		/merged reads failed due to too short:/p
		"| sed -e "
		s/^\(.\):[ \t]\([0-9,]\)[ \t]$/$nf;len_filter;\1;\2/")>>${DirOut}report_essai_part2.csv

(fastp -i $DirDupl${number}_unmerged_rmduplicate.fq -o $DirLength_unmerged${number}_unmerged_rmduplicate_sup35bp.fq -l 34 --disable_quality_filtering --disable_adapter_trimming -h $DirLength_unmerged${number}_unmerged_rmduplicate_sup35bp.html>$DirOut${number}_unmerged_rmduplicate_sup35bp.txt 2>&1

		nf=${number}_unmerged_rmduplicate_sup35bp.txt
		cat $DirOut/$nf | sed -n -e "
		/unmerged reads passed filter:/p
		/unmerged reads failed due to too short:/p
		"| sed -e "
		s/^\(.\):[ \t]\([0-9,]\)[ \t]$/$nf;len_filter;\1;\2/")>>${DirOut}report_essai_part2.csv

(fastp -i $DirDupl${number}_singleton_rmduplicate.fq -o $DirLength_singletons${number}_singleton_rmduplicate_sup35bp.fq -l 34 --disable_quality_filtering --disable_adapter_trimming -h $DirLength_singletons${number}_singleton_rmduplicate_sup35bp.html>$DirOut${number}_singleton_rmduplicate_sup35bp.txt 2>&1

		nf=${number}_singleton_rmduplicate_sup35bp.txt
		cat $DirOut/$nf | sed -n -e "
		/unmerged reads passed filter:/p
		/unmerged reads failed due to too short:/p
		"| sed -e "
		s/^\(.\):[ \t]\([0-9,]\)[ \t]$/$nf;len_filter;\1;\2/")>>${DirOut}report_essai_part2.csv

	
done

date