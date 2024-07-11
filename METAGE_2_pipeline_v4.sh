#!/bin/bash
#SBATCH --job-name=metage2

####################################################
#       Pipeline premapping part 2             #
####################################################

date

mkdir -p 4_Repair 5_Merging_FastP 6_Ddup_FastP 7_Length_filter 8_Merged 8_Unmerged 9_Singletons 10_Align_BWA

DirProjet="/path/to/the/your/project/"
DirRaw="/path/to/the/raw/files/" #Change here the path to the raw files
DirGen=$DirProjet/3_Geneious/
DirRep=$DirProjet/4_Repair/
DirMerg=$DirProjet/5_Merging_FastP/
DirLength_merged=$DirLength/6_Merged/
DirLength_unmerged=$DirLength/6_Unmerged/
DirLength_singletons=$DirLength/6_Singletons/
DirDupl=$DirProjet/7_Ddup_FastP/
DirLength=$DirProjet/8_Length_filter/
DirOut=$DirProjet/OUT/

for file in $DirRaw*_R1_001.fastq.gz ; do
	name=${file##*/}  # retain the part after the last slash
	number=${name%_R1*} # keep the sample number before the "_"

	echo ""
	echo "_________ Sample "$number" _________"

	echo ""
	echo "____  GENEIOUS ______"
	# On traite ici la sortie de l'étape 2 - Fastp pour enlever les derniers adaptateurs, ainsi que les untrimmed de cutadapt qui contiennet encore des adaptateurs
	# ces fichiers sont nommes sous la forme:
			# 1804_R1_adrestUnpaired.fastq pour les sequences sortie de fastp sans adaptateurs
			# 1804_R1_untrimadaptUnpaired.fastq: les sequences untrimmed de cutadapt extraites avec geneious qui avaient des adaptateurs
			# 1804_R1_unusedUnpaired.fastq: les sÃ©quences en sortie de fastp qui avaient des adaptateurs
	# on fusionne les sequences sortie de fastp sans adaptateurs (fichier unused), les séquences untrimmed qui avaient des adaptateurs extraites avec geneious (fichier untrimadapt), les séquences en sortie de fastp qui avaient des adaptateurs (fichier adrest)
	# ATTENTION : Geneious ajoute parfois des mots Ã  la fin des noms des sequences, on supprime avec awk le mot "extraction" prÃ©sent dans les sequences extraites
	
	#____Voici la commande awk:
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
	echo "____ 4 -  REPAIR ____"
	source /local/env/envconda.sh
	conda activate /ENV/bbtools
	repair.sh in=$DirFastAdapt${number}_R1_wtout_qul_len.fq in2=$DirFastAdapt${number}_R2_wtout_qul_len.fq --overwrite=true out=$DirRep${number}_R1_rep.fq out2=$DirRep${number}_R2_rep.fq outs=$DirRep${number}_outfile_singleton.fq >${DirOut}${number}_report_repair.txt 2>&1

	echo ""
	echo "____ 5 - MERGING - FASTP ____"
	fastp_env=/my_env/Fastp/fastp_0_23_1

	fastp --in1 $DirRep${number}_R1_rep.fq --in2 $DirRep${number}_R2_rep.fq -m --merged_out $DirMerg${number}_merged.fq --out1 $DirMerg${number}_R1_unmerged.fq --out2 $DirMerg${number}_R2_unmerged.fq --failed_out $DirMerg${number}_merged_failed.fq --disable_length_filtering --average_qual 20 --trim_poly_g --trim_poly_x --cut_right --cut_right_mean_quality 20 --disable_adapter_trimming --overlap_len_require 5 >$DirOut${number}_report_merge.txt 2>&1

	echo ""
	echo "____ 6 - QUALITE SINGLETONS - FASTP ____"

	fastp --in1 $DirRep${number}_outfile_singleton.fq --out1 $DirRep${number}_outfile_singleton_qual_filter.fq --disable_length_filtering --average_qual 20 --trim_poly_g --trim_poly_x --cut_right 20 --disable_adapter_trimming > ${DirOut}${number}_outfile_singleton.txt 2>&1

	#_____ concat unmerged R1 et R2
		cat $DirMerg${number}_R1_unmerged.fq $DirMerg${number}_R2_unmerged.fq >$DirMerg${number}_all_unmerged.fq


	echo ""
	echo "____ 7 - DDUP - FASTP ____"

	fastp -i $DirMerg${number}_merged.fq --dedup -o $DirDupl${number}_merged_rmduplicate.fq --disable_length_filtering --disable_quality_filtering --disable_adapter_trimming -h $DirDupl${number}_merged_rmduplicate.html>$DirOut${number}_report_merged_rmduplicates.txt 2>&1
	
	fastp -i $DirMerg${number}_all_unmerged.fq --dedup -o $DirDupl${number}_unmerged_rmduplicate.fq --disable_length_filtering --disable_quality_filtering --disable_adapter_trimming -h $DirDupl${number}_unmerged_rmduplicate.html>$DirOut${number}_report_unmerged_rmduplicates.txt 2>&1

	fastp -i $DirRep${number}_outfile_singleton_qual_filter.fq --dedup -o $DirDupl${number}_singleton_rmduplicate.fq --disable_length_filtering --disable_quality_filtering --disable_adapter_trimming -h $DirDupl${number}_singleton_rmduplicate.html>$DirOut${number}_report_singleton_rmduplicates.txt 2>&1


	echo ""
	echo "____ 8 - FILTRE DE LONGUEUR ____"

	fastp -i $DirDupl${number}_merged_rmduplicate.fq -o $DirLength_merged${number}_merged_rmduplicate_sup35bp.fq -l 34 --disable_quality_filtering --disable_adapter_trimming -h $DirLength_merged${number}_merged_rmduplicate_sup35bp.html>$DirOut${number}_merged_rmduplicate_sup35bp.txt 2>&1

	fastp -i $DirDupl${number}_unmerged_rmduplicate.fq -o $DirLength_unmerged${number}_unmerged_rmduplicate_sup35bp.fq -l 34 --disable_quality_filtering --disable_adapter_trimming -h $DirLength_unmerged${number}_unmerged_rmduplicate_sup35bp.html>$DirOut${number}_unmerged_rmduplicate_sup35bp.txt 2>&1

	fastp -i $DirDupl${number}_singleton_rmduplicate.fq -o $DirLength_singletons${number}_singleton_rmduplicate_sup35bp.fq -l 34 --disable_quality_filtering --disable_adapter_trimming -h $DirLength_singletons${number}_singleton_rmduplicate_sup35bp.html>$DirOut${number}_singleton_rmduplicate_sup35bp.txt 2>&1

done

date