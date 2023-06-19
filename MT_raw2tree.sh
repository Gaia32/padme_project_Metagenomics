#!/bin/bash
#SBATCH --job-name=paleo
# Thx Deborah Diquelou for some parts of her script.

date

##################################################################
#       	 This is dogs pipeline ALL VF              #
##################################################################

#_____________ --> ATTENTION: Il faut changer ici le path du dossier DirProjet (mettre l'endroit où l'on veut faire son analyse)
DirProjet=(/scratch/mahautgr/VFinal) # <------------- + faire attention à ne pas mettre de / à la fin
cd $DirProjet

mkdir -p OUT 0_Raw 1_Cutadapt 2_Rmadapt_FastP 2_bis_Geneious 3_Merging_FastP 4_Filt_FastP 5_Dedup_FastP 6_Mapping 7_Samtools 7bis_Map_Damage 8_Fasta 9_Htsbox 10_Mitotoolpy

DirRaw=(/groups/Paleogenomics/DOG/REGIS_april2022.dir/Nextseq_DEBRUYNE_010422/fastq/) # <----------- Ici le path des données brutes
DirCut=$DirProjet/1_Cutadapt/
DirRmad=$DirProjet/2_Rmadapt_FastP/
DirMerg=$DirProjet/3_Merging_FastP/
DirFilt=$DirProjet/4_Filt_FastP/
DirDdup=$DirProjet/5_Dedup_FastP/
DirMap=$DirProjet/6_Mapping/
DirSam=$DirProjet/7_Samtools/
DirMamD=$DirProjet/7bis_Map_Damage/
DirFast=$DirProjet/8_Fasta/
DirHTS=$DirProjet/9_Htsbox/
DirOut=$DirProjet/OUT/ #Pour mettre les reports de chaque échantillons dans un dossier
Mitotoolpy=/groups/Paleogenomics/ENV/MitoToolsPy/mitotoolpy-seq.py

for file in ${DirRaw}*_R1*.fastq.gz; do 
	name=${file##*/} # on garde tout le nom du fichier après le dernier /
	sample_tag=${name%_R*} # on garde le tag de l'échantillon (avant le _R)

	echo $sample_tag
	echo "_________________________ Sample "$sample_tag" _________________________"

	gunzip -c $DirRaw${sample_tag}"_R1_001.fastq.gz" > $Dir0${sample_tag}_R1.fq
	gunzip -c $DirRaw${sample_tag}"_R2_001.fastq.gz" > $Dir0${sample_tag}_R2.fq


####################################################
# ETAPE 1 - Elimination des adaptateurs sur R1 et R2
## ETAPE 1.1 - CUTADAPT

	echo "___________ 1 - RM ADAPT - CUTADAPT ___________"
############### Environnement CUTADAPT
	source /local/env/envcutadapt-1.15.sh

############### Script
	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -O 3 -e 0.5 -m 1 -o $DirCut${sample_tag}_R1_cutadapt.fq --untrimmed-output=$DirCut${sample_tag}_R1_untrim.fq $Dir0${sample_tag}_R1.fq > $DirOut${sample_tag}_R1_cutadapt_report.txt 2>&1
	cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -O 3 -e 0.5 -m 1 -o $DirCut${sample_tag}_R2_cutadapt.fq --untrimmed-output=$DirCut${sample_tag}_R2_untrim.fq $Dir0${sample_tag}_R2.fq > $DirOut${sample_tag}_R2_cutadapt_report.txt 2>&1


####################################################
## ETAPE 2 - FASTP
	echo "___________ 2 - RM ADAPT - FASTP ___________"
############### Environnement FastP
	source /local/env/envconda.sh
	conda activate /groups/Paleogenomics/ENV/Fastp_0_23/

############### Script
	fastp -i $DirCut${sample_tag}_R1_cutadapt.fq -o $DirRmad${sample_tag}_R1_rmadapt.fq --adapter_fasta "/groups/Paleogenomics/DOG/Adapters_&_SNIP/"adapters.fasta --disable_quality_filtering --disable_length_filtering --dont_eval_duplication -h $DirRmad${sample_tag}_R1_rmadapt.html>$DirOut${sample_tag}_R1_rmadapt_report.txt 2>&1
	fastp -i $DirCut${sample_tag}_R2_cutadapt.fq -o $DirRmad${sample_tag}_R2_rmadapt.fq --adapter_fasta "/groups/Paleogenomics/DOG/Adapters_&_SNIP/"adapters.fasta --disable_quality_filtering --disable_length_filtering --dont_eval_duplication -h $DirRmad${sample_tag}_R2_rmadapt.html>$DirOut${sample_tag}_R2_rmadapt_report.txt 2>&1


####################################################
# ETAPE 3 - Merging et filtration
	echo -e "3"
	echo "___________ 3 - PAIRING - MERGING - BBTOOLS / FASTP ___________"

## OutilL - BBTOOLS
############### Script
	source /local/env/envconda.sh
	conda activate /home/genouest/cnrs_umr6553/ddiquelou/bbtools
    
	echo -e "REPAIR"
    repair.sh in1=$DirRmad${sample_tag}_R1_rmadapt.fq in2=$DirRmad${sample_tag}_R2_rmadapt.fq --overwrite=true out=$DirMerg${sample_tag}_R1_repair.fq out2=$DirMerg${sample_tag}_R2_repair.fq outs=$DirMerg${sample_tag}_repair_singleton.fq>$DirOut${sample_tag}_repair_report.txt 2>&1

## OutilL - FASTP
############### Script
    source /local/env/envconda.sh
	conda activate /groups/Paleogenomics/DOG/ENV/Fastp_0_23/
    
    echo -e "\nMERGE"
    fastp --in1 $DirMerg${sample_tag}_R1_repair.fq --in2 $DirMerg${sample_tag}_R2_repair.fq -A --disable_adapter_trimming --disable_length_filtering --disable_quality_filtering --dont_eval_duplication --overlap_len_require 25 -m --merged_out $DirMerg${sample_tag}_merged.fq --out1 $DirMerg${sample_tag}_R1_unmerged.fq --out2 $DirMerg${sample_tag}_R2_unmerged.fq --unpaired1 $DirMerg${sample_tag}_R1_singleton.fq --unpaired2 $DirMerg${sample_tag}_R2_singleton.fq --failed_out $DirMerg${sample_tag}_merged_failed.fq>$DirOut${sample_tag}_merging_report.txt 2>&1

# on concatène tout, y compris les singletons issus du pairing
	cat $DirMerg${sample_tag}_merged.fq $DirMerg${sample_tag}_R1_unmerged.fq $DirMerg${sample_tag}_R2_unmerged.fq $DirMerg${sample_tag}_R1_singleton.fq $DirMerg${sample_tag}_R2_singleton.fq $DirMerg${sample_tag}_repair_singleton.fq > $DirMerg${sample_tag}_merging_all.fq


####################################################
# ETAPE 4 - Filtration + longueur
	echo "___________ 4 - FILTRATION - FASTP ___________"
    fastp -i $DirMerg${sample_tag}_merging_all.fq -o $DirFilt${sample_tag}_filt.fq -A --disable_adapter_trimming --dont_eval_duplication -l 25 --average_qual 20 --trim_poly_g --trim_poly_x --low_complexity_filter --cut_right --cut_right_mean_quality 20 >$DirOut${sample_tag}_filt_report.txt 2>&1


####################################################
# ETAPE 5 - Suppression des duplicats de PCR
	echo "___________ 5 - DEDUP - FASTP ___________"

## Outil - FASTP
############### Script
	fastp -i $DirFilt${sample_tag}_filt.fq --dedup -o $DirDdup${sample_tag}_rmduplicates.fq -A --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering -h $DirDdup${sample_tag}_rmduplicate.html>$DirOut${sample_tag}_rmduplicates_report.txt 2>&1


####################################################
# ETAPE 6 - Alignement sur le génome de référence
	echo "___________ 6 - BWA ___________"
## Outil - BWA
############### Script

# on met au format fasta
	sed -n '1~4s/^@/>/p;2~4p' $DirDdup${sample_tag}_rmduplicates.fq > $DirDdup${sample_tag}_rmduplicates.fasta

# genome de ref
	GENOME=/groups/dog/data/canFam3/sequence/bwa_index/Canis_familiaris.CanFam3.1.72.dna_sm.toplevel.fa

	source /local/env/envbwa-0.7.17.sh
# ********** MAPPING ALN **********'
	bwa aln -t 8 -l 1024 -o 2 -n 0.01 $GENOME $DirDdup${sample_tag}_rmduplicates.fasta > $DirMap${sample_tag}_aln.sai
# ********** MAPPING SAMSE **********'
	bwa samse $GENOME $DirMap${sample_tag}_aln.sai $DirDdup${sample_tag}_rmduplicates.fq > $DirMap${sample_tag}_aln_samse.sam


####################################################
# ETAPE 7 - Extraction des reads mitochondriaux
## Outil - Samtools
	echo "___________ 7 - SAMTOOLS ___________"

    #en in on prend le fichier .bam des mapped à q25
    #en out on obtient une séquence consensus du génome mitochondrial pour chaque échantillon

    source /local/env/envsamtools-1.15.sh

    # on filtre les unmapped
	samtools view -bT $DirDdup${sample_tag}_rmduplicates.fasta -F 4 $DirMap${sample_tag}_aln_samse.sam > $DirSam${sample_tag}_aln_samse_mapped.sam 
    
    # on filtre les reads d'une qualité < 25
    samtools view -h -q 25 $DirSam${sample_tag}_aln_samse_mapped.sam > $DirSam${sample_tag}_aln_samse_mappedQ25.sam 
    
    # sam to bam mapped Q25
    samtools view -S -b $DirSam${sample_tag}_aln_samse_mappedQ25.sam > $DirSam${sample_tag}_aln_samse_mappedQ25.bam 
    
    #sort par read_name
    samtools sort -n --threads 6 $DirSam${sample_tag}_aln_samse_mappedQ25.bam -o $DirSam${sample_tag}_sort.bam

    #fixmate
    samtools fixmate -m --threads 6 $DirSam${sample_tag}_sort.bam $DirSam${sample_tag}_sort_fixmate.bam

    #sort per coordinates
    samtools sort --threads 6 -o $DirSam${sample_tag}_sort_fixmate.sorted.bam $DirSam${sample_tag}_sort_fixmate.bam

    ##flagstat du markdup
    samtools flagstat $DirSam${sample_tag}_sort_fixmate.sorted.bam > $DirSam${sample_tag}_sort_fixmate_flagstat.sorted.bam ### attention à l'extention là

    ##index le fichier markdup.bam
    samtools index -@ 8 $DirSam${sample_tag}_sort_fixmate.sorted.bam > $DirSam${sample_tag}_sort_fixmate.sorted.bam.bai

    ##extraire les MT
    samtools view $DirSam${sample_tag}_sort_fixmate.sorted.bam MT -o $DirSam${sample_tag}_sort_fixmate_MT.sorted.bam
    samtools view $DirSam${sample_tag}_sort_fixmate_MT.sorted.bam > $DirSam${sample_tag}_sort_fixmate_MT_view.sorted.bam
    
    samtools flagstat $DirSam${sample_tag}_sort_fixmate_MT.sorted.bam > $DirSam${sample_tag}_sort_fixmate_MT_flagstat.sorted.txt

    samtools coverage $DirSam${sample_tag}_sort_fixmate_MT.sorted.bam > $DirOut${sample_tag}_samcoverage.txt
    samtools depth $DirSam${sample_tag}_sort_fixmate_MT.sorted.bam > $DirOut${sample_tag}_samdepth.txt


####################################################
# ETAPE 7 bis - Quantification des motifs de dommages de l'ADN 
## Outil - Mapdamage
    echo "_______________________ Mapdamage _______________________ "
    source /local/env/envconda.sh
    conda activate /groups/Paleogenomics/ENV/MapDamage-2.2.1

    mapDamage -i $DirSamt${sample_tag}_sort_fixmate_MT.sorted.bam -r $GENOME --title=${sample_tag} --no-stats >> MapDamage.pdf

####################################################
# ETAPE 8 - Génération des séquences consensus
## Outil - Htsbox
    echo "_______________________ 8 - Htsbox _________________________"

    source /local/env/envconda.sh
    conda activate /groups/Paleogenomics/ENV/htsbox

    htsbox pileup -cCf $GENOME $DirSam${sample_tag}_sort_fixmate_MT.sorted.bam | less -S 
    htsbox_env pileup -f $GENOME -q20 -q30 -Fs3 $DirSam${sample_tag}_sort_fixmate_MT.sorted.bam > $DirHTS${sample_tag}_htsbox.fasta

    source /local/env/envpython-3.9.5.sh
    python htsbox_correction.py $DirHTS${sample_tag}_htsbox.fasta $DirHTS${sample_tag}_htsbox_corrected.fasta

####################################################
# ETAPE 9 - Génération des halplogroupes
## Outil - Mitotoolpy
    echo "_______________________ 9 - Mitotoolpy _________________________"
    #en in on a la séquence consensus du génome mitochondrial pour chaque échantillon
    #en out on a un tableau avec toutes les données (haplogroupes, variants)

    source /local/env/envpython-3.9.5.sh
    python $Mitotoolpy -s dog -r whole -i $DirHTS${sample_tag}_htsbox_corrected.fasta -o $DirMito${sample_tag}_MitoToolPy.txt

####################################################
# ETAPE 10 - Multialignement
## Outil - Mafft
	echo " _______________________ 10 - Mafft _________________________"

	source /local/env/envmafft-7.475.sh

	mafft --add Data_mito_triees.fasta --reorder Old_alignement_3x.fst > Alignementginsi.fst
	ginsi Alignementginsi.fst > New_alignement_ginsi.fst

####################################################
# ETAPE 10bis - Cleaning de l'alignement
## Outil - Hmmcleaner
	echo " _________________________ 10bis - Hmmcleaner / trimal / prequal_______________________"

    source /local/env/envconda.sh
    conda activate /groups/Paleogenomics/ENV/trimal 
    trimal -in New_alignement_ginsi.fst -out New_alignement_ginsi_curred_spurious_10%.fst -gt 0.3 # <-- on prend celui là

# Options choisies:
# * - gt 0.9 --> si il y a des gap dans 90% des positions, on supprime la colonne
# * -resoverlap 0.60 --> trimAI calcul un score à chaque position par rapport à l'élément le plus commun à hauteur de 60%
# * -seqoverlap 60 --> Si la séquence a 60% de ses positions qui ont un score en dessous du resoverlap, la séquence est supprimée.
# --> il faut que 60% de ses positiosn passent le score pour que la seq soit gardée
# -gt -gapthreshold <n>    1 - (fraction of sequences with a gap allowed).
# trimal -in <inputfile> -out <outputfile> -gt 0.9 -cons 60 
# Removes all positions in the alignment with gaps in 10% or more of the sequences, unless this leaves less than 60%. In such case, print the 60% best (with less gaps) positions. 


####################################################
# ETAPE 11 - Génération des arbres
## Outil - Iqtree
	echo " _________________________ 11 - Iqtree _______________________"

    source /local/env/envconda.sh
    conda activate /groups/Paleogenomics/ENV/iqtree

    -m GTR2+G4+F -alrt 1000

    iqtree -s example.phy -m MFP

# Options choisies:
# * -s --> pour spécifier le nom du fichier d'alignement en input
# * -m MFP --> model finder as a substitution model
# * -B 1000 --> nbr of boostrap replicates for UFboot
# * -b 100 --> nbr of nonparametric boostrap 
# * -m TIM2+I+G --> run UFBoot ultrafastbootstrap
# * -T 8 specifies the nbr og GPU 

done

date
