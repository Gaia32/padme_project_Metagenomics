#!/bin/bash
#SBATCH --job-name=paleo
mkdir -p 7bis_Map_Damage

DirProjet=(/scratch/mahautgr/VFinal)
DirSamt=$DirProjet/7_Samtools/
DirMamD=$DirProjet/7bis_Map_Damage/

for file in ${DirSamt}*_sort_fixmate_MT.sorted.bam; do
    name=${file##*/}
	sample_tag=${name%_sort_fixmate*}

	echo "_______________________ Sample "$sample_tag" _______________________"

    source /local/env/envconda.sh 
    conda activate /groups/Paleogenomics/ENV/MapDamage-2.2.1
    GENOME=/groups/dog/data/canFam3/sequence/bwa_index/Canis_familiaris.CanFam3.1.72.dna_sm.toplevel.fa

    mapDamage -i $DirSamt${sample_tag}_sort_fixmate_MT.sorted.bam -r $GENOME --title=${sample_tag} --no-stats >> MapDamage.pdf

done
