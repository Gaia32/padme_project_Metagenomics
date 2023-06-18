#!/bin/bash
#SBATCH --job-name=htsbox

# This script allows to use htsbox to built consensus fasta from bam files
# /!\ be carefull, the parametrisation might not be good, you might want to use geneous or angsd to build consensus

Dir=(PATH/TO/THE/PROJECT) 

DirFast=$DirProjet/8_Fasta/
DirHTS=$DirProjet/9_Htsbox/
GENOME=/groups/dog/data/canFam3/sequence/bwa_index/Canis_familiaris.CanFam3.1.72.dna_sm.toplevel.fa

for file in ${Dir}*_aln_*; do
	name=${file##*/} # on garde tout le nom du fichier après le dernier /
	sample_tag=${name%_aln*} # on garde le tag de l'échantillon (avant le _aln)
	echo "_______________________ Sample "$sample_tag" _______________________"
    echo "_______________________ HTSBOX_________________________"
    #en in on prend le fichier .bam des mapped à q25
    #en out on obtient une séquence consensus du génome mitochondrial pour chaque échantillon

    source /local/env/envconda.sh
    conda activate /groups/Paleogenomics/ENV/htsbox
    htsbox pileup -f $GENOME -q20 -q30 -Fs3 $Dir${sample_tag}_sort_fixmate_MT.sorted.bam > $Dir${sample_tag}_htsbox.fasta

    source /local/env/envpython-3.9.5.sh
    python n2gaps.py $Dir${sample_tag}_htsbox.fasta $Dir${sample_tag}_htsbox_corrected.fasta # on applique la correction des N vers les gaps

done
