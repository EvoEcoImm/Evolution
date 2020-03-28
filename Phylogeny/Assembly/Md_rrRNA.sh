#! /bin/bash
#SBATCH -J Md_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=6:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./termite -1 /scratch/shulinhe/Sequence_2nd_run/data/trinity_md_as/Md_S9_R1_001.fastq.gz.P.qtrim.gz -2 /scratch/shulinhe/Sequence_2nd_run/data/trinity_md_as/Md_S9_R1_001.fastq.gz.P.qtrim.gz --al-conc-gz ./Md_rRNA --un-conc-gz ./Md_trim --quiet
