#! /bin/bash
#SBATCH -J LS05_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=6:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./termitidea -1 /scratch/shulinhe/Sequence_5th_run/trinity_LS05_as/LS05_S7_R1_001.fastq.gz.P.qtrim.gz -2 /scratch/shulinhe/Sequence_5th_run/trinity_LS05_as/LS05_S7_R2_001.fastq.gz.P.qtrim.gz --al-conc-gz ./LS05_rRNA --un-conc-gz ./LS05_trim --quiet
