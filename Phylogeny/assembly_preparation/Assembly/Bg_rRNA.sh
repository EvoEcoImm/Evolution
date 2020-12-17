#! /bin/bash
#SBATCH -J remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH -o rmrrna_1.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=12:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./cockroach -1 /scratch/shulinhe/Sequence_2nd_run/data/trinity_bg_as/Bg_S6_R1_001.fastq.gz.P.qtrim.gz -2 /scratch/shulinhe/Sequence_2nd_run/data/trinity_bg_as/Bg_S6_R2_001.fastq.gz.P.qtrim.gz --un-conc-gz ./Bg_trim --quiet
