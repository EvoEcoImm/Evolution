#! /bin/bash
#SBATCH -J LS29_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=25000
#SBATCH --time=16:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./termitidea -1 /scratch/shulinhe/Sequence_5th_run/trinity_LS29_as/LS29_S11_R1_001.fastq.gz.P.qtrim.gz -2 /scratch/shulinhe/Sequence_5th_run/trinity_LS29_as/LS29_S11_R2_001.fastq.gz.P.qtrim.gz --al-conc-gz ./LS29_rRNA --un-conc-gz ./LS29_trim --quiet
