#! /bin/bash
#SBATCH -J LS17_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=30000
#SBATCH --time=18:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./termitidea -1 /scratch/shulinhe/Sequence_5th_run/trinity_LS17_as/LS17_S8_R1_001.fastq.gz.P.qtrim.gz -2 /scratch/shulinhe/Sequence_5th_run/trinity_LS17_as/LS17_S8_R2_001.fastq.gz.P.qtrim.gz --al-conc-gz ./LS17_rRNA --un-conc-gz ./LS17_trim --quiet
