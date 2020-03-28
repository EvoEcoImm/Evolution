#! /bin/bash
#SBATCH -J PG24_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20000
#SBATCH --time=14:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./termitidea -1 /scratch/shulinhe/Sequence_5th_run/trinity_PG24_as/PG24_S12_R1_001.fastq.gz.P.qtrim.gz -2 /scratch/shulinhe/Sequence_5th_run/trinity_PG24_as/PG24_S12_R2_001.fastq.gz.P.qtrim.gz --al-conc-gz ./PG24_rRNA --un-conc-gz ./PG24_trim --quiet
