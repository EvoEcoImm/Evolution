#! /bin/bash
#SBATCH -J Ca_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=6:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./cryptocercus -1 /scratch/shulinhe/Sequence_2nd_run/data/trinity_ca_as/Ca_S8_R1_001.fastq.gz.P.qtrim.gz -2 /scratch/shulinhe/Sequence_2nd_run/data/trinity_ca_as/Ca_S8_R1_001.fastq.gz.P.qtrim.gz --al-conc-gz ./Ca_rRNA --un-conc-gz ./Ca_trim --quiet
