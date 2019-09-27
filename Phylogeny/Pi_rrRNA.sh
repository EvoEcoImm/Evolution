#! /bin/bash
#SBATCH -J Pi_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=6:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./rhinotermitidae -1 /scratch/shulinhe/Sequence_2nd_run/data/trinity_pi_as/Pi_S13_R1_001.fastq.gz.P.qtrim.gz -2 /scratch/shulinhe/Sequence_2nd_run/data/trinity_pi_as/Pi_S13_R2_001.fastq.gz.P.qtrim.gz --al-conc-gz ./Pi_rRNA --un-conc-gz ./Pi_trim --quiet
