#! /bin/bash
#SBATCH -J Cb_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=30000
#SBATCH --time=11:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./termite -1 /scratch/shulinhe/Sequence_5th_run/trinity_cb_as/Cb_S14_R1_001.fastq.gz.P.qtrim.gz -2 /scratch/shulinhe/Sequence_5th_run/trinity_cb_as/Cb_S14_R2_001.fastq.gz.P.qtrim.gz --al-conc-gz ./Cb_rRNA --un-conc-gz ./Cb_trim --quiet
