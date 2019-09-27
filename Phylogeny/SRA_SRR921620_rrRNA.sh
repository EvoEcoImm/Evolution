#! /bin/bash
#SBATCH -J Hc2_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=6:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./mantids -1 /scratch/shulinhe/SRAdata/trinity_SRR921620/SRR921620_1.fastq.gz.P.qtrim.gz -2 /scratch/shulinhe/SRAdata/trinity_SRR921620/SRR921620_2.fastq.gz.P.qtrim.gz --al-conc-gz ./SRA_SRR921620_rRNA --un-conc-gz ./SRA_SRR921620_trim --quiet