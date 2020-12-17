#! /bin/bash
#SBATCH -J Cd_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=10:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./termite -1 /scratch/shulinhe/SRAdata/trinity_cd_as/Cryptotermes_domesticus_1_1_edi.fastq.gz.P.qtrim.gz -2 /scratch/shulinhe/SRAdata/trinity_cd_as/Cryptotermes_domesticus_1_2_edi.fastq.gz.P.qtrim.gz --al-conc-gz ./SRA_Cd_rRNA --un-conc-gz ./SRA_Cd_trim --quiet
