#! /bin/bash
#SBATCH -J SRA_Rf_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=3:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rf_as_1/Reticulitermes_flavipes_2_1_edi.fastq.gz.qtrim.fq --al-gz ./SRA_Rf_rRNA_2 --un-gz ./SRA_Rf_trim_2 --quiet
