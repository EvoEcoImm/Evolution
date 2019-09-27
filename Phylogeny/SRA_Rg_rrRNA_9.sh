#! /bin/bash
#SBATCH -J Rg_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=2:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rg_as/Reticulitermes_grassei_9.fastq.gz.qtrim.fq --al-gz ./SRA_Rg_rRNA_9 --un-gz ./SRA_Rg_trim_9 --quiet

