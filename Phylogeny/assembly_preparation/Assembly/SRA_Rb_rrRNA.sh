#! /bin/bash
#SBATCH -J Rb_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=6:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rb_as/Reticulitermes_banyulensis_1_edi.fastq.gz.qtrim.fq --al-gz ./SRA_Rb_rRNA --un-gz ./SRA_Rb_trim --quiet
