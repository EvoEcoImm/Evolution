#! /bin/bash
#SBATCH -J Rg_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=12:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rg_as/Reticulitermes_grassei_1.fastq.gz.qtrim.fq --al-gz ./SRA_Rg_rRNA_1 --un-gz ./SRA_Rg_trim_1 --quiet
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rg_as/Reticulitermes_grassei_2.fastq.gz.qtrim.fq --al-gz ./SRA_Rg_rRNA_2 --un-gz ./SRA_Rg_trim_2 --quiet
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rg_as/Reticulitermes_grassei_3.fastq.gz.qtrim.fq --al-gz ./SRA_Rg_rRNA_3 --un-gz ./SRA_Rg_trim_3 --quiet
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rg_as/Reticulitermes_grassei_4.fastq.gz.qtrim.fq --al-gz ./SRA_Rg_rRNA_4 --un-gz ./SRA_Rg_trim_4 --quiet
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rg_as/Reticulitermes_grassei_5.fastq.gz.qtrim.fq --al-gz ./SRA_Rg_rRNA_5 --un-gz ./SRA_Rg_trim_5 --quiet
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rg_as/Reticulitermes_grassei_6.fastq.gz.qtrim.fq --al-gz ./SRA_Rg_rRNA_6 --un-gz ./SRA_Rg_trim_6 --quiet
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rg_as/Reticulitermes_grassei_7.fastq.gz.qtrim.fq --al-gz ./SRA_Rg_rRNA_7 --un-gz ./SRA_Rg_trim_7 --quiet
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rg_as/Reticulitermes_grassei_8.fastq.gz.qtrim.fq --al-gz ./SRA_Rg_rRNA_8 --un-gz ./SRA_Rg_trim_8 --quiet
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rg_as/Reticulitermes_grassei_9.fastq.gz.qtrim.fq --al-gz ./SRA_Rg_rRNA_9 --un-gz ./SRA_Rg_trim_9 --quiet

