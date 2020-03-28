#! /bin/bash
#SBATCH -J Rl_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=6:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rl_as/Reticulitermes_lucifugus_1_1_edi.fastq.gz.qtrim.fq --al-gz ./SRA_Rl_rRNA_1 --un-gz ./SRA_Rl_trim_1 --quiet
bowtie2 --very-sensitive-local -x ./rhinotermitidae -U /scratch/shulinhe/SRAdata/trinity_rl_as/Reticulitermes_lucifugus_2_1_edi.fastq.gz.qtrim.fq --al-gz ./SRA_Rl_rRNA_2 --un-gz ./SRA_Rl_trim_2 --quiet

