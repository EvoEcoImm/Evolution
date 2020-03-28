#! /bin/bash
#SBATCH -J Pa2015_1_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=6:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./cockroach -1 /scratch/shulinhe/SRAdata/trinity_pa20151_as/Periplaneta_americana_2015_1_1_edi.fastq.gz.P.qtrim.gz -2 /scratch/shulinhe/SRAdata/trinity_pa20151_as/Periplaneta_americana_2015_1_2_edi.fastq.gz.P.qtrim.gz --al-conc-gz ./SRA_Pa2015_1_rRNA --un-conc-gz ./SRA_Pa2015_1_trim --quiet
