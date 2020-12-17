#! /bin/bash
#SBATCH -J SRA_Rf_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o SRA_Rf_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=18000
#SBATCH --time=14:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 18G --single SRA_Rf_trim_1.gz,SRA_Rf_trim_2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_SRA_Rf_as &>SRA_Rf_log01.txt

