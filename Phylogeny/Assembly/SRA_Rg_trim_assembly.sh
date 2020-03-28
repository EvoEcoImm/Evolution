#! /bin/bash
#SBATCH -J SRA_Rg_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o SRA_Rg_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=41000
#SBATCH --time=16:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --single SRA_Rg_trim_1.gz,SRA_Rg_trim_2.gz,SRA_Rg_trim_3.gz,SRA_Rg_trim_4.gz,SRA_Rg_trim_5.gz,SRA_Rg_trim_6.gz,SRA_Rg_trim_7.gz,SRA_Rg_trim_8.gz,SRA_Rg_trim_9.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_SRA_Rg_as &>SRA_Rg_log01.txt

