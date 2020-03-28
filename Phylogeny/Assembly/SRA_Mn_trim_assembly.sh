#! /bin/bash
#SBATCH -J SRA_Mn_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o SRA_Mn_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=42000
#SBATCH --time=14:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --left SRA_Mn_trim_1.gz --right SRA_Mn_trim_2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_SRA_Mn_as &> SRA_Mn_log01.txt
