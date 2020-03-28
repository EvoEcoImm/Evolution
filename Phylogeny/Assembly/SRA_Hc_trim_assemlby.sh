#! /bin/bash
#SBATCH -J SRA_Hc_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o SRA_Hc_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=42000
#SBATCH --time=14:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --left SRA_Hc_trim_1.gz --right SRA_Hc_trim_2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_SRA_Hc_as &> SRA_Hc_log01.txt
