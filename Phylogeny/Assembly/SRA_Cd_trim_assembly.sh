#! /bin/bash
#SBATCH -J SRA_Cd_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o SRA_Cd_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=42000
#SBATCH --time=10:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --left SRA_Cd_trim_1.gz --right SRA_Cd_trim_2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_SRA_Cd_as &> SRA_Cd_log05.txt
