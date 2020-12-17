#! /bin/bash
#SBATCH -J SRA_Ba_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o SRA_Ba_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=42000
#SBATCH --time=16:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --left SRA_Ba_trim_1.gz --right SRA_Ba_trim_2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_SRA_Ba_as &> SRA_Ba_log01.txt
