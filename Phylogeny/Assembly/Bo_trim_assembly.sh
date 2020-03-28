#! /bin/bash
#SBATCH -J Bo_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o Bo_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=42000
#SBATCH --time=10:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --left Bo_trim_1.gz --right Bo_trim_2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_Bo_as &> Bo_log02.txt
