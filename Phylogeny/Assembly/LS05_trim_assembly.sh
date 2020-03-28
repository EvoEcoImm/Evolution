#! /bin/bash
#SBATCH -J LS05t_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o LS05t_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=42000
#SBATCH --time=1-12:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --left LS05_trim.1.gz --right LS05_trim.2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_LS05_as &> LS05_log01.txt
