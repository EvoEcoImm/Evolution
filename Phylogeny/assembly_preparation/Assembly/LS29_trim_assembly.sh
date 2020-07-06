#! /bin/bash
#SBATCH -J LS29t_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o LS29t.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=42000
#SBATCH --time=1-16:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --left LS29_trim.1.gz --right LS29_trim.2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_LS29_as &> LS29_log01.txt
