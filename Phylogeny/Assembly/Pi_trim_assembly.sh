#! /bin/bash
#SBATCH -J Pi_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o Pi_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=42000
#SBATCH --time=06:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --left Pi_trim_1.gz --right Pi_trim_2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_Pi_as &> Pi_log03.txt
