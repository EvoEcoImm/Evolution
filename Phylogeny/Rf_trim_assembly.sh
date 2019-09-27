#! /bin/bash
#SBATCH -J Rf_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o Rf_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=42000
#SBATCH --time=06:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --left Rf_trim_1.gz --right Rf_trim_2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_Rf_as &> Rf_log02.txt
