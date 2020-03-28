#! /bin/bash
#SBATCH -J Nc_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o Nc_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=88000
#SBATCH --time=10:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 89G --left Nc_trim_1.gz --right Nc_trim_2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_Nc_as &> Nc_log03.txt
