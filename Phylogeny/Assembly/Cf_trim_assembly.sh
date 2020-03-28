#! /bin/bash
#SBATCH -J Cf_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o Cf_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=42000
#SBATCH --time=14:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --left Cf_trim_1.gz --right Cf_trim_2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_Cf_as &> Cf_log01.txt
