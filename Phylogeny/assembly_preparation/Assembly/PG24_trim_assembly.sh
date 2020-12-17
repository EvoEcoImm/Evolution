#! /bin/bash
#SBATCH -J PG24t_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o PG24t.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=42000
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --left PG24_trim.1.gz --right PG24_trim.2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_PG24_as &> PG24_log01.txt
