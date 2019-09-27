#! /bin/bash
#SBATCH -J LS19t_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o LS19t.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=42000
#SBATCH --time=18:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --left LS19_trim.1.gz --right LS19_trim.2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_LS19_as &> LS19_log02.txt
