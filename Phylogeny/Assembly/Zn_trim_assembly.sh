#! /bin/bash
#SBATCH -J Znt_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o Znt.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=42000
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 42G --left Zn_trim.1.gz --right Zn_trim.2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_Zn_as &> Zn_log02.txt
