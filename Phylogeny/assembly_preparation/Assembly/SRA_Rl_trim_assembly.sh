#! /bin/bash
#SBATCH -J SRA_Rl_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o SRA_Rl_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=18000
#SBATCH --time=06:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 18G --single SRA_Rl_trim_1.gz,SRA_Rl_trim_2.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_SRA_Rl_as &>SRA_Rl_log02.txt

