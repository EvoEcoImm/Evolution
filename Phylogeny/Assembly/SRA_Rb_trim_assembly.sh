#! /bin/bash
#SBATCH -J SRA_Rb_as
#SBATCH -D /scratch/shulinhe/phylogeny/
#SBATCH -o SRA_Rb_as.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=18000
#SBATCH --time=01:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
Trinity --seqType fq --max_memory 18G --single SRA_Rb_trim.gz --CPU 12 --normalize_reads --normalize_by_read_set --output trinity_SRA_Rb_as &>SRA_Rb_log02.txt

