#! /bin/bash
#SBATCH -J rrRNA1
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH -o filter.%j.out
#SBATCH --nodes=1-1
#SBATCH --ntasks=4
#SBATCH --mem=4000
#SBATCH --time=30:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
for i in `cat sample2`; do cd trinity_"$i"_as; bowtie2 --very-sensitive-local -x ../cryptocercus -f "$i"_Trinity_fil.fasta --al-gz ./"$i"_Trinity_fil_rRNA --un-gz ./"$i"_Trinity_fil_trim --quiet ; cd ../; done
