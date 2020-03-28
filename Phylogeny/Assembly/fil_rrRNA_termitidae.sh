#! /bin/bash
#SBATCH -J rrRNA5
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH -o filter.%j.out
#SBATCH --nodes=1-1
#SBATCH --ntasks=4
#SBATCH --mem=4000
#SBATCH --time=2:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
for i in {"LS05","LS17","LS19","LS26","LS29","PG24"}; do cd trinity_"$i"_as; bowtie2 --very-sensitive-local -x ../termitidea -f "$i"_Trinity_fil.fasta --al-gz ./"$i"_Trinity_fil_rRNA --un-gz ./"$i"_Trinity_fil_trim --quiet ; cd ../; done
