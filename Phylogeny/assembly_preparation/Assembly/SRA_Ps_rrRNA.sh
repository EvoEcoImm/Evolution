#! /bin/bash
#SBATCH -J Ps_remove
#SBATCH -D /scratch/shulinhe/phylogeny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --time=6:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=shulinhe@zedat.fu-berlin.de
bowtie2 --very-sensitive-local -x ./rhinotermitidae -1 /scratch/shulinhe/SRAdata/trinity_ps_as/Prorhinotermes_simplex_1_edi.fastq.gz.P.qtrim.gz -2 /scratch/shulinhe/SRAdata/trinity_ps_as/Prorhinotermes_simplex_2_edi.fastq.gz.P.qtrim.gz --al-conc-gz ./SRA_Ps_rRNA --un-conc-gz ./SRA_Ps_trim --quiet
