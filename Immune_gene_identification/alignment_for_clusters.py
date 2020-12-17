import re
import collections
from Bio import SeqIO,AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import MafftCommandline
from StringIO import StringIO
import subprocess
import os

"""with open("all_species",'rt') as file:
	species=[entry.strip('\n') for entry in file.readlines() if entry!="cm1\n" and entry!="nc1\n"]"""


with open("all_db6_immune_checkid_family_manrefine.csv",'rt') as file:
	allimmune=[entry.strip('\n').split('\t') for entry in file.readlines()]
	species=list(set([entry[0] for entry in allimmune if entry[0]!="." and "quant" not in entry[0]]))


for speciesname in species:
	specimmune=[entry for entry in allimmune if entry[0]== speciesname and entry[1]!= "evidential"]
	speciso=[entry[6] for entry in specimmune]
	specfasta=[seq for seq in SeqIO.parse("/media/shulinhe/DATA/QuanL/Assembly/"+speciesname+"_Trinity.fasta","fasta") if seq.id in speciso]
	specpro=[entry[10] for entry in specimmune]
	specprofasta=[seq for seq in SeqIO.parse("/media/shulinhe/DATA/QuanL/Transdecoder/"+speciesname+"_Trinity.fasta.transdecoder.pep","fasta") if seq.id in specpro]
	SeqIO.write(specprofasta,"linshipro.pep","fasta")
	cluster=[re.sub(r'(TRINITY_DN\d+_c\d+)(\_g\d+)',r'\1', entry[2]) for entry in specimmune]
	subclureplicate=[item for item,count in collections.Counter(cluster).items() if count >1]
	for entry in subclureplicate:
		iso=[ent[6] for ent in specimmune if entry in ent[2]]
		pro=[ent[10] for ent in specimmune if entry in ent[2]]
		fastaiso=[seq for seq in specfasta if seq.id in iso]
		SeqIO.write(fastaiso,"linshi.fasta","fasta")
		subprocess.call("mafft --clustalout linshi.fasta >> "+speciesname+"Align",shell=True)
		fastapro=[seq for seq in specprofasta if seq.id in pro]
		SeqIO.write(fastapro,"linshi.fasta","fasta")
		subprocess.call("mafft --clustalout linshi.fasta >> "+speciesname+"Align",shell=True)


	subprocess.call("cd-hit -i linshipro.pep -o "+speciesname+"cd_hit -c 0.95 -n 5 -d 0 -T 4" ,shell=True)
