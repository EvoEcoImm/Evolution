#! /home/she/opt/Python-3.4/bin//python3


import re
import os
sample2fullname={'Bg':'Blattella_germanica','Bo':'Blatta_orientalis','SRA_Ba':'Blaberus_atropos','SRA_Pa2016':'Periplaneta_americana South Korea','SRA_Es':'Eupolyphaga_sinensis','SRA_Cw':'Cryptocercus_wrighti','Cm':'Cryptocercus_meridianus','Ca':'Cryptocercus_pudacoensis','Md':'Mastotermes_darwiniensis','Cb':'Cryptotermes_brevis','SRA_Cd':'Cryptotermes_domesticus','Znev':'Zootermopsis nevadasis_nuttingi','Zns':'Zootermopsis_nevadasis','Hs':'Hodotermopsis_sjostedti','Nc':'Neotermes_castaneus','Kf':'Kalotermes_flavicollis','Rf':'Reticulitermes flavipes','Rs':'Reticulitermes_speratus','SRA_Rb':'Reticulitermes_banyulensis','SRA_Rf':'Reticulitermes_flavipes_SRA','SRA_Rg':'Reticulitermes_grassei','SRA_Rl':'Reticulitermes_lucifugus','Cf':'Coptotermes_formosanus','Pi':'Prorhinotermes_inopiinatus','SRA_Ps':'Prorhinotermes_simplex','Ms':'Macrotermes_subhyalinus','Nt':'Nasutitermes_takasagoensis','SRA_Of':'Odontotermes_formosanus','LS05':'Pericapri','LS17':'Indotermes','LS19':'Mirocapri','LS26':'Globitermes','LS29':'Nasutitermes','Mnat':'Macrotermes_natalensis','PG24':'Promirotermes','SRA_Hcb':'Empusa_pennata','SRA_Hcc':'Metallyticus_splendidus','SRA_Hc':'Hymenopus_coronatus'}	

"""filelinshi=input("gene_count_file_location_name")"""

with open("Orthogroups.GeneCount.csv",'rt') as file:
	ntotal=[entry.strip('\r\n').split('\t') for entry in file.readlines()]


num=[ntotal[0]]
for entry in ntotal[1:]:
	abc=[entry[0]]
	abc.extend([int(ab) for ab in entry[1:]])
	if max(abc[1:-1]) <2:
		abcd=abc
		num.append(abcd)

"""cover rate species"""
"""cover=float(input("Cover rate?"))
ncover=int(cover*(len(num[0])-2))
numID=[entry[0] for entry in num[1:] if entry[-1] >= ncover]
print ('No. of orthologus groups', len(numID))"""
	
"""num50ID=[entry[0] for entry in num if entry[-1] >= 20]

"""'''cover 75% species'''"""
num75ID=[entry[0] for entry in num if entry[-1] >= 30]"""

"""cover main group"""
mainspecies=['Md_he','Cb_he','SRA_Cd_he', 'Zns_he','Hs_454_he','Nc_he','Kf_he','Rf_he','Cf_he','Rs_454_he','Pi_he','SRA_Ps_he']
mainindexdict={e:i for i,e in enumerate(num[0]) if e in mainspecies}
maingroup=[entry for entry in num[1:] if entry[mainindexdict['Md_he']]==1 and (entry[mainindexdict['Cb_he']]==1 or entry[mainindexdict['SRA_Cd_he']]==1) and (entry[mainindexdict['Hs_454_he']]==1 or entry[mainindexdict['Zns_he']]==1) and (entry[mainindexdict['Nc_he']]==1 or entry[mainindexdict['Kf_he']]==1) and entry[mainindexdict['Cf_he']]==1 and (entry[mainindexdict['Rf_he']]==1 or entry[mainindexdict['Rs_454_he']]==1) and (entry[mainindexdict['Pi_he']]==1 or entry[mainindexdict['SRA_Ps_he']]==1)]

for i in range(len(maingroup[0])-2):
  a=[entry[i+1] for entry in maingroup]
  print (num[0][i+1], sum(a))

maingroupID=[entry[0] for entry in maingroup]
print ('No. of crucial orthologus groups', len(maingroupID))

"""orthofile=input("gene_orthologus_groups_file_locate_name")"""

with open("Orthogroups.csv",'rt') as file:
	maingroup2gene=[list(filter(None,entry)) for entry in (entry.strip('\r\n').split('\t') for entry in file.readlines()) if entry[0] in maingroupID]
	'''covergroup2gene=[list(filter(None,entry)) for entry in (entry.strip('\r\n').split('\t') for entry in file.readlines()) if entry[0] in numID]'''



"""allfasta=input("All fasta file")"""
from Bio import SeqIO
allfasta=list(SeqIO.parse("../all35.fasta",'fasta'))

os.mkdir("Orthogroup")
os.chdir("./Orthogroup/")

import re

for ID in maingroup2gene:
	orthogroupmain=[abcdre for abcdre in allfasta if abcdre.id in ID]
	for reads in orthogroupmain:
		rname=sample2fullname[re.sub('\_\d+','',reads.id)]
		reads.id=rname
		reads.name=reads.description=""
	SeqIO.write(orthogroupmain,ID[0],"fasta")



os.chdir("../")

"""os.mkdir("Cover_orthogroup")
os.chdir("./Cover_orthogroup/")
for aID in covergroup2gene:
	orthogroupcover=[re for re in allfasta if re.id in aID]
	for reads in orthogroupcover:
		rname=sample2fullname[re.sub('\_\d+','',reads.id)]
		reads.id=rname
		reads.name=reads.description=""
	SeqIO.write(orthogroupcover,aID[0],"fasta")





os.chdir("../")"""

"""'''---'Bg_he', 'Bo_he', 'SRA_Ba_he', '#SRA_Pa2015_he#', 'SRA_Pa2016_he','SRA_Es_he'
---'SRA_Cw_he','Cm_he','Ca_he'
---'Md_he','Cb_he','SRA_Cd_he','Zn_he', 'Zns_he','Hs_454_he','Nc_he','Kf_he'
---'Rf_he', 'Rs_454_he','SRA_Rb_he', 'SRA_Rf_he', 'SRA_Rg_he', 'SRA_Rl_he' ---'Cf_he','Pi_he', 'SRA_Ps_he' 
---'Ms_he',  'Nt_454_he',  'SRA_Mn_he', 'SRA_Of_he', 'LS05_he', 'LS17_he' ---'LS19_he', 'LS26_he', 'LS29_he', 'Mng_he', 'PG24_he' 
---'SRA_Hc_he', 'SRA_Hca_he', 'SRA_Hcb_he', 'SRA_Hcc_he','''"""
