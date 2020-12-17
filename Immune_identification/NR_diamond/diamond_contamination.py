#! /usr/bin/python


import argparse
import re
from collections import Counter
from ete3 import NCBITaxa
ncbi=NCBITaxa()

def parser():
	args=argparse.ArgumentParser()
	args.add_argument('-i','--input',help='The name of the Hmmsearch output_file.')
	args.add_argument('-o','--output',help='The name of the results output file')
	args=args.parse_args()
	return args

taxfile = open(parser().input,'rt')
outfile = open(parser().output,'wt')

def tax(file):
	isonr=[entry.strip('\n').split('\t') for entry in file.readlines()]
	iso_Metazoa=[]
	iso_Other=[]
	for a,b,e in isonr:
		if b!="0":
			name=ncbi.get_lineage(b)
			c=ncbi.get_taxid_translator(name)
			d=[c[iname] for iname in name]
			if "Metazoa" in d:
				iso_Metazoa.append([a,d[-1],b])
			else:
				iso_Other.append([a,d[-1],b])
	
	return iso_Other


def count(gene_Metazoa, gene_Other):
	gene_all=gene_Metazoa+gene_Other
	gene_tax=[entry[1].split()[0] for entry in gene_all]
	genuscount=[[a,b] for a,b in Counter(gene_tax).items()]
	genuscount.sort(key=lambda genuscount:genuscount[1],reverse=True)
	return genuscount

def output_write_out(file, gene_Other):
	for a,b,i in gene_Other:
		file.write("%s\t%s\t%s\n" % (a,b,i))


contamination=tax(taxfile)
output_write_out(outfile, contamination)

"""iso_Metazoa=[]
iso_Other=[]
for a,b,e in isonr:
	if b!="0":
		name=ncbi.get_lineage(b)
		c=ncbi.get_taxid_translator(name)
		d=[c[iname] for iname in name]
		if "Metazoa" in d:
			iso_Metazoa.append([a,d[-1],b])
		else:
			iso_Other.append([a,d[-1],b])



isonredi=[[re.sub("_i\d*","", entry[0]),entry[1],float(entry[2])] for entry in isonr if entry[1]!="0"]
isonredi.sort(key=lambda isonredi:(isonredi[0],isonredi[2]),reverse=True)
genedict={entry[0]:(entry[1].split(';')[0],entry[2]) for entry in isonredi}
gene_Metazoa=[]
gene_Other=[]
for a,b in genedict.items():
	if b[0]!="0":
		name=ncbi.get_lineage(b[0])
		c=ncbi.get_taxid_translator(name)
		d=[c[iname] for iname in name]
		if "Metazoa" in d:
			gene_Metazoa.append([a,d[-1],b[0]])
		else:
			gene_Other.append([a,d[-1],b[0]])

Other_not_bac_protist=[a for a in gene_Other if "2759" in a or "131567" in a]
bac_protist=[a for a in gene_Other if "2759" not in a and "131567" not in a] 

with open("nctax_nr_gene_contamination",'wt') as file:
	for a,b,i in gene_Other:
		file.write("%s\t%s\t%s\n" % (a,b,i))




with open("nctax_nr_genus_statistic",'wt') as file:
	for a,b in genuscount:
		file.write("%s\t%d\n" % (a,b))"""


