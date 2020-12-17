#! /usr/bin/python

import argparse
import re
from collections import Counter
from ete3 import NCBITaxa
ncbi=NCBITaxa()

def parser():
	args=argparse.ArgumentParser()
	args.add_argument('-i', '--input', help='The Diamond taxonomy output file.')
	args.add_argument('-r', '--contamination',help='The name of contamination file')
	args.add_argument('-o','--output',help='The name of metazoa file')
	args.add_argument('-n', '--number',help='The name of count in genus')
	args=args.parse_args()
	return args	

taxonomy_file=open(parser().input, 'rt')
contamination_file=open(parser().contamination,'wt')
output_file=open(parser().output,'wt')
count_file=open(parser().number,'wt')

def taxon_dict(file):
	isonr=[entry.strip('\n').split('\t') for entry in file.readlines()]
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

	gene_all=gene_Metazoa+gene_Other
	gene_tax=[entry[1].split()[0] for entry in gene_all]
	genuscount=[[a,b] for a,b in Counter(gene_tax).items()]
	genuscount.sort(key=lambda genuscount:genuscount[1],reverse=True)
	return gene_Other,genuscount,gene_Metazoa


gene_contamination, genus_count, gene_t_Metazoa=taxon_dict(taxonomy_file)



def gene_contamination_out(file,gene_Other):
	for a,b,i in gene_Other:
		file.write("%s\t%s\t%s\n" % (a,b,i))
	file.close()


def count_out(file, genus_count):
	for a,b in genus_count:
		file.write("%s\t%d\n" % (a,b))
	file.close()


def gene_metazoa_out(file,gene_Metazoa):
	for a,b,i in gene_Metazoa:
		file.write("%s\t%s\t%s\n" % (a,b,i))
	file.close()



gene_contamination_out(contamination_file,gene_contamination)

gene_metazoa_out(output_file,gene_t_Metazoa)

count_out(count_file, genus_count)
	
'''
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

'''

""" Other_not_bac_protist=[a for a in gene_Other if "2759" in a or "131567" in a]
bac_protist=[a for a in gene_Other if "2759" not in a and "131567" not in a] """
