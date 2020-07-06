#! /usr/bin/python


import argparse
import re
from Bio import SeqIO
import subprocess
import os

def parser():
	args=argparse.ArgumentParser()
	args.add_argument('-i','--input',help='The name of the Hmmsearch output_file.')
	args.add_argument('-p','--blastp',help='The name of the blastp output_file')
	args.add_argument('-t','--trinotate',help='The name of the trinotate output_file')
	args.add_argument('-g2f','--gene2family',help="The gene to family file from immune database")
	args.add_argument('-f','--fasta',help='The name of the fasta file')
	args.add_argument('-o','--output',help='The name of the results output file')
	args=args.parse_args()
	return args


hmmfile = open(parser().input,'rt')
blastpfile= open(parser().blastp, 'rt')
trinotate=open(parser().trinotate,'rt')
dgene2family=open(parser().gene2family,'rt')
outfile = open(parser().output,'wt')


def hmmfil_dict(file):
	entries= [re.split("  *", entry) for entry in file.read().split('\n') if "#" not in entry and entry]
	entriesamp=[entry for entry in entries if "DEFs" in entry or "DROs" in entry or "Termicin" in entry or "Attacin" in entry]
	specamp=[(entry[18].split('~~')[0],entry[0],entry[2],float(entry[4]),float(entry[7]),entry[-1].split(':')[0],entry[20],entry[21],entry[22]) for entry in (entry for entry in entriesamp) if float(entry[4])<=0.00001 and float(entry[7]) <=0.001]
	entriesre=[entry for entry in entries if not ("DEFs" in entry or "DROs" in entry or "Termicin" in entry or "Attacin" in entry)]
	specres=[(entry[18].split('~~')[0],entry[0],entry[2],float(entry[4]),float(entry[7]),entry[-1].split(':')[0],entry[20],entry[21],entry[22]) for entry in (entry for entry in entriesre) if float(entry[4])<=0.00001 and float(entry[7]) <=0.001 and int(entry[-3].replace("len:",""))>=100]
	spec=specamp+specres
	"""spec=[(entry[18].split('~~')[0],entry[0],entry[2],float(entry[4]),float(entry[7]),entry[-1].split(':')[0],entry[20],entry[21],entry[22]) for entry in (entry for entry in entries) if float(entry[4])<=0.00001 and float(entry[7]) <=0.001]"""
	spec.sort(key=lambda spec:(spec[3],spec[4]),reverse=True)
	adict={entry[1]:(entry[0],)+entry[2:] for entry in spec}
	return adict

def g2f(file):
	g2f_dict={entry[0]:entry[1] for entry in (entry.strip('\n').split('\t') for entry in file.readlines())}
	return g2f_dict


def blastp_fil(file,g2f_dict):
	blastplist=[entry.strip('\n').split('\t') for entry in file.readlines()]
	blastp_re=[(entry[0],g2f_dict[entry[1]],float(entry[10])) for entry in blastplist]
	blastp_re.sort(key=lambda blastp_re:blastp_re[2],reverse=True)
	dict={entry[0]:(entry[1],entry[2]) for entry in blastp_re if entry[2]<=0.00001}
	return dict



"combine hmm and blastp"
def combhp(dict1,dict2):
	ph=[b+dict2[a]+(a,) for a,b in (entry for entry in dict1.items() if entry[0] in dict2.keys()) if dict2[a][0]==b[1]]
	ph.sort(key=lambda ph:(ph[2],ph[3]),reverse=True)
	fdict={entry[0]:entry[1:] for entry in ph}
	return fdict


"""def blastx_fil(file,g2f_dict):
	blastxlist=[entry.strip('\n').split('\t') for entry in file.readlines()]
	blastx_re=[(entry[0],g2f_dict[entry[1]],float(entry[-2])) for entry in blastxlist if float(entry[-2])<=0.001]
	blastx_re.sort(key=lambda blastx_re:blastx_re[2],reverse=True)
	xdict={entry[0]:(entry[1],entry[2]) for entry in blastx_re}
	return xdict


def combhpx(phlist, dict3):
	phx=[x+dict3[x[4]] for x in (entry for entry in phlist if entry[4] in dict3.keys()) if x[1]==dict3[x[4]][0]]
	phx.sort(key=lambda spec:(spec[2],spec[3]),reverse=True)
	fdict={entry[0]:entry[1:] for entry in phx}
	return fdict"""


	
def combiso(fdict):
	final=[(a,)+b for a,b in fdict.items()]
	isofinal=[entry[-1] for entry in final]
	return isofinal


def gettrinotate(file, isofinal1):
	ctrinotate=[[entry[0],entry[2],entry[6],entry[7],entry[4]] for entry in (entry.strip('\n').split('\t') for entry in file.readlines()[1:]) if entry[4] in isofinal1]
	trinotate_edi=[]
	for trans,blastx,blastp,pfam,transid in ctrinotate:
		d=[trans,transid]
		if blastx !=".":
			blastx1="|".join([blastx.split('^')[0],"_".join(blastx.split('^')[2:5]),blastx.split('^')[5].replace('RecName: Full=','').replace(';',''),blastx.split('^')[6].split(';')[1].replace(' ','')])
		if blastx ==".":
			blastx1="."
		if blastp !=".":
			blastp1=[blastp.split('^')[0],"_".join(blastp.split('^')[2:5]).replace(',','_'),blastp.split('^')[5].replace('RecName: Full=','').replace(';',''),blastp.split('^')[6].split(';')[1].replace(' ','')]
		if blastp ==".":
			blastp1=['.']*4
		if pfam !=".":
			pfam1=[entry.split('^') for entry in pfam.replace('E:','').split('`')]
			pfam1.sort(key=lambda pfam1:float(pfam1[-1]))
			pfams=[pfam1[0][0],pfam1[0][2],pfam1[0][4]]
		if pfam ==".":
			pfams=['.']*3
		d.extend(blastp1+pfams+[blastx1])
		trinotate_edi.append(d)
	return trinotate_edi

def combhpxt(tedi,fdict):
	tfinal=[list((entry[0],)+fdict[entry[0]]+tuple(entry[2:9])) for entry in tedi]
	outputlist=[]
	for entry in tfinal:
		if "Apaf" in entry[1]:
			if "apoptotic" in entry[13].lower():
				outputlist.append(entry)
		elif "ATG16" in entry[1]:
			if ("autophagy" in entry[13].lower() or "autophagy" in entry[16].lower() or "atg" in entry[13]):
				outputlist.append(entry)
		elif "Cactus" in entry[1]:
			if "cactus" in entry[13].lower() or "inhibitor" in entry[13].lower() or "cactus" in entry[13].lower() or "inhibitor" in entry[16].lower() or "relish" in entry[13].lower():
				outputlist.append(entry)
		elif "ATG14" in entry[1]:
			if "beclin" in entry[13].lower() or "radiation" in entry[13] or "autophagy" in entry[16].lower():
				outputlist.append(entry)
		elif "ATG2" in entry[1]:
			if "autophagy" in entry[13].lower():
				outputlist.append(entry)
		elif "ATG18" in entry[1]:
			if "repeat" in entry[13].lower() or "autophagy" in entry[13].lower():
				outputlist.append(entry)
		elif "ATG7" in entry[1]:
			if "ubiquitin" in entry[13].lower() or "atg" in entry[13].lower():
				outputlist.append(entry)
		elif "CLIPs" in entry[1]:
			if "allergen" not in entry[13].lower():
				outputlist.append(entry)
		elif "Caspar" in entry[1]:
			if "factor" in entry[13].lower():
				outputlist.append(entry)
		elif "CASPs" in entry[1]:
			if "caspase" in entry[13].lower():
				outputlist.append(entry)
		elif "Domeless" in entry[1]:
			if "receptor" in entry[13].lower():
				outputlist.append(entry)
		elif "Hopscoth" in entry[1]:
			if "hopscotch" in entry[13].lower() or "jak" in entry[13].lower():
				outputlist.append(entry)
		elif "HEP" in entry[1]:
			if "jnk" in entry[13].lower() or "hemip" in entry[13].lower():
				outputlist.append(entry)
		elif "JNK" in entry[1]:
			if "interacting" in entry[13].lower() or "stress" in entry[13].lower() or "jnk" in entry[13].lower() or "basket" in entry[13].lower():
				outputlist.append(entry)
		elif "Ird5" in entry[1]:
			if "inhibitor" in entry[13].lower() or "kappa" in entry[13].lower() or "relish" in entry[13].lower() or "factor" in entry[13].lower():
				outputlist.append(entry)
		elif "IAPs" in entry[1]:
			if "inhibitor" in entry[13].lower() or "iap" in entry[13].lower():
				outputlist.append(entry)
		elif "MLs" in entry[1]:
			if "allergen" not in entry[13]:
				outputlist.append(entry)
		elif "RELs" in entry[1]:
			if "kappa" in entry[13].lower() or "nuclear" in entry[13].lower():
				outputlist.append(entry)
		elif "PPOs" in entry[1]:
			if "allergen" not in entry[13].lower():
				outputlist.append(entry)
		elif "Pelle" in entry[1]:
			if "pelle" in entry[13]:
				outputlist.append(entry)
		elif "Stam" in entry[1]:
			if "adapter" in entry[13].lower():
				outputlist.append(entry)
		elif "SPZ" in entry[1]:
			if "spaetzle" in entry[13].lower() or "spz" in entry[13].lower():
				outputlist.append(entry)
		elif "SCRBs" in entry[1]:
			if "peste" not in entry[13].lower():
				outputlist.append(entry)
		elif "Stam" in entry[1]:
			if "signal" in entry[13].lower() or "stam" in entry[13].lower:
				outputlist.append(entry)
		elif "Tube" in entry[1]:
			if "synthetase" not in entry[13].lower() and "synthase" not in entry[13].lower():
				outputlist.append(entry)
		elif "TAK1" in entry[1]:
			if "mitogen" in entry[13].lower() or "kinase kinase kinase" in entry[13].lower():
				outputlist.append(entry)
		elif "Tube" in entry[1]:
			if "tube" in entry[13].lower() or "interleukin" in entry[13].lower():
				outputlist.append(entry)
		elif "TLR" in entry[1]:
			if "toll" in entry[13].lower():
				outputlist.append(entry)
		elif "ULK" in entry[1]:
			if "ulk" in entry[13].lower() or "atg" in entry[13].lower() or "autophagy" in entry[13].lower() or "unc" in entry[13].lower():
				outputlist.append(entry)
		else:
			outputlist.append(entry)
	output=[entry[:7]+entry[9:] for entry in outputlist]
	return output

def clust(outputlist,fastafile):
	iso=[entry[4] for entry in outputlist]
	fastaiso=[seq for seq in SeqIO.parse(fastafile,"fasta") if seq.id in iso]
	SeqIO.write(fastaiso,"test.fa","fasta")
	subprocess.call("perl ~/opt/iAssembler-v1.3.3.x64/iAssembler.pl -i test.fa",shell=True)
	with open("./test.fa_output/contig_member") as file:
		contig=[entry.strip('\n').split('\t') for entry in file.readlines()]
	isocls=[entry[1:] for entry in contig if len(entry)>2]
	if isocls==[]:
		for gene in outputlist:
			gene.insert(0,0)
	if isocls!=[]:
		dupiso=[entry for iso in isocls for entry in iso]
		for gene in outputlist:
			if gene[4] in dupiso:
				gene.insert(0,1)
			else:
				gene.insert(0,0)	
	os.system("rm -r mira_assembly test.fa_output")
	os.system("rm test.fa")
	return outputlist


def output_write_out(file, outputl):
	for entry in outputl:
		file.write('%s\n' % "\t".join(str(ite) for ite in entry))
	file.close()

hmmfil_dict=hmmfil_dict(hmmfile)
gene2family_dict=g2f(dgene2family)
blastp_dict=blastp_fil(blastpfile, gene2family_dict)
hmmblastp=combhp(hmmfil_dict,blastp_dict)
#blastx_dict=blastx_fil(blastxfile, gene2family_dict)
#hmmblastpx=combhpx(hmmblastp,blastx_dict)
hmmblastpiso=combiso(hmmblastp)
trinotatelist=gettrinotate(trinotate,hmmblastpiso)
finalhmmtrinotate=combhpxt(trinotatelist,hmmblastp)
dfinallist=clust(finalhmmtrinotate, parser().fasta)
output_write_out(outfile,dfinallist)










