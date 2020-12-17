
with open("all_diamond_id_len",'rt') as file:
	id2len={entry[0]:float(entry[2]) for entry in (abc.strip('\n').split(' ') for abc in file.readlines())}

with open("all_immune_edi.csv",'rt') as file:
	allimm=[entry.strip('\n').split('\t') for entry in file.readlines()]

spec=list(set(entry[0] for entry in allimm))


for i in spec:
	specimm=[entry for entry in allimm if entry[0]==i]
	specimmgene=[entry[6] for entry in specimm]
	if i not in {"cm","nc"}:
		with open(i +".match",'rt') as file:
			specimmnr=[entry for entry in (entry.strip('\n').split('\t') for entry in file.readlines()) if entry[0] in specimmgene ]



	if i in {"cm","nc"}:
		with open("/media/shulinhe/DATA/QuanT/NR_Diamond/"+ i +".match",'rt') as file:
			specimmnr=[entry for entry in (entry.strip('\n').split('\t') for entry in file.readlines()) if entry[0] in specimmgene ]




	immgenenr=[]
	coverdict={}
	for iso in specimmgene:
		groupimmgene=[entry for entry in specimmnr if entry[0]==iso]
		if len(groupimmgene):
			bestimmgene=groupimmgene[0]
			cover=(float(bestimmgene[9])-float(bestimmgene[8]))/id2len[bestimmgene[1]]*100
			coverdict[iso]=str(cover)
			bestimmgene.extend([str(cover),str(id2len[bestimmgene[1]])])
			immgenenr.append(bestimmgene)




	for entry in specimm:
		if entry[6] in coverdict.keys(): 
			entry.append(coverdict[entry[6]])



	with open("all_Immunenrdiamondt.csv",'a+') as file:
		for entry in immgenenr:
			file.write("%s\t%s\n" % (i,'\t'.join(entry)))





	with open("all_Immunenr_cover.csv",'a+') as file:
		for entry in specimm:
			file.write("%s\n" % '\t'.join(entry))



