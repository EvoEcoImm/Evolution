
import collections

with open("all_db6_immune_checkid_family_manrefine.csv") as file:
	alllist0=[entry.strip('\n').split('\t') for entry in file.readlines()]

alllist=[]
for entry in alllist0:
	if entry[0]==".":
		entry[0]=entry[2][:4]
		alllist.append(entry)
	elif entry[1]=="evidential":
		entry[0]=entry[0]+entry[1]
	else:
		alllist.append(entry)



speciesname=list(set(entry[0] for entry in alllist))

family=set([entry[3] for entry in alllist])

countdict={entry:["0"] for entry in family}
countdict["family"]=["species"]
countdict["total"]=["0"]
for name in speciesname:
	countdict["family"].append(name)
	spec=[entry for entry in alllist if entry[0]==name]
	countdict["total"].append(str(len(spec)))
	for fam in family:
		countdict[fam].append(str(len([entry for entry in spec if entry[3]==fam])))



for entry in countdict.items():
	entry[1].insert(0,entry[0])



allcount=[]
for i in ("family","total","Attacin","Termicin","DEFs","DROs","Transferrin","TEPs","LYSs","destabilase","FREPs","PPOs","HPXs","TPXs","GPXs","SODs","CATs","MLs","CTLs","GALEs","PGRPs","GNBP","SPZs_Toll","TLR_Toll","SCRAs","SCRBs","SCRCs","Myd88_Toll","Cactus_Toll","Traf_Toll","JNK_ip_Toll","Dif_Toll","Pelle_Toll","RELs_Toll","Tube_Toll","Pellino-Toll","Key_IMD","Caspar_IMD","Imd_IMD","TAB2_IMD","Ird5_IMD","Fadd_IMD","TAK1_IMD","HEP","JNK","Stam_JAK-STAT","Hopscoth_JAK-STAT","Domeless_JAK-STAT","STAT_JAK-STAT","IAPs","Apaf-caspas","CASPs","ATG13","ATG18B","ATG2","ATG3","ATG6","ATG7","ATG5","ATG8","ATG9","ATG12","ATG14","ATG16","ULK_ATG","PepC54_ATG","ATG4b","CLIPs","SRPNs"):
	allcount.append(countdict[i])


with open("allmanrefine_count.csv",'wt') as file:
	for entry in allcount:
		file.write("%s\n" % '\t'.join(entry))




