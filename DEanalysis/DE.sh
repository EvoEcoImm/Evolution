##Gene_trans_map to tx2gene##
awk '{ print $2 " " $1}' nc_full_trinity.Trinity.fasta.gene_trans_map | sed '1 i\TXNAME\tGENEID' | sed 's/ /\t/g' > tx2gene.tsv
awk '{ print $2 " " $1}' cm_full_trinity.Trinity.fasta.gene_trans_map | sed '1 i\TXNAME\tGENEID' | sed 's/ /\t/g' > tx2gene.tsv
awk '{ print $2 " " $1}' bo_full_trinity.Trinity.fasta.gene_trans_map | sed '1 i\TXNAME\tGENEID' | sed 's/ /\t/g' > tx2gene.tsv


##Contamination gene script##
Diamond_nr.py -i /media/shulinhe/DATA/Full_analysis/Nc_full/nc/nc_full.tax -o nc_contamination -n nc_full_genus_statistic
Diamond_nr.py -i /media/shulinhe/DATA/Full_analysis/Cm_full/cm/cm_full.tax -o cm_contamination -n cm_full_genus_statistic
Diamond_nr.py -i /media/shulinhe/DATA/Full_analysis/Bo_full/bo/bo_full.tax -o bo_contamination -n bo_full_genus_statistic


##Differential expression analysis##

Rscript ./nc_injected.R #Termites Injected DE analysis
Rscript ./nc_caste.R    #Termites caste and social infection-DE analysis
Rscript ./Cm_de.R       #Cryptocercus DE analysis
Rscript ./Bo_de.R       #Blatta orientalis DE analysis
