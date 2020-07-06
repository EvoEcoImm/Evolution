This folder contains all the codes for differential expression analysis:
1.sample files:
	Neotermes_samples_file.txt
	Bo_samples_file.txt
	Cm_samplename.txt
2. all raw read files in each species were assembled seperately for further analysis, according to file "Description.txt".
3. all 3 assemblies were annotated by following Trinotate pipeline, according to file "Annotation.sh". The immune genes were identified using the same protocols of "immune gene identification" in the home directory. 
4. quantification of gene expression were conducted by using salmon, according to file "Quantification.sh".
5. Further differential expression anlaysis was done with DESeq2, according to file "DE.sh".
6. The further Gene Ontology analysis were conducted based on Differential expression output, according to files "Cm/BO/NC_metazoa_GO_analysis.sh" and GO_plot.R.
7. Other plotted figures were based on output of DEanalysis at step 5 and annotation from step 3.
