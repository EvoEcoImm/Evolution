#### This folder contain all the inforamtion about phylogenetic tree inference and dating.

Note: The file location and file name should be changed when necessary.

The matrix was prepared as described in the manuscript.
-The raw reads for each species at first were filtered with rRNA and mitochondrial sequences(sequences in *rRNA_mtDNA.zip* in this folder).
-Subsequently, the filtered raw reads were used for assembling with Trinity, redudancy reducing with CD-HIT-EST, and low expression transcript filtering with script in Trinity.
-Afterwards, the assemblies were filtered again with rRNA and mtDNA.
-Finally, the cleaned and filtered assemblies were used to predict proteins with Transdecoder(minium length 60 aa).
-The proteomes were used to infer orthology groups and the orthology groups were selected as described in the paper.
-The selected orthology groups were aligned by MAFFT with L-INS-i, masked by trimAl, and concatenated with Phytility. The matrix is available in Dryad and here, named as orthov4.phylip and orthov4.fasta.
-**To inference bayes tree with exabayes**
> `mpirun -np 48 ~/opt/exabayes-1.4.1/exabayes -f orthov4.phylip -m PROT -n TREEv4 -s 258 -c config.nex`
-**To inference maximum likelihood tree with RAxML**
>`~/opt/standard-RAxML/raxmlHPC-PTHREADS-AVX2 -f a -x 12345 -p 12345 -# 1000 -s orthov4.fasta -T 24 -m PROTGAMMAAUTO -n v4raxm1000b`
-**To date the phylogeny from former two methods** This had been run four time in our analysis, the calibration nodes were indicated in *Fossil selection.xlsx*.
>`~/opt/phylobayes4.1c/data/pb -d orthov4.phylip -T Treeforphylobayes -jtt -ugam -bd -sb -cal calibration jtt05071`

The output files are in the folder. The final time tree is *Time_Tree.newick*.
