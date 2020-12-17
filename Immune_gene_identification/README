This folder contain the procedures described in the manuscript about the immune gene identification and curation.

Note: the full names of abbrevations in samplename are indicated in Dryad README file. The file names and locations in the following commands were the personal working file names and locations, which should be modified for suitable use. 

- __clustalo database building__, *familynames* contains the gene families we predicted in the whole study, including a few families were splitted into more groups but we combine the results at last to report in the manuscript. *Immune_database_* contains the immune sequences downloaded to build up the searching database.
>     for i in `cat familynames`; do clustalo -i "$i".fa -o ../clustaloalign/"$i".clustalo1.alin --outfmt=st;sed "1a\#=GF ID $i" "$i".clustalo1.alin >> orthodbclusto.sto; done
>     hmmbuild orthodbclusto.hmm orthodbclusto.sto

- __protein blast database building__
>     cat *.fasta > orthodb.fa
>     makeblastdb -in orthodb.fa -dbtype prot

- __the raw reads were assembled with Trinity and annotated by following the Trinotate pipelines[https://github.com/Trinotate/Trinotate.github.io/wiki].__

- __the assemblies were also queried against nr database with diamond for later overlap detection in clusters__

- __Hmmsearch for protein set__
>     for i in `cat samplename`; do hmmsearch -o "$i".out --domtblout "$i".dom_out.tab --tblout "$i".target_out.tab --noali --notextw -E 1e-5 --domE 1 --incE 0.001 --cpu 4 /media/shulinhe/DATA/immune_database/orthodb_v6/orthodbclusto.hmm /media/shulinhe/DATA/QuanL/Transdecoder/"$i"_Trinity.fasta.transdecoder.pep; done

- __Blastp for protein set__
>     for i in `cat samplename`; do sed '/^#/d' "$i".target_out.tab |sed 's/  */ /g'|cut -f1 -d " "|sort -u > linshi.id; xargs samtools faidx /media/shulinhe/DATA/QuanL/Transdecoder/"$i"_Trinity.fasta.transdecoder.pep < linshi.id >> "$i".pep;
blastp -query "$i".pep -db /media/shulinhe/DATA/immune_database/orthodb_v6/Immune_v6.fa -num_threads 3 -evalue 0.00001 -max_target_seqs 1 -outfmt 6 -out "$i".db6.blastp.outfmt6; rm "$i".pep; done

- __prepared trinotate output, gene2family, fasta files__
>     for i in `cat samplename`; do
../Immune_curate.py -i "$i".target_out.tab -p "$i".db6.blastp.outfmt6 -t /media/shulinhe/DATA/QuanL/Trinotate_out/"$i"_trinotate_annotation_report.xls -g2f /media/shulinhe/DATA/immune_database/orthodb_v6/Immunedbv6gene2family -f /media/shulinhe/DATA/QuanL/Assembly/"$i"_Trinity.fasta -o "$i"_immune_curate; done

- __concatenate immune genes from each species together__
>     for i in `cat samplename`; do awk -v a="$i" '{print a "\t" $0}' "$i"_immune_curate >> ../all_db6_immune_precheck;done


- __replace uniprot name with uniport gene name, and manually check the uniprot gene names and predicted immune protein names__
>     cut -f12 all_db6_immune_precheck |sort -u |sed '/^\.$/d' > linshiuniprotid
>     python uniprotparse.py 
>     awk -F'\t' 'BEGIN{OFS="\t"}FNR==NR{a[$1]=$2"|"$3;next}{$14=a[$12];print}' uniprotfamilydescription all_db6_immune_precheck >all_db6_immune_precheckid; rm linshiuniprotid

- __name the checked file as *all_db6_immune_checkid_family.csv*, following the further curation steps with alignment of sequences in subcluster in trinity header and check of overlap of trinity__ subclusters in cluster.
> run *Refine_based_on_diamond_nr_blastp.R* in a Rstudio or other R interacitve window
>     python alignment_for_clusters.py

- __manually check the overlap of subclusters in a cluster and alignment of sequences in subcluster,name the file as *all_db6_immune_checkid_family_manrefine.csv*__

- __the last step is to check the gene expression with the predicted immune genes, the expession data are in the *Experssion* folder, and manually check identity of the genes with low experssion__

- __Finally, the summarize the total gene number__
>     python sum.py

All the predicted immune genes are in *Predicted_immune_genes.xlsx*, and the count numbers in *Data for Figure 2 and dots in Figure 1 and for CAFE Phylosignal analysis.
