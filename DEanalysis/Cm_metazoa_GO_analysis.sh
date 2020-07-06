##---GO analysis---##
mkdir ./GOseq/
#extract go
head -n1 ../cm_full_trinotate_annotation_report.xls > ./GOseq/cm_full_metazoa_trinotation_report
cut -f1 ../cm_full_tax_nr_metazoa| grep -F -f - ../cm_full_trinotate_annotation_report.xls -w >> ./GOseq/cm_full_metazoa_trinotation_report

~/opt/Trinotate-v3.1.1/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls ./GOseq/cm_full_metazoa_trinotation_report -G --include_ancestral_terms > ./GOseq/cm_full_metazoa_go_annotations.txt

#factor labelling
cut -f1 cmres_TvsC_sig_up.csv -d ","|sed '1d'|awk '{print "Treatii" "\t" $1}' > ./GOseq/cmii_factor_labelling.txt
cut -f1 cmres_TvsC_sig_down.csv -d ","| sed '1d'|awk '{print "Controlii" "\t" $1}' >> ./GOseq/cmii_factor_labelling.txt

#gene length file
~/opt/trinityrnaseq-Trinity-v2.5.1/util/misc/fasta_seq_length.pl ../cm_full_trinity.Trinity.fasta > ./GOseq/cm_full_Trinity.fasta.seq_lens

~/opt/trinityrnaseq-Trinity-v2.5.1/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map ../cm_full_trinity.Trinity.fasta.gene_trans_map --out_prefix cm --name_sample_by_basedir cm_c_1_S1/quant.sf cm_c_2_S9/quant.sf cm_t_1_S2/quant.sf cm_t_2_S10/quant.sf

~/opt/trinityrnaseq-Trinity-v2.5.1/util/misc/TPM_weighted_gene_length.py --gene_trans_map ../cm_full_trinity.Trinity.fasta.gene_trans_map --trans_lengths ./GOseq/cm_full_Trinity.fasta.seq_lens --TPM_matrix cm.isoform.TMM.EXPR.matrix > ./GOseq/cm_Trinity.gene_lengths.txt

head -n1 ./GOseq/cm_Trinity.gene_lengths.txt > ./GOseq/cm_metazoa_Trinity.gene_lengths.txt
cut -f1 ../cm_full_tax_nr_metazoa| grep -F -f - ./GOseq/cm_Trinity.gene_lengths.txt -w >> ./GOseq/cm_metazoa_Trinity.gene_lengths.txt

#background file
cd ./GOseq/
cut -f1 cm_metazoa_Trinity.gene_lengths.txt > cm_metazoa_backgroup.txt

#goseq in Trinity

~/opt/trinityrnaseq-Trinity-v2.5.1/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling cmii_factor_labelling.txt --GO_assignments cm_full_metazoa_go_annotations.txt --lengths cm_metazoa_Trinity.gene_lengths.txt --background cm_metazoa_backgroup.txt

##---plot_GO---##

cut -f1-9 Treatii.GOseq.enriched |grep BP |awk 'BEGIN{FS="\t"}$8<=0.05' > treatii.goseq_BP.cut
cut -f1-9 Treatii.GOseq.enriched |grep MF |awk 'BEGIN{FS="\t"}$8<=0.05' > treatii.goseq_MF.cut
cut -f1-9 Controlii.GOseq.enriched |grep BP |awk 'BEGIN{FS="\t"}$8<=0.05' > controlii.goseq_BP.cut
cut -f1-9 Controlii.GOseq.enriched |grep MF |awk 'BEGIN{FS="\t"}$8<=0.05' > controlii.goseq_MF.cut
#copy GOs to REVIGO in order to reduce GO redanduncy
#edit redanduncy GOs and export file "REVIGOs"

#'''grep -F -f REVIGOs Treat.GOseq.enriched > b'''

#cat treat.goseq*revigo.csv |cut -f1 -d "," | grep -F -f - Treat.GOseq.enriched > tgoplot
#cut -f1 treat.goseq_MF.cut|grep -F -f - Treat.GOseq.enriched >> tgoplot

#cat control.goseq_*_revigo.csv|cut -f1 -d ","|grep -F -f - Control.GOseq.enriched >> cgoplot

#R perform the plot
#Rscript ./GO_plot.R
