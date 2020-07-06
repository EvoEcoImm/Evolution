##---GO analysis---##
#extract go
cut -f1 cmtax_NR_gene_contamination| grep -F -f - /media/shulinhe/DATA/QuanT/Annotation/cm_trinotate_annotation_report.xls -w -v > cm_host_trinotation_report

~/opt/Trinotate-v3.1.1/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls cm_host_trinotation_report -G --include_ancestral_terms > cm_go_annotations.txt

#factor labelling
cut -f1 cmres_TvsC_sig_up.csv -d ","|sed '1d'|awk '{print "Treat" "\t" $1}' > factor_labelling.txt

cut -f1 cmres_TvsC_sig_down.csv -d ","| sed '1d'|awk '{print "Control" "\t" $1}' >> factor_labelling.txt

#gene length file
~/opt/trinityrnaseq-Trinity-v2.5.1/util/misc/fasta_seq_length.pl /media/shulinhe/DATA/QuanT/Annotation/cm_Trinity.fasta > cm_Trinity.fasta.seq_lens

~/opt/trinityrnaseq-Trinity-v2.5.1/util/abundance_estimates_to_matrix.pl --est_method kallisto --gene_trans_map /media/shulinhe/DATA/QuanT/Annotation/cm_Trinity.fasta.gene_trans_map --out_prefix Cm --name_sample_by_basedir cm_c_1/abundance.tsv cm_c_2/abundance.tsv cm_t_1/abundance.tsv cm_t_2/abundance.tsv

~/opt/trinityrnaseq-Trinity-v2.5.1/util/misc/TPM_weighted_gene_length.py --gene_trans_map /media/shulinhe/DATA/QuanT/Annotation/cm_Trinity.fasta.gene_trans_map --trans_lengths cm_Trinity.fasta.seq_lens --TPM_matrix Cm.isoform.TMM.EXPR.matrix > cm_Trinity.gene_lengths.txt

cut -f1 cmtax_nr_gene_contamination| grep -F -f - cm_Trinity.gene_lengths.txt -w -v > cm_nbp_Trinity.gene_lengths.txt

#background file

cut -f1 cm_nbp_Trinity.gene_lengths.txt > nbp_backgroup.txt

#goseq in Trinity

~/opt/trinityrnaseq-Trinity-v2.5.1/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling factor_labelling.txt --GO_assignments cm_go_annotations.txt --lengths cm_nbp_Trinity.gene_lengths.txt --background nbp_backgroup.txt

##---plot_GO---##
cut -f1-9 Treat.GOseq.enriched |grep BP |awk 'BEGIN{FS="\t"}$8<=0.05' > treat.goseq_BP.cut

cut -f1-9 Treat.GOseq.enriched |grep MF |awk 'BEGIN{FS="\t"}$8<=0.05' > treat.goseq_MF.cut

cut -f1-9 Control.GOseq.enriched |grep BP |awk 'BEGIN{FS="\t"}$8<=0.05' > control.goseq_BP.cut

cut -f1-9 Control.GOseq.enriched |grep MF |awk 'BEGIN{FS="\t"}$8<=0.05' > control.goseq_MF.cut

#copy GOs to REVIGO in order to reduce GO redanduncy
#edit redanduncy GOs and export file "REVIGOs"

#'''grep -F -f REVIGOs Treat.GOseq.enriched > b'''

cat treat.goseq*revigo.csv |cut -f1 -d "," | grep -F -f - Treat.GOseq.enriched > tgoplot
cut -f1 treat.goseq_MF.cut|grep -F -f - Treat.GOseq.enriched >> tgoplot

cat control.goseq_*_revigo.csv|cut -f1 -d ","|grep -F -f - Control.GOseq.enriched >> cgoplot

#R perform the plot
Rscript ./GO_plot.R
