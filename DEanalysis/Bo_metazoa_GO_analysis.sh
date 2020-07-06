##---GO analysis---##
mkdir ./GOseq/
#extract go
head -n1 ../bo_full_trinotate_annotation_report.xls > ./GOseq/bo_full_metazoa_trinotation_report
cut -f1 ../bo_full_tax_nr_metazoa| grep -F -f - ../bo_full_trinotate_annotation_report.xls -w >> ./GOseq/bo_full_metazoa_trinotation_report

~/opt/Trinotate-v3.1.1/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls ./GOseq/bo_full_metazoa_trinotation_report -G --include_ancestral_terms > ./GOseq/bo_full_metazoa_go_annotations.txt

#factor labelling
cut -f1 bon_res_TvsC_sig_up.csv -d ","|sed '1d'|awk '{print "Treatn" "\t" $1}' > ./GOseq/bon_factor_labelling.txt
cut -f1 bon_res_TvsC_sig_down.csv -d ","| sed '1d'|awk '{print "Controln" "\t" $1}' >> ./GOseq/bon_factor_labelling.txt
cut -f1 boi_res_TvsC_sig_up.csv -d ","|sed '1d'|awk '{print "Treatis" "\t" $1}' > ./GOseq/boi_factor_labelling.txt
cut -f1 boi_res_TvsC_sig_down.csv -d ","| sed '1d'|awk '{print "Controlis" "\t" $1}' >> ./GOseq/boi_factor_labelling.txt
cut -f1 bo_ii_TvsC_sig_up.csv -d ","|sed '1d'|awk '{print "Treatii" "\t" $1}' > ./GOseq/boii_factor_labelling.txt
cut -f1 bo_ii_TvsC_sig_down.csv -d ","| sed '1d'|awk '{print "Controlii" "\t" $1}' >> ./GOseq/boii_factor_labelling.txt

#gene length file
cut -f1 ./BoC1_S1/quant.sf |grep -F -f - ../bo_full_trinity.Trinity.fasta.gene_trans_map > ./GOseq/bo_goseq.gene_trans_map 
~/opt/trinityrnaseq-Trinity-v2.5.1/util/misc/fasta_seq_length.pl ../bo_full_trinity.Trinity.fasta > ./GOseq/bo_full_Trinity.fasta.seq_lens
cut -f2 ./GOseq/bo_goseq.gene_trans_map|grep -F -f - ./GOseq/bo_full_Trinity.fasta.seq_lens > ./GOseq/bo_goseq.seq_lens

~/opt/trinityrnaseq-Trinity-v2.5.1/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map ./GOseq/bo_goseq.gene_trans_map --out_prefix bo --name_sample_by_basedir Cn1_S7/quant.sf Cn2_S8/quant.sf Cn3_S9/quant.sf Tn1_S10/quant.sf Tn2_S11/quant.sf Tn3_S12/quant.sf BoC1_S1/quant.sf BoC2_S2/quant.sf BoT1_S3/quant.sf BoT2_S4/quant.sf Ci1_S1/quant.sf Ci2_S2/quant.sf Ci3_S3/quant.sf Ti1_S4/quant.sf Ti2_S5/quant.sf Ti3_S6/quant.sf

~/opt/trinityrnaseq-Trinity-v2.5.1/util/misc/TPM_weighted_gene_length.py --gene_trans_map ./GOseq/bo_goseq.gene_trans_map --trans_lengths ./GOseq/bo_goseq.seq_lens --TPM_matrix bo.isoform.TMM.EXPR.matrix > ./GOseq/bo_Trinity.gene_lengths.txt

head -n1 ./GOseq/bo_Trinity.gene_lengths.txt > ./GOseq/bo_metazoa_Trinity.gene_lengths.txt
cut -f1 ../bo_full_tax_nr_metazoa| grep -F -f - ./GOseq/bo_Trinity.gene_lengths.txt -w >> ./GOseq/bo_metazoa_Trinity.gene_lengths.txt

#background file
cd ./GOseq/
cut -f1 bo_metazoa_Trinity.gene_lengths.txt > bo_metazoa_backgroup.txt

#goseq in Trinity

~/opt/trinityrnaseq-Trinity-v2.5.1/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling bon_factor_labelling.txt --GO_assignments bo_full_metazoa_go_annotations.txt --lengths bo_metazoa_Trinity.gene_lengths.txt --background bo_metazoa_backgroup.txt

~/opt/trinityrnaseq-Trinity-v2.5.1/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling boi_factor_labelling.txt --GO_assignments bo_full_metazoa_go_annotations.txt --lengths bo_metazoa_Trinity.gene_lengths.txt --background bo_metazoa_backgroup.txt

~/opt/trinityrnaseq-Trinity-v2.5.1/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling boii_factor_labelling.txt --GO_assignments bo_full_metazoa_go_annotations.txt --lengths bo_metazoa_Trinity.gene_lengths.txt --background bo_metazoa_backgroup.txt

##---plot_GO---##
cut -f1-9 Treatn.GOseq.enriched |grep BP |awk 'BEGIN{FS="\t"}$8<=0.05' > treatn.goseq_BP.cut
cut -f1-9 Treatn.GOseq.enriched |grep MF |awk 'BEGIN{FS="\t"}$8<=0.05' > treatn.goseq_MF.cut
cut -f1-9 Controln.GOseq.enriched |grep BP |awk 'BEGIN{FS="\t"}$8<=0.05' > controln.goseq_BP.cut
cut -f1-9 Controln.GOseq.enriched |grep MF |awk 'BEGIN{FS="\t"}$8<=0.05' > controln.goseq_MF.cut

cut -f1-9 Treatis.GOseq.enriched |grep BP |awk 'BEGIN{FS="\t"}$8<=0.05' > treatis.goseq_BP.cut
cut -f1-9 Treatis.GOseq.enriched |grep MF |awk 'BEGIN{FS="\t"}$8<=0.05' > treatis.goseq_MF.cut
cut -f1-9 Controlis.GOseq.enriched |grep BP |awk 'BEGIN{FS="\t"}$8<=0.05' > controlis.goseq_BP.cut
cut -f1-9 Controlis.GOseq.enriched |grep MF |awk 'BEGIN{FS="\t"}$8<=0.05' > controlis.goseq_MF.cut

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
