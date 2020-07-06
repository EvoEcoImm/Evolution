##---GO analysis---##
mkdir ./GOseq/
#extract go
head -n1 ../nc_full_trinotate_annotation_report.xls > ./GOseq/nc_full_metazoa_trinotation_report
cut -f1 ../nc_full_tax_nr_metazoa| grep -F -f - ../nc_full_trinotate_annotation_report.xls -w >> ./GOseq/nc_full_metazoa_trinotation_report

~/opt/Trinotate-v3.1.1/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls ./GOseq/nc_full_metazoa_trinotation_report -G --include_ancestral_terms > ./GOseq/nc_full_metazoa_go_annotations.txt


#gene length file
cut -f1 ./TW3_S14/quant.sf |grep -F -f - ../nc_full_trinity.Trinity.fasta.gene_trans_map > ./GOseq/nc_goseq.gene_trans_map 
~/opt/trinityrnaseq-Trinity-v2.5.1/util/misc/fasta_seq_length.pl ../nc_full_trinity.Trinity.fasta > ./GOseq/nc_full_Trinity.fasta.seq_lens
cut -f2 ./GOseq/nc_goseq.gene_trans_map|grep -F -f - ./GOseq/nc_full_Trinity.fasta.seq_lens > ./GOseq/nc_goseq.seq_lens

~/opt/trinityrnaseq-Trinity-v2.5.1/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map ./GOseq/nc_goseq.gene_trans_map --out_prefix ncn --name_sample_by_basedir CR1_S5/quant.sf CS1_S3/quant.sf CW1_S1/quant.sf TR1_S6/quant.sf TS1_S4/quant.sf TW1_S2/quant.sf CR2_S11/quant.sf CS2_S9/quant.sf CW2_S7/quant.sf TR2_S12/quant.sf TS2_S10/quant.sf TW2_S8/quant.sf CR3_S17/quant.sf CS3_S15/quant.sf CW3_S13/quant.sf TR3_S18/quant.sf TS3_S16/quant.sf TW3_S14/quant.sf
~/opt/trinityrnaseq-Trinity-v2.5.1/util/misc/TPM_weighted_gene_length.py --gene_trans_map ./GOseq/nc_goseq.gene_trans_map --trans_lengths ./GOseq/nc_goseq.seq_lens --TPM_matrix ncn.isoform.TMM.EXPR.matrix > ./GOseq/ncn_Trinity.gene_lengths.txt
head -n1 ./GOseq/ncn_Trinity.gene_lengths.txt > ./GOseq/ncn_metazoa_Trinity.gene_lengths.txt
cut -f1 ../nc_full_tax_nr_metazoa| grep -F -f - ./GOseq/ncn_Trinity.gene_lengths.txt -w >> ./GOseq/ncn_metazoa_Trinity.gene_lengths.txt

~/opt/trinityrnaseq-Trinity-v2.5.1/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map ./GOseq/nc_goseq.gene_trans_map --out_prefix nci --name_sample_by_basedir C1_S1/quant.sf nc_cp_1_S3/quant.sf nc_cr_2_S15/quant.sf nc_tp_1_S4/quant.sf nc_tr_2_S16/quant.sf T1_S4/quant.sf C2_S2/quant.sf nc_cp_2_S11/quant.sf nc_cs_1_S5/quant.sf nc_tp_2_S12/quant.sf nc_ts_1_S6/quant.sf T2_S5/quant.sf C3_S3/quant.sf nc_cr_1_S7/quant.sf nc_cs_2_S13/quant.sf nc_tr_1_S8/quant.sf nc_ts_2_S14/quant.sf T3_S6/quant.sf
~/opt/trinityrnaseq-Trinity-v2.5.1/util/misc/TPM_weighted_gene_length.py --gene_trans_map ./GOseq/nc_goseq.gene_trans_map --trans_lengths ./GOseq/nc_goseq.seq_lens --TPM_matrix nci.isoform.TMM.EXPR.matrix > ./GOseq/nci_Trinity.gene_lengths.txt
head -n1 ./GOseq/nci_Trinity.gene_lengths.txt > ./GOseq/nci_metazoa_Trinity.gene_lengths.txt
cut -f1 ../nc_full_tax_nr_metazoa| grep -F -f - ./GOseq/nci_Trinity.gene_lengths.txt -w >> ./GOseq/nci_metazoa_Trinity.gene_lengths.txt

#factor labelling
cut -f1 IS_TvsC_sig_up.csv -d ","|sed '1d'|awk '{print "Treatis" "\t" $1}' > ./GOseq/ncis_factor_labelling.txt
cut -f1 IS_TvsC_sig_down.csv -d ","| sed '1d'|awk '{print "Controlis" "\t" $1}' >> ./GOseq/ncis_factor_labelling.txt
cut -f1 wres_TvsC_sig_up.csv -d ","|sed '1d'|awk '{print "wTreatii" "\t" $1}' > ./GOseq/wii_factor_labelling.txt
cut -f1 wres_TvsC_sig_down.csv -d ","| sed '1d'|awk '{print "wControlii" "\t" $1}' >> ./GOseq/wii_factor_labelling.txt
cut -f1 sres_TvsC_sig_up.csv -d ","|sed '1d'|awk '{print "sTreatii" "\t" $1}' > ./GOseq/sii_factor_labelling.txt
cut -f1 sres_TvsC_sig_down.csv -d ","| sed '1d'|awk '{print "sControlii" "\t" $1}' >> ./GOseq/sii_factor_labelling.txt
cut -f1 rres_TvsC_sig_up.csv -d ","|sed '1d'|awk '{print "rTreatii" "\t" $1}' > ./GOseq/rii_factor_labelling.txt
cut -f1 rres_TvsC_sig_down.csv -d ","| sed '1d'|awk '{print "rControlii" "\t" $1}' >> ./GOseq/rii_factor_labelling.txt

cut -f1 RvsS_sig_up.csv -d ","|sed '1d'|awk '{print "RvsSup" "\t" $1}' > ./GOseq/RvsS_factor_labelling.txt
cut -f1 RvsS_sig_down.csv -d ","| sed '1d'|awk '{print "RvsSdown" "\t" $1}' >> ./GOseq/RvsS_factor_labelling.txt
cut -f1 RvsW_sig_up.csv -d ","|sed '1d'|awk '{print "RvsWup" "\t" $1}' > ./GOseq/RvsW_factor_labelling.txt
cut -f1 RvsW_sig_down.csv -d ","| sed '1d'|awk '{print "RvsWdown" "\t" $1}' >> ./GOseq/RvsW_factor_labelling.txt
cut -f1 WvsS_sig_up.csv -d ","|sed '1d'|awk '{print "WvsSup" "\t" $1}' > ./GOseq/WvsS_factor_labelling.txt
cut -f1 WvsS_sig_down.csv -d ","| sed '1d'|awk '{print "WvsSdown" "\t" $1}' >> ./GOseq/WvsS_factor_labelling.txt

cut -f1 TC_R_sig_up.csv -d ","|sed '1d'|awk '{print "RTreatn" "\t" $1}' > ./GOseq/rn_factor_labelling.txt
cut -f1 TC_R_sig_down.csv -d ","| sed '1d'|awk '{print "RControln" "\t" $1}' >> ./GOseq/rn_factor_labelling.txt
cut -f1 TC_W_sig_up.csv -d ","|sed '1d'|awk '{print "WTreatn" "\t" $1}' > ./GOseq/wn_factor_labelling.txt
cut -f1 TC_W_sig_down.csv -d ","| sed '1d'|awk '{print "WControln" "\t" $1}' >> ./GOseq/wn_factor_labelling.txt
cut -f1 TC_S_sig_up.csv -d ","|sed '1d'|awk '{print "STreatn" "\t" $1}' > ./GOseq/sn_factor_labelling.txt
cut -f1 TC_S_sig_down.csv -d ","| sed '1d'|awk '{print "SControln" "\t" $1}' >> ./GOseq/sn_factor_labelling.txt

#background file
cd ./GOseq/
cut -f1 ncn_metazoa_Trinity.gene_lengths.txt > ncn_metazoa_backgroup.txt
cut -f1 nci_metazoa_Trinity.gene_lengths.txt > nci_metazoa_backgroup.txt

#goseq in Trinity
##ncn gene_length file##
for i in {"ncis_factor_labelling.txt","wii_factor_labelling.txt","sii_factor_labelling.txt","rii_factor_labelling.txt"};
do ~/opt/trinityrnaseq-Trinity-v2.5.1/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling $i --GO_assignments nc_full_metazoa_go_annotations.txt --lengths nci_metazoa_Trinity.gene_lengths.txt --background nci_metazoa_backgroup.txt; done

##nci gene_length file##
for i in {"RvsS_factor_labelling.txt","RvsW_factor_labelling.txt","WvsS_factor_labelling.txt","rn_factor_labelling.txt","wn_factor_labelling.txt","sn_factor_labelling.txt"}; do
~/opt/trinityrnaseq-Trinity-v2.5.1/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling $i --GO_assignments nc_full_metazoa_go_annotations.txt --lengths ncn_metazoa_Trinity.gene_lengths.txt --background ncn_metazoa_backgroup.txt; done

##---plot_GO---##
ls *GOseq.enriched |cut -f1 -d "." > gosamplelist
for i in `cat gosamplelist`; do
cut -f1-9 "$i".GOseq.enriched |grep BP |awk 'BEGIN{FS="\t"}$8<=0.05' > "$i".goseq_BP.cut;
cut -f1-9 "$i".GOseq.enriched |grep MF |awk 'BEGIN{FS="\t"}$8<=0.05' > "$i".goseq_MF.cut;
done
#copy GOs to REVIGO in order to reduce GO redanduncy
#edit redanduncy GOs and export file "REVIGOs"

#'''grep -F -f REVIGOs Treat.GOseq.enriched > b'''

#cat treat.goseq*revigo.csv |cut -f1 -d "," | grep -F -f - Treat.GOseq.enriched > tgoplot
#cut -f1 treat.goseq_MF.cut|grep -F -f - Treat.GOseq.enriched >> tgoplot

#cat control.goseq_*_revigo.csv|cut -f1 -d ","|grep -F -f - Control.GOseq.enriched >> cgoplot

#R perform the plot
#Rscript ./GO_plot.R
