"""Dimonad blast against nr database"""
git clone https://github.com/bbuchfink/diamond.git
wget ftp://ftp.ncbi.nlm.nih.gov/blast//db/FASTA/nr.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
diamond makedb --in nr.gz -d nr --taxonmap prot.accession2taxid.gz --taxonnodes nodes.dmp
'''this takes 42G and 12CPU, protein matches, it takes around 1day 12hour'''
diamond blastx -d nr -q Bo_full_trinity.Trinity.fasta -o Bo_full.match -f 6 qseqid sseqid pident length mismatch gapopen qlen qstart qend qcovhsp slen sstart send stitle staxids evalue bitscore --more-sensitive -c 1 &> Bo_full_log.txt
diamond blastx -d nr -q Cm_full_trinity.Trinity.fasta -o Cm_full.match -f 6 qseqid sseqid pident length mismatch gapopen qlen qstart qend qcovhsp slen sstart send stitle staxids evalue bitscore --more-sensitive -c 1 &> Cm_full_log.txt
diamond blastx -d nr -q Nc_full_trinity.Trinity.fasta -o Nc_full.match -f 6 qseqid sseqid pident length mismatch gapopen qlen qstart qend qcovhsp slen sstart send stitle staxids evalue bitscore --more-sensitive -c 1 &> Nc_full_log.txt
'''this is for taxonomy composition analysis'''
diamond blastx -d nr -q Bo_full_trinity.Trinity.fasta -o Bo_full.tax -f 102 --more-sensitive -c 1 &> Bo_tax_log.txt
diamond blastx -d nr -q Cm_full_trinity.Trinity.fasta -o Cm_full.tax -f 102 --more-sensitive -c 1 &> Cm_tax_log.txt
diamond blastx -d nr -q Nc_full_trinity.Trinity.fasta -o Nc_full.tax -f 102 --more-sensitive -c 1 &> Nc_tax_log.txt


'''Trinotate annotation'''
wget https://github.com/Trinotate/Trinotate/archive/Trinotate-v3.1.1.tar.gz
tar zxvf Trinotate-v3.1.1.tar.gz
wget https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v5.5.0.tar.gz
tar zxvf TransDecoder-v5.5.0.tar.gz
'''preparation database'''
$TRINOTATE_HOME/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate
makeblastdb -in uniprot_sprot.pep -dbtype prot
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
'''Blastx'''
blastx -query Bo_full_trinity.Trinity.fasta -db /path/to/uniprot_sprot.pep -num_threads 12 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > Bo_full_blastx.outfmt6
blastx -query Cm_full_trinity.Trinity.fasta -db /path/to/uniprot_sprot.pep -num_threads 12 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > Cm_full_blastx.outfmt6
blastx -query Nc_full_trinity.Trinity.fasta -db /path/to/uniprot_sprot.pep -num_threads 12 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > Nc_full_blastx.outfmt6
'''ORF prediction'''
TransDecoder.LongOrfs -t Bo_full_trinity.Trinity.fasta -m 60 -O Bo_full_predict
TransDecoder.LongOrfs -t Cm_full_trinity.Trinity.fasta -m 60 -O Cm_full_predict
TransDecoder.LongOrfs -t Nc_full_trinity.Trinity.fasta -m 60 -O Nc_full_predict
'''Blastp'''
blastp -query Bo_full_predict/longest_orfs.pep -db /path/to/uniprot_sprot.pep -num_threads 12 -max_target_seqs 1 -evalue 1e-5 -outfmt 6 > Bo_full_blastp.outfmt6
blastp -query Cm_full_predict/longest_orfs.pep -db /path/to/uniprot_sprot.pep -num_threads 12 -max_target_seqs 1 -evalue 1e-5 -outfmt 6 > Cm_full_blastp.outfmt6
blastp -query Nc_full_predict/longest_orfs.pep -db /path/to/uniprot_sprot.pep -num_threads 12 -max_target_seqs 1 -evalue 1e-5 -outfmt 6 > Nc_full_blastp.outfmt6
'''Protein domain search'''
hmmscan --cpu 12 --domtblout Bo_full_TrinotatePFAM.out /path/to/Pfam-A.hmm Bo_full_predict/longest_orfs.pep > Bo_full_pfam.log
hmmscan --cpu 12 --domtblout Cm_full_TrinotatePFAM.out /path/to/Pfam-A.hmm Cm_full_predict/longest_orfs.pep > Cm_full_pfam.log
hmmscan --cpu 12 --domtblout Nc_full_TrinotatePFAM.out /path/to/Pfam-A.hmm Nc_full_predict/longest_orfs.pep > Nc_full_pfam.log
'''Protein prediction'''
Transdecoder.Predict -t Bo_full_trinity.Trinity.fasta --retain_pfam_hits Bo_full_TrinotatePFAM.out --retain_blastp_hits Bo_full_blastp.outfmt6 -O Bo_full_predict
Transdecoder.Predict -t Cm_full_trinity.Trinity.fasta --retain_pfam_hits Cm_full_TrinotatePFAM.out --retain_blastp_hits Cm_full_blastp.outfmt6 -O Cm_full_predict
Transdecoder.Predict -t Nc_full_trinity.Trinity.fasta --retain_pfam_hits Nc_full_TrinotatePFAM.out --retain_blastp_hits Nc_full_blastp.outfmt6 -O Nc_full_predict
