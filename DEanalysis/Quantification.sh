Salmon

wget https://github.com/COMBINE-lab/salmon/releases/download/v0.12.0/salmon-0.12.0_linux_x86_64.tar.gz
tar zxvf salmon-0.12.0_linux_x86_64.tar.gz

salmon index -t Bo_full_trinity.Trinity.fasta -i Bo_full_index
while IFS= read -r i; do
salmon quant -i Bo_full_index -l a -1 ${i}_R1_001.fastq.gz -2 ${i}_R2_001.fastq.gz -p 12 -o ./Bo/$i
done < Bo_kallisto.txt

salmon index -t Cm_full_trinity.Trinity.fasta -i Cm_full_index
while IFS= read -r i; do
salmon quant -i Cm_full_index -l a -1 ${i}_R1_001.fastq.gz -2 ${i}_R2_001.fastq.gz -p 12 -o ./Cm/$i
done < Cm_kallisto.txt

salmon index -t Nc_full_trinity.Trinity.fasta -i Nc_full_index
while IFS= read -r i; do
salmon quant -i Nc_full_index -l a -1 ${i}_R1_001.fastq.gz -2 ${i}_R2_001.fastq.gz -p 12 -o ./Nc/$i
done < Nc_kallisto.txt


'''Bo  transcripts qualification'''
'''prepare reference'''
~/opt/trinityrnaseq-Trinity-v2.5.1/util/align_and_estimate_abundance.pl --transcripts bo_full_trinity.Trinity.fasta --est_method kallisto --trinity_mode --prep_reference --output_dir ./
'''align library'''
while IFS= read -r i; do
~/opt/trinityrnaseq-Trinity-v2.5.1/util/align_and_estimate_abundance.pl --transcripts bo_full_trinity.Trinity.fasta --seqType fq --left ${i}_R1_001.fastq.gz --right ${i}_R2_001.fastq.gz --est_method kallisto --trinity_mode --output_dir ./${i}
done < Bo_kallisto.txt

'''Nc  transcripts qualification'''
'''prepare reference'''
~/opt/trinityrnaseq-Trinity-v2.5.1/util/align_and_estimate_abundance.pl --transcripts Nc_full_trinity.Trinity.fasta --est_method kallisto --trinity_mode --prep_reference --output_dir ./
'''align library'''
while IFS= read -r i; do
~/opt/trinityrnaseq-Trinity-v2.5.1/util/align_and_estimate_abundance.pl --transcripts Nc_full_trinity.Trinity.fasta --seqType fq --left ${i}_R1_001.fastq.gz --right ${i}_R2_001.fastq.gz --est_method kallisto --trinity_mode --output_dir ./${i}
done < Nc_kallisto.txt

'''Cm  transcripts qualification'''
'''prepare reference'''
~/opt/trinityrnaseq-Trinity-v2.5.1/util/align_and_estimate_abundance.pl --transcripts Cm_full_trinity.Trinity.fasta --est_method kallisto --trinity_mode --prep_reference --output_dir ./
'''align library'''
while IFS= read -r i; do
~/opt/trinityrnaseq-Trinity-v2.5.1/util/align_and_estimate_abundance.pl --transcripts Cm_full_trinity.Trinity.fasta --seqType fq --left ${i}_R1_001.fastq.gz --right ${i}_R2_001.fastq.gz --est_method kallisto --trinity_mode --output_dir ./${i}
done < Cm_kallisto.txt
