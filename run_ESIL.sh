zcat *1.fastq.gz | bgzip -c > Pool1_all.1.gz
zcat *2.fastq.gz | bgzip -c > Pool1_all.2.gz
echo "Pooling complete"

bwa mem -t 4 ESIL_seq.fa Pool1_all.1.gz Pool1_all.2.gz -o ESIL.coller.sam
samtools view -bS -o ESIL.coller.bam ESIL.coller.sam
samtools sort -@ 4 -o ESIL.coller.sort.bam ESIL.coller.bam
samtools index ESIL.coller.sort.bam
samtools mpileup -vu -f ESIL_seq.fa ESIL.coller.sort.bam | bcftools call -mv -Oz > ESIL.coller.calls.vcf
bcftools index ESIL.coller.calls.vcf
bcftools query -f '%CHROM %POS %REF %ALT %QUAL [ %INDEL %DP %DP4]\n' ESIL.coller.calls.vcf -o ESIL.coller.calls.txt
echo "SNPs indexed"

for i in $(seq 1 5); do
	
	bwa mem -t 4 ESIL_seq.fa ESIL.${i}_1.fastq.gz ESIL.${i}_2.fastq.gz | samtools view -bS - > ESIL.${i}.bam 
	echo "BWA finished ${i}"
	samtools sort -@ 4 -o ESIL.${i}.sort.bam ESIL.${i}.bam
	echo "Sort finished ${i}"
	samtools index ESIL.${i}.sort.bam
	samtools mpileup -vu -f ESIL_seq.fa ESIL.${i}.sort.bam | bcftools call -m -T ESIL.coller.calls.vcf -Oz > ESIL.${i}.targcalls.vcf.gz
	echo "Mpileup finished ${i}"
	bcftools index ESIL.${i}.targcalls.vcf.gz
	bcftools query -f '%CHROM %POS %REF %ALT %QUAL [ %INDEL %DP %DP4]\n' ESIL.${i}.targcalls.vcf.gz -o ESIL.${i}targcalls.txt
	echo "All done ${1}"

done
