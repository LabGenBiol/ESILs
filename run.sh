#!/bin/bash
snps="ChP_curated_snps.txt" # Path to tab separated file with the SNPs used for calling the COs
seq="ChP_seq.fa" # Path to indexed sequence, used as a reference in alignment. This should be indexed with bwa index.
prefix="JMJ14"  # Set prefix (name of your samples) here
samples=$(seq 1 192) # Set how many samples you have
start=10934455 # Set the first position of your interval
suffix1="targcalls.txt" # Do not change this
suffix2="mod.txt" # Do not change this
suffix3="allele_targ.complete.txt" # Do not change this

for i in ${samples}; do
  
  bwa mem -t 16 ChP_seq.fa ${prefix}_${i}_1.fastq.gz ${prefix}_${i}_2.fastq.gz | samtools view -bS - > ${prefix}.${i}.bam 
  echo "BWA finished ${i}"
  samtools sort -@ 16 -o ${prefix}.${i}.sort.bam ${prefix}.${i}.bam 
  echo "Sort finished ${i}"
  samtools index ${prefix}.${i}.sort.bam
  samtools mpileup -vu -f ChP_seq.fa ${prefix}.${i}.sort.bam | bcftools call -m -T ChP_coller.calls.vcf -Oz > ${prefix}.${i}.targcalls.vcf.gz
  echo "Mpileup finished ${i}"
  bcftools index ${prefix}.${i}.targcalls.vcf.gz
  bcftools query -f '%CHROM %POS %REF %ALT %QUAL [ %INDEL %DP %DP4]\n' ${prefix}.${i}.targcalls.vcf.gz -o ${prefix}.${i}.${suffix1}
  awk '{
      if ($1 != "chloroplast" && $1 != "mitochondria" && $6 != 1) {
          $1 = "3";
          $2 += 10934455;
          $3 = $3;
          n = split($8, dp4, ",");
          $4 = dp4[1] + dp4[2];
          new_col5 = $5;
          $6 = dp4[3] + dp4[4];
          print $1 "\t" $2 "\t" $3 "\t" $4 "\t" new_col5 "\t" $6;
      }
  }' "${prefix}.${i}.${suffix1}" > "${prefix}.${i}.${suffix2}"

  awk 'BEGIN { FS = OFS = "\t" }
     NR == FNR {                
         file[$1,$2] = $4 "\t" $6
         next
     }
     {                         
         key = $1 SUBSEP $2
         if (key in file) {
             split(file[key], values, "\t")
             $4 = values[1]
             $6 = values[2]
         }
         print;
     }' "${prefix}.${i}.${suffix2}" "$snps" > "${prefix}.${i}.${suffix3}"


    echo "All done ${i}"
done
ls *.${suffix3} > ${prefix}.samples.txt
R --slave --vanilla --args ${prefix}.samples.txt < AUTO.R
  
