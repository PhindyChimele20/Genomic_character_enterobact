WKD="/data"
cd $WKD

#quality check
fastqc -o "/nlustre/users/tshongwe/Prac/Microbial_Genomics_PGRG5" PGRG5_R1.fq.1 PGRG5_R2.fq.1                                          

#Run Trimmomatic on paired-end reads
java -jar $TRIMMOMATIC \
PE -phred33 -validatePairs PGRG5_R1.fq.1 PGRG5_R2.fq.1 PGRG5_R1.trim.fq.1  PGRG5_R1.unpaired.fq.1 PGRG5_R2.trim.fq.1 PGRG5_R2.unpaired.fq.1 LEADING:25 TRAILING:25 SLIDINGWINDOW:4:20 MINLEN:100
       
#mapping
bwa mem -t 16 Reference.fna PGRG5_R1.trim.fq.1 PGRG5_R2.trim.fq.1 > Mapped.sam

#converting sam file,indexing bam and sorting
samtools view -b Mapped.sam > Mapped.bam
samtools sort -@ 4 -o Mapped.sorted.bam Mapped.bam
samtools index Mapped.sorted.bam

#Generate consensus
samtools consensus -f fasta -o consensus.fa Mapped.sorted.bam

#Variant calling
REFERENCE="Reference.fna"
BAM="Mapped.sorted.bam"

bcftools mpileup -f $REFERENCE $BAM | bcftools call --ploidy 1 -mv -Oz -o variants.vcf.gz 

bcftools index variants.vcf.gz

bcftools filter -i 'AF>0.25' variants.vcf.gz -Oz -o filtered_variants.vcf.gz

#Variant annotation
bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\n' filtered_variants.vcf.gz > variants.bed
bedtools intersect -a variants.bed -b features.bed -wa -wb > annotated_variants.txt

#De-novo Genome assembly
spades.py -1 PGRG5_R1.fq -2 PGRG5_R2.fq --isolate -o spades_output -t 4
prokka spades_output/contigs.fasta --outdir prokka_output --prefix PGRG5

