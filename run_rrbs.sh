###first do qc step
for i in *.fastq.gz
do 
    trim_galore -q 30 --phred33 --fastqc --illumina --clip_R1 4 --length 30 --rrbs -o /extData/NGS/KG_methy/ $i
done
###prepare the genome index
bismark_genome_preparation --bowtie2 --verbose /data/NGS/genomes/hg19refseq/  
###map the fq file to the genome with bismark
for i in *.fq.gz
do 
    bismark -l 32 -q --unmapped --ambiguous --phred33-quals --genome_folder /data/NGS/genomes/hg19refseq/ -p 5 --multicore 5 $i 
done
###if you need bedgraph use the following step
for i in *.sam
do 
    bismark_methylation_extractor -s --bedGraph --counts --buffer_size 30G --cytosine_report --genome_folder . --multicore 18 $i
done
##sort the bam file

for i in *.sam
do
        samtools sort -@ 40 -o ${i%%.bam}.sort.bam $i
done
### transfer the bam to the sam files as the methylkit need to read the file without header
for i in *.sort.bam
do
        samtools view -h $i >${i%%.bam}.sort.sam
done
####
####remove the header from sam file
for i in *.sort.sam
do
    grep -v '^@' $i > ${i%%.sort.sam}.trim.sam; 
done
###now input to the R by methylKit package

