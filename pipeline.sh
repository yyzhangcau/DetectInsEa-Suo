# 1. prepare genome and insert fasta file
#mkdir genome
#wget -c -P $PWD/genome/ -O EA.genome.fa https://toxodb.org/common/downloads/Current_Release/EacervulinaHoughton/fasta/data/ToxoDB-59_EacervulinaHoughton_Genome.fasta
REF=$PWD/genome/ToxoDB-59_EacervulinaHoughton_Genome.fasta
#bwa index $REF
#samtools faidx $REF
#makeblastdb -dbtype nucl -parse_seqids -in $REF
insert=$PWD/genome/genome_insert.seq
bwa index $insert


# 2 align ngs data against the insert sequence
R1=$PWD/data/R1.fq
R2=$PWD/data/R2.fq
bwa mem -t 4 -k 31 $insert $R1 $R2 |samtools view -h -q 5 - |samtools sort - |awk '{print $1}' | sort |uniq > reads.id.txt
seqtk subseq $R1 reads.id.txt > extract.reads.R1.fq
seqtk subseq $R2 reads.id.txt > extract.reads.R2.fq

# 3 mask homologous sequences in ref genome
blastn -query $insert -db $REF -perc_identity 0.98 -outfmt 6 -windown_size 31 |awk 'OFS="\t"{if($9<$10)print $2,$9,$10;else print $2,$10,$9}' > mask.bed
maskFastaFromBed -fi $REF -bed mask.bed -fo ${REF}.masked
# 4 prepare concatenate sequences
awk 1 ${REF}.masked $insert > concatenateSeq.fa 


makeblastdb -dbtype nucl -parse_seqids -in concatenateSeq.fa

# re-assembly insert and flanking sequences
spades.py -1 extract.reads.R1.fq -2 extract.read.R2.fq -o reassembled.out
blastn -query reassembled.out/scaffolds.fasta -db concatenateSeq.fa -outfmt 6 -out results.out

