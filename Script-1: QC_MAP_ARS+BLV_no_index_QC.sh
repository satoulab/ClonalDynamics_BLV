# ************ Main Script Begins Here ************ #
#
	echo
	echo "Begin analysis"
	echo "Task started at"
		date
	begin=$(date +%s)

# Prepare fastq files
	echo "Unzip"
		gunzip *gz
	echo "convert fastq file name"
		mv  *R1* R1.fastq
		mv  *R2* R2.fastq

# Make analysis directories
	echo
	echo "Create analysis directories"
		mkdir 1_qc
		mkdir 2_map

	echo "Directories 1_qc & 2_map created"
	echo "Begin step 1 in 1_qc directory"
		cd 1_qc

# Remove adaptor sequence
	echo
	echo "Remove adaptor from read 1"
	echo
		cutadapt -m 1 -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 5 ../R1.fastq -o R1_step1.fastq
	echo
	echo "Remove adaptor from read 2"
	echo
		cutadapt -m 1 -b GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 ../R2.fastq -o R2_step1.fastq

# Checking quality of the fastq files
	echo "Doing FastQC check"
	echo
	fastqc R1_step1.fastq
	fastqc R2_step1.fastq
	echo	

# Read QC
	echo "Reads quality check"
	echo
		/mnt/nas/programs/prinseq-lite-0.20.4/prinseq-lite.pl -min_len 10 -trim_qual_right 20 -min_qual_mean 20 -fastq R1_step1.fastq -fastq2 R2_step1.fastq -out_good clean	
# Cleaning intermediated files

	rm R1_step1.fastq
	rm R2_step1.fastq

	echo "Proceed to step 2 in 2_map directory"
		cd ../2_map		

# Mapping to reference genome
	echo 
	echo "Mapping using bwa"
	echo "Mapping against ARS-UCD1.2+BLV"
		bwa mem -t 4 -Y -L 0 -M -R "@RG\tID:sample\tSM:sample\tPL:Illumina" /mnt/nas/reference_genome/BWA/BLV/ARS+BLV/genome.fa ../1_qc/clean_1.fastq ../1_qc/clean_2.fastq > ARS-UCD1.2+virus.sam
	echo 
	echo "Mapping complete - proceed to extract reads"
	echo 
# Reads extraction and selection
	echo "Begin reads extraction and selection"
	echo 
		samtools view -Sb ARS-UCD1.2+virus.sam > ARS-UCD1.2+virus.bam
	echo 
	echo "Filter and accept only first mapped reads"
		samtools view -bh -F 256 -o ARS-UCD1.2+virus_uniq.bam ARS-UCD1.2+virus.bam 
		samtools sort ARS-UCD1.2+virus_uniq.bam -o ARS-UCD1.2+virus_uniq_sort.bam

	echo "Duplicates removal"
                java -jar /mnt/nas/programs/picard.jar MarkDuplicates INPUT=./ARS-UCD1.2+virus_uniq_sort.bam OUTPUT=./ARS-UCD1.2+virus_uniq_sort_duprmv.bam METRICS_FILE=marked_dup_metrics.txt REMOVE_DUPLICATES=TRUE MAX_RECORDS_IN_RAM=null ASSUME_SORTED=true TMP_DIR=./
	echo "Convert bam file to sam file"
		samtools view -h ARS-UCD1.2+virus_uniq_sort_duprmv.bam > ARS-UCD1.2+virus_uniq_sort_duprmv.sam
	echo 
	echo "Filter reads which map to BLV"
		awk -F"\t" '($3 ~ /BLV/ || $7 ~/BLV/ || $1~/^@/) {print}' ARS-UCD1.2+virus_uniq_sort_duprmv.sam > virus_uniq_map_duprmv.sam
	echo 
	echo "Extract reads with soft-clipping"
        	awk -F"\t" '($6 ~/S/ || $1~/^@/) {print}' virus_uniq_map_duprmv.sam > total_softclipping.sam
	echo "sam-to-bam"
		samtools view -Sb virus_uniq_map_duprmv.sam > virus_uniq_map_duprmv.bam

# Cleaning intermediated files
		rm ARS-UCD1.2+virus.sam
		rm ARS-UCD1.2+virus.bam
		rm ARS-UCD1.2+virus_uniq_sort.bam
		rm ARS-UCD1.2+virus_uniq.bam

# Unzip fastq file

		gzip ../R*
		
# Count number of reads 
	echo "Analysis complete - displaying read counts"

	echo "Mapped reads"
		samtools view -F 0x4 ../2_map/ARS-UCD1.2+virus_uniq_sort_duprmv.bam | cut -f 1 | sort | uniq | wc -l

	echo "Reads mapped to BLV"
		samtools view -F 0x4 ../2_map/virus_uniq_map_duprmv.bam | cut -f 1 | sort | uniq | wc -l

	echo "Soft-clipping Reads mapped to total viral reads before cleaning"
                samtools view -S -F 0x4 ../2_map/total_softclipping.sam | cut -f 1 | sort | uniq | wc -l

	 "Congratulations! Data processing complete!"
	echo "Task completed on"
		date
	end=$(date +%s)
	duration=$(($end-$begin))
	echo "Script run time : $(($duration / 3600)) hours $((($duration % 3600) / 60)) minutes $(($duration % 60)) seconds"
	echo
#
# ************ Main Script Ends Here ************ #
