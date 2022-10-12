#!/bin/bash

####################################################
#                                                  # 
#           2_map_BLV+LTR+noLTR          	   #
#                                                  #
####################################################
#
#
# 
# ************ Main Script Begins Here ************ #

 	echo
        echo "Begin analysis"
        echo "Task started at"
                date
        begin=$(date +%s)

# make analysis directories
	echo
	echo "Create analysis directories"
		mkdir 3_virus_host_junction
	echo "Directories  3_virus_host_junction created"

	echo "Proceed to mapping step in 3_virus_host_junction directory"
		cd ./3_virus_host_junction

# Mapping to reference genome (ARS-UCD1.2+AB513134+LTR)
	echo
	echo "Mapping using bwa"
	echo "Mapping against ARS-UCD1.2+BLV+LTR"
		bwa mem -t 4 -M -R "@RG\tID:sample\tSM:sample\tPL:Illumina" /mnt/nas/reference_genome/BWA/BLV/ARS-UCD1.2+BLV+LTR_2/genome.fa ../1_qc/clean_1.fastq ../1_qc/clean_2.fastq > ARS-UCD1.2+virus+LTR.sam
	echo
	echo "Mapping complete - proceed to extract reads"
	echo

# Reads extraction and selection
	echo "Begin reads extraction and selection"
	echo
		samtools view -Sb ARS-UCD1.2+virus+LTR.sam > ARS-UCD1.2+virus+LTR.bam
	echo
	echo "Filter and accept only first mapped reads"
		samtools view -bh -F 256 -o ARS-UCD1.2+virus+LTR_uniq.bam ARS-UCD1.2+virus+LTR.bam
	echo
		samtools sort ARS-UCD1.2+virus+LTR_uniq.bam > ARS-UCD1.2+virus+LTR_uniq_sort.bam
	echo
     	echo "Duplicates removal"
                java -jar /mnt/nas/programs/picard.jar MarkDuplicates INPUT=./ARS-UCD1.2+virus+LTR_uniq_sort.bam OUTPUT=./ARS-UCD1.2+virus+LTR_uniq_sort_duprmv.bam METRICS_FILE=marked_dup_metrics.txt TMP_DIR=./tmp REMOVE_DUPLICATES=TRUE MAX_RECORDS_IN_RAM=null ASSUME_SORTED=true

	echo "Convert bam to sam file"
		samtools view -h ARS-UCD1.2+virus+LTR_uniq_sort_duprmv.bam > ARS-UCD1.2+virus+LTR_uniq_sort_duprmv.sam
	echo
	echo "Filter reads containing BLV"

awk -F "\t" '($3 ~ /chrBLV_LTR/ || $7 ~/chrBLV_LTR/ || $3 ~/chrBLV_noLTR/ || $7 ~/chrBLV_noLTR/ || $1~/^@/) {print}' ARS-UCD1.2+virus+LTR_uniq_sort_duprmv.sam > virus+LTR_duprmv_uniq_sort.sam

# Remove imtermediate reads 
        echo
        echo
# making total_chimera_paired_sort.bam file

awk -F "\t" '($3 ~ /chrBLV_LTR/ && $7 ~/chr/ || $3 ~ /chr/ && $7 ~/chrBLV_LTR/ || $3 ~ /chrBLV_noLTR/ && $7 ~/chr/ || $3 ~ /chr/ && $7 ~/chrBLV_noLTR/ || $1~/^@/) {print}' virus+LTR_duprmv_uniq_sort.sam > total_chimera_paired.sam

		samtools view -Sb total_chimera_paired.sam > total_chimera_paired.bam

# IS_analysis by using python
        echo "activate python2.7"
               source activate py27_bioinfo
        echo
	echo "perform IS analysis"
        

          python2.7 /mnt/nas/previous_members/no/shell_scripts/detect_Is.py total_chimera_paired.bam /mnt/nas/previous_members/no/shell_scripts/ARS-UCD1.2+BLV_LTR+noLTR_chrome_size.txt 1000 20


# extract Total_IS.sam
	  python2.7 /mnt/nas/previous_members/no/shell_scripts/grep_reads.py total_chimera_paired.bam ./total_chimera_paired.bam.out.txt ./
	echo
	awk -F"\t" ' $1~/^@/ {print}' total_chimera_paired.sam > header
	echo
		cat header IS* > total_IS.sam 
	echo	
		samtools view -Sb total_IS.sam > total_IS.bam
		samtools sort total_IS.bam -o total_IS_sort.bam
		samtools view -h total_IS_sort.bam > total_IS_sort.sam

# remove intermediated files
		rm ARS-UCD1.2+virus+LTR.*
		rm ARS-UCD1.2+virus+LTR_uniq.bam
		rm ARS-UCD1.2+virus+LTR_uniq_sort.bam
		rm -r tmp
		rm ARS-UCD1.2+virus+LTR_uniq_sort_duprmv.sam
		rm total_IS.*
		rm total_IS_sort.bam

	echo "Congratulations! Data processing complete!"
	echo "Task completed on"
		date
	end=$(date +%s)
	duration=$(($end-$begin))
 	echo "Script run time : $(($duration / 3600)) hours $((($duration % 3600) / 60)) minutes $(($duration % 60)) seconds"
	echo
#
# ************ Main Script Ends Here ************ #

