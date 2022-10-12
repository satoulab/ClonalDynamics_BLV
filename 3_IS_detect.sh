#!/bin/bash

# ************ Main Script Begins Here ************ #

 	echo
        echo "Begin analysis"
        echo "Task started at"
                date
        begin=$(date +%s)

# make analysis directories
	echo
	echo "Create analysis directories"
		mkdir 4_virus_host_junction
	echo "Directories  4_virus_host_junction created"

    echo "Bam file copy"
        cp ./3_virus_host_junction/total_chimera_paired.bam ./4_virus_host_junction/total_chimera_paired.bam

	echo "Proceed to mapping step in 4_virus_host_junction directory"
		cd ./4_virus_host_junction



# IS_analysis by using python
        echo "activate python2.7"
               source activate py27_bioinfo
        echo
	echo "perform IS analysis"
        

          python2.7 /mnt/nas/previous_members/no/shell_scripts/detect_Is.py total_chimera_paired.bam /mnt/nas/previous_members/no/shell_scripts/ARS-UCD1.2+BLV_LTR+noLTR_chrome_size.txt 1000 0


# extract Total_IS.sam
	  python2.7 /mnt/nas/previous_members/no/shell_scripts/grep_reads.py total_chimera_paired.bam ./total_chimera_paired.bam.out.txt ./

	echo "Congratulations! Data processing complete!"
	echo "Task completed on"
		date
	end=$(date +%s)
	duration=$(($end-$begin))
 	echo "Script run time : $(($duration / 3600)) hours $((($duration % 3600) / 60)) minutes $(($duration % 60)) seconds"
	echo
#
# ************ Main Script Ends Here ************ #
