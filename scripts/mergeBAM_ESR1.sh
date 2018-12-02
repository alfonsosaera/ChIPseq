#! /bin/bash

# Alfonso Saera Vila

# V 1.0
# 01/12/2018

###############################################################################
# TIMING the SCRIPT
###############################################################################
start=`date +%s`

echo " "
echo "************************************************"
echo "*               MERGING BAM FILES              *"
echo "************************************************"
echo " "
echo `date +"%T"`" Running file: $0"
echo " "

###############################################################################
# Functions
###############################################################################
function tiempo {
  min=$(( $1 / 60 ))
  seconds=$(( $1 - $(( $min * 60 )) ))
  minutes=$min
  hours=0
  if [ "$min" -gt 60 ]
  then
    hours=$(( $min / 60 ))
    minutes=$(( $min - $(( $hours * 60 )) ))
  fi
  echo $hours" h., "$minutes" min., "$seconds" secs."
}

###############################################################################
# set variables
###############################################################################
Working_dir_path="./"
result_path="${Working_dir_path}/results"

###############################################################################
# MERGE BAM files
###############################################################################
# WT samples: SRR540188 SRR540190 SRR540192
echo `date +"%T"`" WT BAM files"
echo " "
start_WT=`date +%s`
samtools merge -@ 2 ${result_path}/siNT.merged.bam \
               ${result_path}/bowtie2_alignment/SRR540188.sorted.bam \
               ${result_path}/bowtie2_alignment/SRR540190.sorted.bam\
               ${result_path}/bowtie2_alignment/SRR540192.sorted.bam
samtools index -@ 2 ${result_path}/siNT.merged.bam
end_WT=`date +%s`

# GATA samples: SRR540189 SRR540191 SRR540193
echo " "
echo `date +"%T"`" GATA BAM files"
echo " "
start_GATA=`date +%s`
samtools merge -@ 2 ${result_path}/siGATA.merged.bam \
               ${result_path}/bowtie2_alignment/SRR540189.sorted.bam \
               ${result_path}/bowtie2_alignment/SRR540191.sorted.bam\
               ${result_path}/bowtie2_alignment/SRR540193.sorted.bam
samtools index -@ 2 ${result_path}/siGATA.merged.bam
end_GATA=`date +%s`

###############################################################################
# TIMING the SCRIPT
###############################################################################
runtime_WT=$(( end_WT - start_WT ))
runtime_GATA=$(( end_GATA - start_GATA ))
end=`date +%s`
runtime=$(( end - start ))

# show timing in terminal
echo " "
echo `date +"%T"`" $0 runtime:"
echo "    WT files merging: "`tiempo $runtime_WT`
echo "    GATA files merging: "`tiempo $runtime_GATA`
echo "        total time: "`tiempo $runtime`
