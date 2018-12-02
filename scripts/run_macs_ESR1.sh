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
echo "*                  RUN MACS2                   *"
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
# RUN MACS2
###############################################################################
# WT with input control
echo `date +"%T"`" siNT samples"
echo " "
start_WT=`date +%s`
mkdir -p ${result_path}/siNT
cd ${result_path}/siNT
macs2 callpeak -t ${result_path}/siNT.merged.bam \
               -c ${result_path}/bowtie2_alignment/SRR540220.sorted.bam \
               -n siNT_vs_input_macs2 \
               -f BAM -g hs -B -q 0.01
end_WT=`date +%s`

# GATA with input control
echo " "
echo `date +"%T"`" siGATA samples"
echo " "
start_GATA=`date +%s`
mkdir -p ${result_path}/siGATA
cd ${result_path}/siGATA
macs2 callpeak -t ${result_path}/siGATA.merged.bam \
               -c ${result_path}/bowtie2_alignment/SRR540220.sorted.bam \
               -n siGATA_vs_input_macs2 \
               -f BAM -g hs -B -q 0.01
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
echo "    siNT samples: "`tiempo $runtime_WT`
echo "    siGATA samples: "`tiempo $runtime_GATA`
echo "        total time: "`tiempo $runtime`
