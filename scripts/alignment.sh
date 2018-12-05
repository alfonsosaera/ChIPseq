#! /bin/bash

# Alfonso Saera Vila

# V 1.0
# 27/11/2018

###############################################################################
# TIMING the SCRIPT
###############################################################################
start=`date +%s`

echo " "
echo "************************************************"
echo "*   BOWTIE 2 ALIGMENT OF TRIMMED FASTQ FILES   *"
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
bowtie2_path="/mnt/518D6BCF3ECC578E/bowtie2-2.3.4.3-linux-x86_64"
index_files="/mnt/518D6BCF3ECC578E/indexes/Arabidopsis_thaliana_Ensembl_TAIR10/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome"
Working_dir_path="/mnt/518D6BCF3ECC578E/ChIP-seq"
result_path="${Working_dir_path}/results"
TR_path="${result_path}/tr_fastq"

###############################################################################
# make list of sample names
###############################################################################
declare -A sample_names

for file in $TR_path/*.fastq
do
  sample=${file##*/}
  name=${sample%_*.fastq}
  sample_names["$name"]=1
done

###############################################################################
# Process samples
###############################################################################
for sample in "${!sample_names[@]}"
do
  echo `date +"%T"`" - Processing sample $sample"

  # create result folder
  mkdir -p $result_path/$sample/bowtie2_alignment

  # align with bowtie2
  cd $result_path/$sample/bowtie2_alignment
  $bowtie2_path/bowtie2 -p 3 \
                        -x $index_files \
                        -1 $TR_path/${sample}_1.fastq \
                        -2 $TR_path/${sample}_2.fastq \
                        -S ${sample}.sam
  # improvement
  # bowtie2 (...parameters...) | samtools sort -o sorted.bam

  # convert to BAM, sort and index
  samtools view --threads 2 -bS ${sample}.sam > ${sample}.bam
  samtools sort --threads 2 ${sample}.bam -o ${sample}.sorted.bam
  samtools index -@ 2 ${sample}.sorted.bam
done

###############################################################################
# TIMING the SCRIPT
###############################################################################
end=`date +%s`
runtime=$(( end - start ))

# show timing in terminal
echo " "
echo $0" runtime:"
echo "    "`tiempo $runtime`
