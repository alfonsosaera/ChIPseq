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
index_files="/mnt/518D6BCF3ECC578E/indexes/Homo_sapiens_NCBI_GRCh38/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/genome"
Working_dir_path="./"
result_path="${Working_dir_path}/results"
TR_path="${result_path}/tr_fastq"

###############################################################################
# make list of sample names
###############################################################################
declare -A sample_names

for file in $TR_path/*.fastq
do
  sample=${file##*/}
  name=${sample%.trim.fastq}
  sample_names["$name"]=1
done

###############################################################################
# Process samples
###############################################################################
for sample in "${!sample_names[@]}"
do
  echo `date +"%T"`" - Processing sample $sample"
  start_sample=`date +%s`

  # create result folder
  mkdir -p $result_path/bowtie2_alignment

  # align with bowtie2
  cd $result_path/bowtie2_alignment
  echo `date +"%T"`" - bowtie2_alignment step"
  $bowtie2_path/bowtie2 -p 3 \
                        -x $index_files \
                        -U $TR_path/${sample}.trim.fastq | \
                        samtools sort --threads 2 -o ${sample}.sorted.bam

  # index bam file
  echo `date +"%T"`" - indexing bam file"
  samtools index -@ 2 ${sample}.sorted.bam

  # sample processing time
  end_sample=`date +%s`
  runtime_sample=$(( end_sample - start_sample ))
  echo " "
  echo ${sample}" processing time:"
  echo "    "`tiempo $runtime_sample`
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
