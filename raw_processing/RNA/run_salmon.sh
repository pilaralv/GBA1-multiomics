#!/bin/sh
module load salmon

# Main variables
cd $MAIN/IPSC_PROCESSED_DATA/salmon/

SAMPLE=$1
OUT_NAME=$2
salmon quant -i salmon_index/ -l A -1 $MAIN/IPSC_RAW_DATA/${SAMPLE}_1.fastq.gz -2 $MAIN/IPSC_RAW_DATA/${SAMPLE}_2.fastq.gz --validateMappings -o ${OUT_NAME}

