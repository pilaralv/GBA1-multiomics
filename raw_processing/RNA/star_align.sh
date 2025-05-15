#!/bin/sh

module load STAR
INDEX=$1
READ1=$2
READ2=$3
PREFIX=$4
STAR --genomeDir ${INDEX} --readFilesIn ${READ1} ${READ2} --readFilesCommand gunzip -c --outFileNamePrefix ${PREFIX} --outSAMunmapped None --outSAMtype BAM SortedByCoordinate 
