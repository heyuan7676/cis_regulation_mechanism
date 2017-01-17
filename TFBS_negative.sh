#!/bin/bash

# CELLTYPE=Mon
# eQTL_TYPE=eSNP

CELLTYPE=$1
eQTL_TYPE=$2
chr=$3


## get the fasta data
FA_DIR='/Users/Yuan/Documents/BLab/Predict_target_genes/data/assembly'
FA_FILE=GRCh37.p13.genome.fa 

INTERM_DIR='/Users/Yuan/Documents/BLab/Predict_target_genes/data/intermediate'/${CELLTYPE}
BED_FILE=PCHiC_peak_matrix_cutoff5_${chr}_${eQTL_TYPE}_negativeset.bed

bedtools getfasta -fi ${FA_DIR}/${FA_FILE} -bed ${INTERM_DIR}/${BED_FILE} -name -fo ${INTERM_DIR}/${BED_FILE%.bed}.fa

## scan for known motifs
## mkdir ${INTERM_DIR}/motif/

MOTIF_DIR='/Users/Yuan/Documents/tools/motif_databases'
MOTIF_F='JASPAR/JASPAR_CORE_2016_vertebrates.meme'

MOTIF_LIST_F=${INTERM_DIR}/motif_list.txt
MOTIF_LIST=`cat ${MOTIF_LIST_F}`

## from Corces, M. Ryan, et al. "Lineage-specific and single-cell chromatin accessibility charts human hematopoiesis and leukemia evolution." 
## Nature genetics (2016) Fig4A, Supplementary Table 3.

## GATA PAX5 FOS SPI1 EHF ELF1 RUNX ROR TAL1-TCF3 CTCF
## MA0482.1 MA0014.2 MA0477.1 MA0080.3 MA0598.1 MA0473.1 MA0002.2 MA0071.1 MA0091.1 MA0139.1
## MOTIF_list='MA0482.1 MA0014.2 MA0477.1 MA0080.3 MA0598.1 MA0473.1 MA0002.2 MA0071.1 MA0091.1 MA0139.1'



for mot in ${MOTIF_LIST}
do 
    fimo --motif ${mot} --qv-thresh --thresh 0.05 --oc ${INTERM_DIR}/motif/${mot}_${eQTL_TYPE}_negative ${MOTIF_DIR}/${MOTIF_F} ${INTERM_DIR}/${BED_FILE%.bed}.fa 2>> ${INTERM_DIR}/motif/fimo.output.${eQTL_TYPE}.txt
done


### Re-generate the bed files where TFB motif has been reported

OUTPUT_FILE=${INTERM_DIR}/motif/fimo.output.${CELLTYPE}.${eQTL_TYPE}.${chr}.negativeset.bed
rm $OUTPUT_FILE

for mot in ${MOTIF_LIST}
do
	cat ${INTERM_DIR}/motif/${mot}_${eQTL_TYPE}_negative/fimo.txt | awk '{print $1,$2}' | sed "/#/d" >> ${OUTPUT_FILE}
done

paste <(cut -d":" -f 1 ${OUTPUT_FILE}  | awk '{print $2}')  \
      <(cut -d":" -f 3 ${OUTPUT_FILE}) \
      <(cut -d":" -f 4 ${OUTPUT_FILE} | tr '-' '\t' )\
      <(cut -d" " -f 1 ${OUTPUT_FILE}) > temp

mv temp ${OUTPUT_FILE}


