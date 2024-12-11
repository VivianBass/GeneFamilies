#!/bin/bash
#$ -l mem_free=20G,h_vmem=20G       
#$ -pe smp 1

echo "**** Job starts ****"              
date

echo "**** info ****"                    
echo "User: ${USER}"                     
echo "Job id: ${JOB_ID}"                 
echo "Job name: ${JOB_NAME}"             
echo "Hostname: ${HOSTNAME}"             
echo "Reads Dir: READS_DIR"               

/mnt/bin/kallisto/v0.46.1/kallisto quant \
  -i REFERENCE_TRANSCRIPTOME_INDEX \
  -o OUT_DIR FRWRD BCKWRD

echo "**** Job ends ****"
date
