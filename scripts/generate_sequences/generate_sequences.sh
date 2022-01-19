#!/bin/bash

PATH_TO_DATA=$1
PATH_TO_OUTDIR=$2
IFS=',' read -ra SEQ_LENGTHS <<< $3
PATH_TO_SEQGEN=$4
MODEL="GTR"
SEED=100

for filename in $PATH_TO_DATA/*/; do

  # Set tree and parameter files
  TREE_ID=$(basename $filename)
  TREE_FILE="${filename}tree_best.newick"
  MODEL_FILE="${filename}log_0.txt"

  # Get model parameters
  while read line; do
    IFS=':' read -ra ADDR <<< $line
      if [[ ${ADDR[0]} == "rate A <-> C" ]]; then
        RATE_A_C=$(echo ${ADDR[1]} | xargs)
      elif [[ ${ADDR[0]} == "rate A <-> G" ]]; then
        RATE_A_G=$(echo ${ADDR[1]} | xargs)
      elif [[ ${ADDR[0]} == "rate A <-> T" ]]; then
        RATE_A_T=$(echo ${ADDR[1]} | xargs)
      elif [[ ${ADDR[0]} == "rate C <-> G" ]]; then
        RATE_C_G=$(echo ${ADDR[1]} | xargs)
      elif [[ ${ADDR[0]} == "rate C <-> T" ]]; then
        RATE_C_T=$(echo ${ADDR[1]} | xargs)
      elif [[ ${ADDR[0]} == "rate G <-> T" ]]; then
        RATE_G_T=$(echo ${ADDR[1]} | xargs)
      elif [[ ${ADDR[0]} == "freq pi(A)" ]]; then
        PI_A=$(echo ${ADDR[1]} | xargs)
      elif [[ ${ADDR[0]} == "freq pi(C)" ]]; then
        PI_C=$(echo ${ADDR[1]} | xargs)
      elif [[ ${ADDR[0]} == "freq pi(G)" ]]; then
        PI_G=$(echo ${ADDR[1]} | xargs)
      elif [[ ${ADDR[0]} == "freq pi(T)" ]]; then
        PI_T=$(echo ${ADDR[1]} | xargs)
      fi
  done < $MODEL_FILE

  # Generate and run seq-gen command
  SEQ_PATH="${PATH_TO_OUTDIR}${TREE_ID}"
  mkdir -p ${SEQ_PATH}
  echo "GENERATE SEQUENCES FOR MODEL ${TREE_ID} OF SIZES $3 BP"
  echo "GENERATED SEQUENCES CAN BE FOUND UNDER ${SEQ_PATH}"
  echo ""
  for SEQ_LEN in "${SEQ_LENGTHS[@]}"; do
    OUTFILE="${SEQ_PATH}/seq_${SEQ_LEN}.fasta"
    GENERATE_SEQ="$PATH_TO_SEQGEN -q -z $SEED -of -m $MODEL -l $SEQ_LEN -f $PI_A,$PI_C,$PI_G,$PI_T -r $RATE_A_C,$RATE_A_G,$RATE_A_T,$RATE_C_G,$RATE_C_T,$RATE_G_T < $TREE_FILE > ${OUTFILE}"
    eval $GENERATE_SEQ
  done
done
#./scripts/generate_sequences/generate_sequences.sh ./data/trees_mini ./data/sequences_mini/ 100,1000,10000 ../Seq-Gen/source/seq-gen




