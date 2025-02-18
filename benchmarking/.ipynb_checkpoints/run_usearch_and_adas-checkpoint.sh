#!/bin/bash

# -----------------------------
# Author: Sherlyn Weng
# Last Updated: 2025-02-17
# This is a script to randomly sample one 16s sequence from greengenes2 
# and run adas-search and usearch to search for similar sequences. 
# This procedure is repeated 100 times.
# This script was edited to remove slurm headers and absolute paths
# -----------------------------

#### STEP 0: Build Adas Database ####

# cd ~/adas_db
# adas-build --input gg2_reference.fna \
#     --max_nb_connection 255 \
#     --sketch-size 256 \
#     --hnsw-ef 800 \
#     -t 16

cpu=32
gg2_fna="gg2_reference.fna"

for i in $(seq 1 100); do
    echo "========== Iteration $i =========="

    #### STEP 1: Sample sequences ####
    echo "Sampling sequences ..."

    # Create a unique directory for this iteration's test data
    TEST_DATA_DIR="~/adas_test/test_data_$i"
    mkdir -p "$TEST_DATA_DIR"
    cd "$TEST_DATA_DIR" || exit 1

    conda activate seqtk
    seqtk sample ~/gg2_test.fna 1 -s $i > test.fa
    conda deactivate

    echo "Sampling sequences done for iteration $i"

    #### STEP 2.1: Search using adas ####
    echo
    echo "Searching using adas ..."

    # Create a unique directory for this iteration's adas results
    ADAS_DIR="$TEST_DATA_DIR/adas"
    mkdir -p "$ADAS_DIR"
    cd "$ADAS_DIR" || exit 1

    # Assuming adas database is built and available at ../adas_db relative to the base directory
    /usr/bin/time -v /home/y1weng/packages-code/ADAS/adas-search --input "$TEST_DATA_DIR/test.fa" \
                             -n 256  \
                             -b /home/y1weng/55_test_adas/adas_db \
                             -t $cpu

    echo "Searching using adas done for iteration $i"
    
    #### STEP 2.2: Chainning using adas ####
    echo
    echo "Extracting reference sequences from adas search results ..."
    
    awk '{print $7}' "$ADAS_DIR/adas.neighbors.txt" | tail -n +3 > "$ADAS_DIR/ids.txt"
    
    conda activate /home/y1weng/mambaforge/envs/wipe_dev3
    seqkit grep -nrf "$ADAS_DIR/ids.txt" $gg2_fna > "$ADAS_DIR/chain_ref.fna"
    conda deactivate
    
    /usr/bin/time -v /home/y1weng/packages-code/ADAS/adas-chain --query "$TEST_DATA_DIR/test.fa" \
        --reference "$ADAS_DIR/chain_ref.fna" \
        --output "$ADAS_DIR/adas.chain.results.txt"
    
    #### STEP 3: Search using usearch ####
    echo
    echo "Searching using usearch ..."

    # Create a unique directory for this iteration's usearch results
    USEARCH_DIR="$TEST_DATA_DIR/usearch"
    mkdir -p "$USEARCH_DIR"
    cd "$USEARCH_DIR" || exit 1

    conda activate /home/y1weng/mambaforge/envs/usearch
    /usr/bin/time -v usearch -usearch_global "$TEST_DATA_DIR/test.fa" \
            -db /home/y1weng/55_test_adas/greengenes2/gg2_reference.fna \
            -id 0.1 \
            -blast6out hits.b6 \
            -strand both \
            -maxaccepts 0 \
            -maxrejects 0 \
            -threads $cpu 
    conda deactivate

    echo "Searching using usearch done for iteration $i"
    echo "========== End of iteration $i =========="
done

conda deactivate