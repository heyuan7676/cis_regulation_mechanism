#!/bin/bash -l
#SBATCH --time 50:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --no-requeue

# =============================


chr=$1

compute_negative_set=True
match_distance=2000000

compute_contacting_degree=True
compute_chromosome_pos_LD=True
SNP_window_for_motif=2000
promoter_window_for_motif=2000

merge_ATAC_data=True
SNP_ATAC_window=2000
gene_ATAC_window=10000


merge_TF_motifs=False



echo chr$chr

export chr
export celltype

export compute_negative_set
export match_distance

export compute_contacting_degree
export compute_chromosome_pos_LD
export SNP_window_for_motif
export promoter_window_for_motif


export merge_ATAC_data
export SNP_ATAC_window
export gene_ATAC_window

export merge_TF_motifs

python Acquire_pairs_2.py
