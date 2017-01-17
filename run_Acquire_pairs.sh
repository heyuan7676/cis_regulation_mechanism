#!/bin/bash -l


# =============================

match_distance=2000000
SNP_window_for_PHiC=5000
promoter_window_for_PHiC=5000
chr=$1
echo $match_distance
echo $SNP_window_for_PHiC
echo $promoter_window_for_PHiC

export chr
export match_distance
export SNP_window_for_PHiC
export promoter_window_for_PHiC


from_beginning=False
compute_negative_set=False
negative_set_flag=False

export from_beginning
export compute_negative_set
export negative_set_flag

python Acquire_pairs_mon.py


# ====== 


compute_negative_set=False
negative_set_flag=True

export compute_negative_set
export negative_set_flag

python Acquire_pairs_mon.py



