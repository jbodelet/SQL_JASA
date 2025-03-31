#!/bin/bash

# Activate conda environment
conda activate tensorflow

export CUDA_VISIBLE_DEVICES=-1

echo "Running all simulations in parallel"
python vae_additive.py -K 5 -q 1 -dir ../simulations/M1/data_n50_var1 &
python vae_additive.py -K 5 -q 1 -dir ../simulations/M1/data_n100_var1 &
python vae_additive.py -K 5 -q 1 -dir ../simulations/M1/data_n500_var1 &
python vae_additive.py -K 5 -q 1 -dir ../simulations/M1/data_n1000_var1 &

# ---------------------
# Run simulations for M2_q1
# ---------------------

python vae_additive.py -K 5 -q 1 -dir ../simulations/M2_q1/data_p20 &
python vae_additive.py -K 7 -q 1 -dir ../simulations/M2_q1/data_p50 &
python vae_additive.py -K 9 -q 1 -dir ../simulations/M2_q1/data_p100 &
python vae_additive.py -K 11 -q 1 -dir ../simulations/M2_q1/data_p200 &
python vae_additive.py -K 15 -q 1 -dir ../simulations/M2_q1/data_p500 &


python vae_additive.py -K 5 -q 3 -dir ../simulations/M2_q3/data_p20 &
python vae_additive.py -K 7 -q 3 -dir ../simulations/M2_q3/data_p50 &
python vae_additive.py -K 9 -q 3 -dir ../simulations/M2_q3/data_p100 &
python vae_additive.py -K 11 -q 3 -dir ../simulations/M2_q3/data_p200 &
python vae_additive.py -K 15 -q 3 -dir ../simulations/M2_q3/data_p500 &

