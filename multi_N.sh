#!/bin/zsh

# check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <a> <N_step> <plot>"
    exit 1
fi

# assign command-line arguments to variables
a=$1
N_step=$2
plot=$3

# compile
make cn

# Clear the frequencies file
echo "N,frequency,error" > data/trapped/frequencies.csv

# run with different parameters
for N in $(seq 99 $N_step 299); do
    echo -ne "Running simulation with N = $N\r"
    ./cn.x $N 0.001 0.0 $a
    python3 plots/probability.py $plot N
done

# plot the frequencies as a function of a
python3 plots/frequencies.py
