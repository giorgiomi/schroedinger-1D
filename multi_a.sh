#!/bin/zsh

# check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <V0> <a_step> <plot>"
    exit 1
fi

# assign command-line arguments to variables
V0=$1
a_step=$2
plot=$3

# compile
make cn

# Clear the frequencies file
echo "a,frequency,error" > data/trapped/frequencies.csv

# run with different parameters
for a in $(seq 0.1 $a_step 0.9); do
    echo -ne "Running simulation with V0 = $V0 and a = $a\r"
    ./cn.x 199 0.001 $V0 $a
    python3 plots/probability.py $plot a
done


# plot the frequencies as a function of a
python3 plots/frequencies.py
