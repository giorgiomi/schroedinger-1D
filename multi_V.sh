#!/bin/zsh

# check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <a> <V_step> <plot>"
    exit 1
fi

# assign command-line arguments to variables
a=$1
V_step=$2
plot=$3

# compile
make cn

# Clear the frequencies file
echo "V,frequency,error" > data/trapped/frequencies.csv

# run with different parameters
for V0 in $(seq 0.0 $V_step 50000.0); do
    echo -ne "Running simulation with V0 = $V0 and a = $a\r"
    ./cn.x $V0 $a
    python3 plots/probability.py $plot V0
done


# plot the frequencies as a function of a
python3 plots/frequencies.py
