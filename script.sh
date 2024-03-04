#!/bin/bash

# Define values for n and p
n_values=(4 16 100 10000 50176 100489 1000000)
p_values=(4 16 100 10000 50176 100489 1000000)

# Create directories ex1 to ex7
for i in {1..7}; do
    mkdir -p ex$i
    cp a.out ex$i/
    cd ex$i || exit
    # Loop through n and p values simultaneously
    for ((j=0; j<${#n_values[@]}; j++)); do
        n=${n_values[$j]}
        p=${p_values[$j]}
        # Run the command
        mpirun -np $p tau_exec a.out $n --oversubscribe
    done
    cd ..
done
