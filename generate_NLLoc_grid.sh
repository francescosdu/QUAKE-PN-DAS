#!/bin/bash
echo "Warning! Generating NLLoc grid requires ~90GB. Do you want to generate grids? (Y/n)"

read -r answer

if [[ "$answer" =~ ^[Yy]$ || -z "$answer" ]]; then
    echo "Generating grids..."
    mkdir ./NLLoc/model_punto_medio
    mkdir ./NLLoc/time_punto_medio
    ./NLLoc/src/Vel2Grid ./NLLoc/run_punto_medio/ISNet_P.run
    ./NLLoc/src/Grid2Time ./NLLoc/run_punto_medio/ISNet_P.run
    ./NLLoc/src/Vel2Grid ./NLLoc/run_punto_medio/ISNet_S.run
    ./NLLoc/src/Grid2Time ./NLLoc/run_punto_medio/ISNet_S.run
else
    echo "Skipping grid generation."
fi
