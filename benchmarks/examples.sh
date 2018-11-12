#!/bin/bash

# See the README before running this script.
# Note: Each benchmark is commented out by default.
# Note: Relative file paths are used.

# add stochkit to path

export PATH="$PATH:$STOCHKIT"

##### BioSimulator.jl Benchmarks #####

##### serial #####

# julia benchmarks_biosimulator.jl kendall mmek autoreg dimer-decay yeast

##### parallel #####

export JULIA_NUM_THREADS=8

# julia benchmarks_biosimulator.jl kendall mmek autoreg dimer-decay yeast

export JULIA_NUM_THREADS=1

##### Gillespie.jl Benchmarks #####

# julia benchmarks_gillespie.jl kendall mmek autoreg dimer-decay yeast

##### StochPy Benchmarks #####

# python3 stochpy_bench.py kendall       4.0 1000 10 5357 1000
# python3 stochpy_bench.py mmek         50.0 1000 10 5357 1000
# python3 stochpy_bench.py dimer-decay 100.0 1000 10 5357 1000
# python3 stochpy_bench.py yeast       100.0 1000 10 5357 1000

# it takes approximately 25 hours to run 1000 samples
# 100 samples is good enough for our purposes

# python3 stochpy_bench.py autoreg     500.0 1000 10 5357 100

##### StochKit Benchmarks #####

# julia benchmarks_stochkit.jl kendall mmek autoreg dimer-decay yeast

##### summarize #####

# julia analysis.jl biosimulator-serial.json   >> summary-biosimulator-serial.txt
# julia analysis.jl biosimulator-parallel.json >> summary-biosimulator-parallel.txt
# julia analysis.jl gillespie.json             >> summary-gillespie.txt
# julia analysis.jl stochkit.json              >> summary-stochkit.txt

# julia analysis.jl ./stochpy/kendall.txt     >> summary-stochpy.txt
# julia analysis.jl ./stochpy/mmek.txt        >> summary-stochpy.txt
# julia analysis.jl ./stochpy/autoreg.txt     >> summary-stochpy.txt
# julia analysis.jl ./stochpy/dimer-decay.txt >> summary-stochpy.txt
# julia analysis.jl ./stochpy/yeast.txt       >> summary-stochpy.txt
