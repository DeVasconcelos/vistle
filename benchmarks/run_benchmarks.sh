#!/bin/bash

# --- Benchmark Setup ---
export VISTLE_CPU_MAP=$PWD/benchmark_iso_cpu.vsl
export VISTLE_GPU_MAP=$PWD/benchmark_iso_gpu.vsl

export ITERATIONS=5
export GRID_SIZE=10

# Make sure Vistle was compiled with CMAKE_BUILD_TYPE=Release!
export VISTLE_BUILD_TO_USE=$HOME/Software/vistle/build_gpu_single_release

# --- Running the benchmarks ---
# FIXME: When running multiple times, not all runs contain benchmark time... 
./benchmark_map.sh $VISTLE_CPU_MAP cpu

# FIXME: GPU benchmark crashes...
#./benchmark_map.sh $VISTLE_GPU_MAP gpu
