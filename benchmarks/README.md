# Benchmarking

This folder provides scripts to compare the run times of two Vistle modules automatically. While they have been developed to compare the performance of an original Vistle module with its Viskores counterpart (recognizable by the "Vtkm" suffix in its name), they can be used to compare any two modules.

This guide explains how to use the scripts and how to prepare a visualization workflow for benchmarking.

## Using the Benchmarking Scripts

In order to use them, the scripts have to be made executable first (this only has to be done once). Navigate to your Vistle root directory and run the following commands:
```bash
$ cd benchmarks
$ chmod +x benchmark_map.sh run_benchmarks.sh
```

### Configuring the Scripts

`run_benchmarks.sh` is the main script for benchmarking and the only one you will have to configure. Open it in your favorite IDE. In the following we will go through its settings one by one.

```bash
export VISTLE_CPU_MAP=$PWD/benchmark_iso_cpu.vsl
export VISTLE_GPU_MAP=$PWD/benchmark_iso_gpu.vsl
```

`VISTLE_CPU_MAP` and `VISTLE_GPU_MAP` are the two visualization workflows that will be benchmarked. You can either use the example workflows provided in this folder or change the paths to point to your own workflows. The example workflows measure the run time of the `IsoSurface` and `IsoSurfaceVtkm` module (for the arbitrary isovalue 1.1), respectively.

```bash
export ITERATIONS=10
export GRID_SIZE=100
```

With `ITERATIONS` you can set how many times each workflow shall be run. `GRID_SIZE` lets you specify the size of the input grid. In the default workflows, an unstructured grid of size `GRID_SIZE` x `GRID_SIZE` x `GRID_SIZE` is created using `Gendat` and then passed to the module computing the isosurface.

```bash
# Make sure Vistle was compiled with CMAKE_BUILD_TYPE=Release!
export VISTLE_BUILD_TO_USE=$HOME/Software/vistle/build_gpu_single_release
```

Finally, make sure to set `VISTLE_BUILD_TO_USE` to the build directory of the Vistle installation you want to use for benchmarking. Make sure that you compile Vistle with "-DCMAKE_BUILD_TYPE=Release".

The runtimes and log files the workflow creates will be stored in an output folder inside the `benchmarking` called `<prefix>-n<GRID_SIZE>-it<ITERATIONS>`. Notice that the prefix is set when calling `benchmark.sh` at the end of the script:

```bash
# --- Running the benchmarks --- 
./benchmark_map.sh $VISTLE_CPU_MAP cpu
./benchmark_map.sh $VISTLE_GPU_MAP gpu
```
 By default, the prefixes are `cpu` and `gpu`, but you can of course change this to any string you want. In this example, the output folders will be called `cpu-n100-it10` and `gpu-n100-it10`.

### Running the Script

Now, you're all set to run the script! Simply execute the following inside the `benchmarking` folder:

```
$ ./run_benchmarks.sh
```

This will run each workflow `ITERATIONS` times, automatically extract all runtimes from the log files and compute the average runtime. The runtimes and the average runtime are stored in the output file `times.txt`. You can, e.g., open it like this in a terminal:

```
$ cat cpu-n100-it10/times.txt

Average: 0.124009s
Run 0: 0.1203s
Run 1: 0.119777s
Run 2: 0.132077s
Run 3: 0.121957s
Run 4: 0.121428s
Run 5: 0.120738s
Run 6: 0.130636s
Run 7: 0.128955s
Run 8: 0.122896s
Run 9: 0.121326s
```

## How to Benchmark Custom Visualization Workflows

Create a visualization workflow as usual using the Vistle GUI. You can, e.g., open it like this from inside the vistle root directory:
```bash
$ vistle
```

Once you are happy with the workflow, click on the Vistle module you would like to benchmark. This will open the module parameter menu on the right. Open the `System Parameters` tab at the top of the parameter module and tick `benchmark`. Now, save your vistle worklow as `.vsl` file.

Open your `.vsl` file in your favorite IDE (it is essentially a python script) and add the following to the end of your file:

```python
barrier()
compute()
barrier()
quit()
```

This will execute your visualization workflow once and then close it. This makes it possible to run the workflow automatically without the GUI or any user input in the benchmark scripts.

If you want to use the `GRID_SIZE` environment variable from `run_benchmarks.sh` in your visualization workflow, e.g., to change the size of your input grid automatically, add the following at the top of your `.vsl` file:

```python
import os

n = int(os.getenv("GRID_SIZE"))
```

You can then use it anywhere in your `.vsl`. Like here to create an unstructured grid of size `GRID_SIZE` x `GRID_SIZE` x `GRID_SIZE`:

```python
mGendat1 = waitForSpawn(umGendat1)
setVectorParam(mGendat1, '_position', 88.36500549316406, -281.88995361328125, True)
setIntParam(mGendat1, 'size_x', n, True)
setIntParam(mGendat1, 'size_y', n, True)
setIntParam(mGendat1, 'size_z', n, True)
```

Feel free to create your own environment variables in `run_benchmark.sh` and use them in your `.vsl` file!

Save your `.vsl` file and change the path accordingly in `run_benchmark.sh`. Of course, there's nothing stopping you from calling `./benchmark_map.sh` on more than two visualization workflows in the script. Just adjust it however you like!
