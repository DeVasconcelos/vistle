# How to Write a GPU Module

In this guide you will learn how to write Vistle modules that can be run on the GPU.    
Vistle makes use of the portable toolkit [VTK-m](https://m.vtk.org/) which allows running scientific visualization algorithms on various devices, including GPUs, and is designed to keep data transfers between devices at a minimum.    

## Overview
- [The VtkmModule Class](#the-vtkmmodule-class)
- [Example 1: A Basic GPU Module](#example-1-a-basic-gpu-module)
- [Example 2: Extending the Core Functionality](#example-2-extending-the-core-functionality)
- [Example 3: A Custom VTK-m Filter](#example-3-a-custom-vtk-m-filter)
- [How to Compile Vistle to Run Code on NVIDIA GPUs](#how-to-compile-vistle-to-run-code-on-nvidia-gpus)


## The VtkmModule Class

The `VtkmModule` class is the base class for VTK-m modules in Vistle. It is designed to make adding new VTK-m algorithms, so-called [filters](https://vtk-m.readthedocs.io/en/v2.2.0/provided-filters.html), as simple as possible by providing core functionality for handling the input data, passing it to the filter implemented by the child class, and writing the filter result to the output ports. At the same time, it is meant to be flexible, allowing the child class to customize and extend these processes, if desired.

### The Constructor

By default, a VTK-m module consists of one input port and one output port, but additional ports can be added by specifying the desired number in the base constructor (`numPorts`). Note that all data on the input ports must be defined on the same grid, as the module will throw an error, otherwise.

Many VTK-m filters work on data fields, so, by default, a VTK-m module expects the data at the input ports to contain mapped data in addition to a grid and will throw an error if there is none. If desired, the child class can remove this requirement by setting `requireMappedData` to `false`.

```cpp
VtkmModule(const std::string &name, int moduleID, mpi::communicator comm, int numPorts = 1, bool requireMappedData = true);
```

### The `setUpFilter` method


### Preparing the input data

### Preparing the output data

## Example 1: A Basic GPU Module
- `setRunFilter`
- add module paramters
- making changes to the construtor

## Example 2: Extending the Core Functionality
- overriding prepareInput and prepareOutput methods

## Example 3: A Custom VTK-m Filter

To learn more about using VTK-m filters check out [VTK-m's user guide](https://vtk-m.readthedocs.io/en/v2.2.0/basic-filter-impl.html).

## How to Compile Vistle to Run Code on NVIDIA GPUs