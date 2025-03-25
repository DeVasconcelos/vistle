# How to Write a GPU Module

In this guide you will learn how to write Vistle modules that can be run on the GPU.    
Vistle makes use of the portable toolkit [VTK-m](https://m.vtk.org/) which allows running scientific visualization algorithms on various devices, including GPUs, and is designed to keep data transfers between devices at a minimum.    

## Overview
- [The VtkmModule Class](#the-vtkmmodule-class)
- [Example 1: A Simple GPU Module](#example-1-a-simple-gpu-module)
- [Example 2: A More Advanced GPU Module](#example-2-a-more-advanced-gpu-module)
- [Example 3: A Custom VTK-m Filter](#example-3-a-custom-vtk-m-filter)
- [How to Compile Vistle to Run Code on NVIDIA GPUs](#how-to-compile-vistle-to-run-code-on-nvidia-gpus)


## The VtkmModule Class

The `VtkmModule` class is the base class for VTK-m modules in Vistle. It is designed to make adding new VTK-m algorithms, so-called [filters](https://vtk-m.readthedocs.io/en/v2.2.0/provided-filters.html), as simple as possible by providing basic functionality for handling the input data, passing it to the filter implemented by the child class, and writing the filter result to the output ports. At the same time, it is designed to be flexible allowing the child class to customize and extend these processes as needed.

### Key Features of `VtkmModule`

- Reads data from input ports, validates it, and transforms it into a VTK-m dataset.
- Applies a VTK-m filter using the `runFilter` method, which must be implemented by child classes.
- Transforms the result back into Vistle objects and writes it to the output ports.

### The Constructor
A typical VTK-m module consists of one input port and one output port, but additional ports can be added by specifying the desired number in the constructor. All data on the ports must be defined on the same grid.

### The `setUpFilter` method

```cpp
VtkmModule(const std::string &name, int moduleID, mpi::communicator comm, int numPorts = 1, bool requireMappedData = true);
```

### Preparing the input data

### Preparing the output data

## Example 1: A Simple GPU Module
- `setRunFilter`
- add module paramters
- making changes to the construtor

## Example 2: A More Advanced GPU Module
- overriding prepareInput and prepareOutput methods

## Example 3: A Custom VTK-m Filter

To learn more about using VTK-m filters check out [VTK-m's user guide](https://vtk-m.readthedocs.io/en/v2.2.0/basic-filter-impl.html).

## How to Compile Vistle to Run Code on NVIDIA GPUs