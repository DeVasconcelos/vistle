# How to Write a GPU Module

In this guide you will learn how to write Vistle modules that can be run on the GPU.    
Vistle makes use of the portable toolkit [VTK-m](https://m.vtk.org/) which allows running scientific visualization algorithms on various devices, including GPUs, and is designed to keep data transfers between devices at a minimum.    


## The `VtkmModule` Class

The `VtkmModule` class is the base class for VTK-m modules in Vistle. A typical VTK-m module consists of one input port and one output port, but additional ports can be added by specifying the desired number in the constructor. All data on the ports must be defined on the same grid.

### Key Features of `VtkmModule`

- Reads data from input ports, validates it, and transforms it into a VTK-m dataset.
- Applies a VTK-m filter using the `runFilter` method, which must be implemented by child classes.
- Transforms the result back into Vistle objects and writes it to the output ports.

### Constructor

```cpp
VtkmModule(const std::string &name, int moduleID, mpi::communicator comm, int numInputPorts = 1, int numOutputPorts = 1);
```

## Example 1: A Simple GPU Module
- `setRunFilter`
- add module paramters
- making changes to the construtor

## Example 2: A More Advanced GPU Module
- overriding prepareInput and prepareOutput methods

## Example 3: A Custom VTK-m Filter

To learn more about using VTK-m filters check out [VTK-m's user guide](https://vtk-m.readthedocs.io/en/v2.2.0/basic-filter-impl.html).

## How to Compile Vistle to Run Code on NVIDIA GPUs