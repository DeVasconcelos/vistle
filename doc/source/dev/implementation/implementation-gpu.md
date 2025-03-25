# How to Write a GPU Module

This guide explains how to write Vistle modules that can be run on the GPU.    
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

The base constructor, which the child class should call in its own constructor, creates the input and output ports of the VTK-m module.

```cpp
VtkmModule(const std::string &name, int moduleID, mpi::communicator comm, int numPorts = 1, bool requireMappedData = true);
```

By default, a VTK-m module consists of one input port and one output port, but additional ports can be added by specifying the desired number in the base constructor (`numPorts`). Note that all data on the input ports must be defined on the same grid, as the module will throw an error, otherwise.

Many VTK-m filters work on data fields, so, by default, a VTK-m module expects the data at the input ports to contain mapped data in addition to a grid and will throw an error if there is none. If desired, the child class can remove this requirement by setting `requireMappedData` to `false`.

### Defining the VTK-m Filter

The `setUpFilter` method, which must be implemented by the child class, creates, sets up and returns the desired VTK-m filter that will be called on the input data.  

```cpp
virtual std::unique_ptr<vtkm::filter::Filter> setUpFilter() const = 0;
```
  
If the VTK-m module has multiple input ports, the filter will only be applied to the data on the first input port, i.e., the filter's active field is set to the field on the first port. The fields on the remaining ports will be mapped to the resulting output grid.

### Preparing the Input Data

The `prepareInputGrid` transforms the input grid into a VTK-m cellset and adds it to the VTK-m dataset `dataset`. Similarly, `prepareOutputGrid`, which is called for each field, transforms the input fields into VTK-m array handles and adds them to `dataset` as well. The filter will, subsequently, be applied to `dataset`.
```cpp
virtual ModuleStatusPtr prepareInputGrid(const vistle::Object::const_ptr &grid, vtkm::cont::DataSet &dataset) const;
virtual ModuleStatusPtr prepareInputField(const vistle::Port *port, const vistle::Object::const_ptr &grid, const vistle::DataBase::const_ptr &field, std::string &fieldName, vtkm::cont::DataSet &dataset) const;
```

A VTK-m module only performs very basic checks on the input ports while reading in the data, i.e., in the `readInPorts` method: It ensures each input port contains data as long as its corresponding output port is connected. Additionally, it makes sure all data fields are defined on the same grid and that at least one input grid provides an input grid.   
Some filters might, however, require additional checks. These can be added by overridding `prepareInputGrid` and/or `prepareInputField` as needed.

### Preparing the Output Data

## Example 1: A Basic GPU Module
- `setRunFilter`
- add module paramters
- making changes to the construtor

## Example 2: Extending the Core Functionality
- overriding prepareInput and prepareOutput methods

## Example 3: A Custom VTK-m Filter

To learn more about using VTK-m filters check out [VTK-m's user guide](https://vtk-m.readthedocs.io/en/v2.2.0/basic-filter-impl.html).

## How to Compile Vistle to Run Code on NVIDIA GPUs