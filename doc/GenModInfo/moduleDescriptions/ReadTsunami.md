
# ReadTsunami
Read cheese tsunami files

This module is designed to read ChEESE tsunami files a new file format for probabilistic tsunami forecasting. Data can be generated with [Tsunami-HYSEA](https://edanya.uma.es/hysea/index.php/models/tsunami-hysea) a numerical model for quake generated tsunami simulation. The raw data output of the simulation is a netCDF file which contains the following mandatory attributes:

- **lon**: Longitude of the sea domain in WGS84 coordinates

- **lat**: Latitdude of the sea domain in WGS84 coordinates

- **grid_lat**: Latitdude of the topography in WGS84 coordinates

- **grid_lon**: Longitude of the topography in WGS84 coordinates

- **time**: Timesteps used for the simulation

- **eta**: Wave amplitude per timestep

Additional attributes like scalar data will be mapped onto the sea surface.

## How the reader works

The reader is able to read netCDF files which uses a [PnetCDF](https://parallel-netcdf.github.io/) supported file type. In general the module will fetch the longitude and latitude values for the sea (**lon**, **lat**) and the batheymetry (**grid_lon**, **grid_lat**) from the ncfile which has to be specified in the parameter browser. Based on each longitude-latitude pair the reader creates a 2D grid, varies the wave height (**eta** - orthogonal to longitude and latitude coordinates) per timestep and creates a 3D polygon surface representing the sea surface out of it. If an attribute **bathymetry** is provided with the netCDF file the reader will build a ground surface as well. 

---

## Install Requirements

- PnetCDF: build with --enable-thread-safe flag as described on [PnetCDF-GitHub](https://github.com/Parallel-NetCDF/PnetCDF)

The pnetcdf install directory should be added to your $PATH environment variable in order to be found by CMake, otherwise it won't be build.


<svg width="2700" height="360" >
<style>.text { font: normal 24.0px sans-serif;}tspan{ font: italic 24.0px sans-serif;}.moduleName{ font: italic 30px sans-serif;}</style>
<rect x="0" y="30" width="270" height="90" rx="5" ry="5" style="fill:#64c8c8ff;" />
<rect x="6.0" y="90" width="30" height="30" rx="0" ry="0" style="fill:#c8c81eff;" >
<title>Sea surface</title></rect>
<rect x="21.0" y="120" width="1.0" height="210" rx="0" ry="0" style="fill:#000000;" />
<rect x="21.0" y="330" width="30" height="1.0" rx="0" ry="0" style="fill:#000000;" />
<text x="57.0" y="333.0" class="text" >Grid Sea (Heightmap/LayerGrid)<tspan> (Sea surface)</tspan></text>
<rect x="42.0" y="90" width="30" height="30" rx="0" ry="0" style="fill:#c8c81eff;" >
<title>Ground surface</title></rect>
<rect x="57.0" y="120" width="1.0" height="180" rx="0" ry="0" style="fill:#000000;" />
<rect x="57.0" y="300" width="30" height="1.0" rx="0" ry="0" style="fill:#000000;" />
<text x="93.0" y="303.0" class="text" >Sea floor (Heightmap/LayerGrid)<tspan> (Ground surface)</tspan></text>
<rect x="78.0" y="90" width="30" height="30" rx="0" ry="0" style="fill:#c8c81eff;" >
<title>Fault wall surface</title></rect>
<rect x="93.0" y="120" width="1.0" height="150" rx="0" ry="0" style="fill:#000000;" />
<rect x="93.0" y="270" width="30" height="1.0" rx="0" ry="0" style="fill:#000000;" />
<text x="129.0" y="273.0" class="text" >Faults (Polygons)<tspan> (Fault wall surface)</tspan></text>
<rect x="114.0" y="90" width="30" height="30" rx="0" ry="0" style="fill:#c8c81eff;" >
<title>Scalar port 0</title></rect>
<rect x="129.0" y="120" width="1.0" height="120" rx="0" ry="0" style="fill:#000000;" />
<rect x="129.0" y="240" width="30" height="1.0" rx="0" ry="0" style="fill:#000000;" />
<text x="165.0" y="243.0" class="text" >Port for scalar number 0<tspan> (Scalar port 0)</tspan></text>
<rect x="150.0" y="90" width="30" height="30" rx="0" ry="0" style="fill:#c8c81eff;" >
<title>Scalar port 1</title></rect>
<rect x="165.0" y="120" width="1.0" height="90" rx="0" ry="0" style="fill:#000000;" />
<rect x="165.0" y="210" width="30" height="1.0" rx="0" ry="0" style="fill:#000000;" />
<text x="201.0" y="213.0" class="text" >Port for scalar number 1<tspan> (Scalar port 1)</tspan></text>
<rect x="186.0" y="90" width="30" height="30" rx="0" ry="0" style="fill:#c8c81eff;" >
<title>Scalar port 2</title></rect>
<rect x="201.0" y="120" width="1.0" height="60" rx="0" ry="0" style="fill:#000000;" />
<rect x="201.0" y="180" width="30" height="1.0" rx="0" ry="0" style="fill:#000000;" />
<text x="237.0" y="183.0" class="text" >Port for scalar number 2<tspan> (Scalar port 2)</tspan></text>
<rect x="222.0" y="90" width="30" height="30" rx="0" ry="0" style="fill:#c8c81eff;" >
<title>Scalar port 3</title></rect>
<rect x="237.0" y="120" width="1.0" height="30" rx="0" ry="0" style="fill:#000000;" />
<rect x="237.0" y="150" width="30" height="1.0" rx="0" ry="0" style="fill:#000000;" />
<text x="273.0" y="153.0" class="text" >Port for scalar number 3<tspan> (Scalar port 3)</tspan></text>
<text x="6.0" y="85.5" class="moduleName" >ReadTsunami</text></svg>

[outputPorts]:<>


## Parameters
|name|description|type|
|-|-|-|
|first_step|first timestep to read|Int|
|last_step|last timestep to read (-1: last)|Int|
|step_increment|number of steps to increment|Int|
|first_rank|rank for first partition of first timestep|Int|
|check_convexity|whether to check convexity of grid cells|Int|
|file_dir|NC File directory|String|
|ghost|Show ghostcells.|Int|
|fault|Show faults.|Int|
|fill|Replace filterValue.|Int|
|VerticalScale|Vertical Scale parameter sea|Float|
|blocks latitude|number of blocks in lat-direction|Int|
|blocks longitude|number of blocks in lon-direction|Int|
|fillValue|ncFile fillValue offset for eta|Float|
|fillValueNew|set new fillValue offset for eta|Float|
|bathymetry |Select bathymetry stored in netCDF)|String|
|Scalar 0|Select scalar.)|String|
|Scalar 1|Select scalar.)|String|
|Scalar 2|Select scalar.)|String|
|Scalar 3|Select scalar.)|String|

## Example Usage

![](../../../module/general/ReadTsunami/example.png)