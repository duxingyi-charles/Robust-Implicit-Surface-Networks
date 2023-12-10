# Robust Computation of Implicit Surface Networks for Piecewise Linear Functions

![](https://user-images.githubusercontent.com/3606672/187550647-900c8e12-bd57-4270-97da-15a0180e0c39.png)

[Xingyi Du](https://duxingyi-charles.github.io/), [Qingnan Zhou](https://research.adobe.com/person/qingnan-zhou/),  [Nathan Carr](https://research.adobe.com/person/nathan-carr/), [Tao Ju](https://www.cse.wustl.edu/~taoju/)
<br/>*ACM Transaction on Graphics (Proceedings of SIGGRAPH 2022)*<br/>

[`Project Page`](https://duxingyi-charles.github.io/publication/robust-computation-of-implicit-surface-networks-for-piecewise-linear-functions/)

## Build

Use the following commands to build on Mac.

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

Program `impl_arrangement` and `material_interface` will be generated in the `build` subdirectory.

## Usage

### `impl_arrangement`

The program creates a polygonal mesh discretization of the implicit arrangement of a collection of implicit functions on a given tetrahedron grid.

Usage: 

    ./impl_arrangement [OPTIONS] config_file

Options:
- -h,--help: Print help message and exit.
- -T,--timing-only: When set, the program only records timing without saving results.
- -R,--robust-test: When set, the program performs robustness test. The program will run twice with different ordering of input functions, and check if the results are consistent.

Positionals:
- `config_file` (REQUIRED): Configuration file, specifying input/output paths and algorithm parameters.

The `config_file` should be a JSON file with the following named parameters:
- `tetMeshFile`: Path to the file storing input tetrahedral grid. It is a JSON file storing vertex coordinates and tetrahedrons. See `examples/tet_mesh/tet5_grid_10k.json` for a concrete example.
- `funcFile`: Path to the file storing input implicit functions. See this [repo](https://github.com/duxingyi-charles/implicit_functions) for details about the file format.
- `outputDir`: Path to the directory to store output files.
- `useLookup`: Whether to use look-up tables to accelerate tetrahedron processing (section 6 of our paper). Default value is "true".
- `useSecondaryLookup`: Whether to use look-up tables for tetrahedral with two active functions (section 6 of our paper). Default value is "true".
- `useTopoRayShooting`: Whether to use topological ray shooting to compute the spatial decomposition induced by the arrangement (section 7 of our paper). Default value is "true".

Note that `tetMeshFile`, `funcFile` and `outputDir` can be either absolute
or relative paths.  In the case of relative paths, they are relative with
respect to the directory containing the configuration file.

An example config file is `examples/implicit_arrangement/config.json`.

Test:

compute implicit arrangement

    ./impl_arrangement ../examples/implicit_arrangement/config.json

record timing without saving results

    ./impl_arrangement -T ../examples/implicit_arrangement/config.json

perform robustness test

    ./impl_arrangement -R ../examples/implicit_arrangement/config.json

### `material_interface`

The program creates a polygonal mesh discretization of the material interfaces of a collection of material functions on a given tetrahedron grid.

Usage:

    ./material_interface [OPTIONS] config_file

Options:
- -h,--help: Print help message and exit.
- -T,--timing-only: When set, the program only records timing without saving results.
- -R,--robust-test: When set, the program performs robustness test. The program will run twice with different ordering of input functions, and check if the results are consistent.

Positionals:
- `config_file` (REQUIRED): Configuration file, specifying input/output paths and algorithm parameters.

The `config_file` should be a JSON file with the following named parameters:
- `tetMeshFile`: Path to the file storing input tetrahedral grid. It is a JSON file storing vertex coordinates and tetrahedrons. See `examples/tet_mesh/tet5_grid_10k.json` for a concrete example.
- `funcFile`: Path to the file storing input material functions. See this [repo](https://github.com/duxingyi-charles/implicit_functions) for details about the file format.
- `outputDir`: Path to the directory to store output files.
- `useLookup`: Whether to use look-up tables to accelerate tetrahedron processing (section 6 of our paper). Default value is "true".
- `useSecondaryLookup`: Whether to use look-up tables for tetrahedral with three active functions (section 6 of our paper). Default value is "true".
- `useTopoRayShooting`: Whether to use topological ray shooting to compute the spatial decomposition induced by the material interfaces (section 7 of our paper). Default value is "true".

Note that `tetMeshFile`, `funcFile` and `outputDir` can be either absolute
or relative paths.  In the case of relative paths, they are relative with
respect to the directory containing the configuration file.

An example config file is `examples/material_interface/config.json`.

Test:

compute material interface

    ./material_interface ../examples/material_interface/config.json

record timing without saving results

    ./material_interface -T ../examples/material_interface/config.json

perform robustness test

    ./material_interface -R ../examples/material_interface/config.json

## Output

The complete set of output files include data files (`mesh.json`, `mesh_patches.msh`, `mesh_chains.msh` and `mesh_cells.msh`)
and information files (`timings.json` and `stats.json`). In timing-only mode (`-T`), the program only generates information files.
In robustness test mode (`-R`), no output files are generated. The result of the robustness test is in command line output.

### data files

`mesh.json`:
JSON file storing the following key-value pairs

| key     | value                         | description                                                                                                                                                                                                                                              |
|---------|-------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| points  | #Vx3 matrix of double numbers | Vertex coordinates of the surface network mesh.                                                                                                                                                                                                          |
| faces   | vector of vector of integers  | Polygonal faces of the surface network mesh. Each face is encoded by a vector of indices of face boundary vertices.                                                                                                                                      |
| edges   | #Ex2 matrix of integers       | Edges of the surface network mesh. Each edge is encoded by a pair of vertex indices.                                                                                                                                                                     |
| chains  | vector of vector of integers  | Chains of non-manifold edges. Each chain is encoded by a vector of edge indices.                                                                                                                                                                         |
| corners | vector of integers            | Indices of non-manifold vertices.                                                                                                                                                                                                                        |
| patches | vector of vector of integers  | A patch is a connected component of faces bounded by non-manifold edges.  Each patch is encoded by a vector of face indices.                                                                                                                             |
| shells  | vector of vector of integers  | A shell is a connected component of the boundary of a 3D region partitioned by the surface network. Each shell is encoded by a vector of oriented patch indices. The positive side of patch i has index 2i. The negative side of patch i has index 2i+1. |
| cells   | vector of vector of integers  | A cell is a 3D region partitioned by the surface network. Each cell is encoded by a vector of shell indices.                                                                                                                                             |

Note: all indices start from 0.

`mesh_patches.msh`, `mesh_chains.msh` and `mesh_cells.msh` store the same data as chains, patches and cells in `mesh.json`, but in MSH format.
We can view these MSH files using [Gmsh](https://gmsh.info/).

### information files

`timing.json`: timing of different stages of our pipeline.

`stats.json`: statistics of intermediate data in our pipeline, e.g., number of surface network mesh vertices.
