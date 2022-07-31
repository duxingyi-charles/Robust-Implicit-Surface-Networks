# Robust Computation of Implicit Surface Networks for Piecewise Linear Functions

![](https://user-images.githubusercontent.com/8947527/182003910-f1e39f30-e67b-498b-b990-c2f5964d5020.jpg)

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

### impl_arrangement

The program creates a polygonal mesh discretization of the implicit arrangement of a collection of implicit functions on a given tetrahedron grid.

Usage: 

    ./impl_arrangement [OPTIONS] config_file

Positionals:
- config_file (REQUIRED): Configuration file, specifying input/output paths and algorithm parameters. See below for details.

Options:
- -h,--help: Print help message and exit.
- -T,--timing-only BOOLEAN: When set to 1 (True), the program only records timing without saving results. Default value is 0 (False).
- -R,--robust-test BOOLEAN: When set to 1 (True), the program performs robustness test. The program will run twice with different ordering of input functions, and check if the results are consistent. Default value is 0 (False).

The `config_file` should be a JSON file with the following named parameters:
- `tetMeshFile`: Absolute path to the file storing input tetrahedral grid. It is a JSON file storing vertex coordinates and tetrahedrons. See `example/tet_mesh/tet5_grid_10k.json` for a concrete example.
- `funcFile`: Absolute path to the file storing input implicit functions. See this [repo](https://github.com/duxingyi-charles/implicit_functions) for details about the file format.
- `outputDir`: Absolute path to the directory to store output files.
- `useLookup`: Whether to use look-up tables to accelerate tet processing (see section 7 of our paper). Default value is "true".
- `use2funcLookup`: Whether to use look-up tables for tetrahedral with two active functions (see section 6 of our paper). Default value is "true".
- `useTopoRayShooting`: Whether to use topological ray shooting to resolve the nesting order of multiple connected components of arrangement mesh (see section 7 of our paper). Default value is "true".

Note that `tetMeshFile`, `funcFile` and `outputDir` take **absolute** paths.

An example config file is `example\implicit_arrangement\config.json`. You should change the paths in the config file according to your own system. 

Test:

compute implicit arrangement

    ./impl_arrangement ../example/implicit_arrangement/config.json

record timing without saving results

    ./impl_arrangement -T 1 ../example/implicit_arrangement/config.json

perform robustness test

    ./impl_arrangement -R 1 ../example/implicit_arrangement/config.json

### material_interface

Usage: 

    ./material_interface [OPTIONS] config_file

Positionals:
- config_file (REQUIRED): Configuration file, specifying input/output paths and algorithm parameters.
Options:
- -h,--help : Print help message and exit.
- -T,--timing-only BOOLEAN    Record timing without output result
- -R,--robust-test BOOLEAN    Perform robustness test.
