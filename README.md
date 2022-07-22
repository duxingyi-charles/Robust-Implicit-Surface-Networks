# Robust Computation of Implicit Surface Networks for Piecewise Linear Functions

![](figure/featured.jpg)

[Xingyi Du](https://duxingyi-charles.github.io/), [Qingnan Zhou](https://research.adobe.com/person/qingnan-zhou/),  [Nathan Carr](https://research.adobe.com/person/nathan-carr/), [Tao Ju](https://www.cse.wustl.edu/~taoju/)
<br/>*ACM Transaction on Graphics (Proceedings of SIGGRAPH 2022)*<br/>

[`Project Page`](https://duxingyi-charles.github.io/publication/robust-computation-of-implicit-surface-networks-for-piecewise-linear-functions/)

## Abstract

Implicit surface networks, such as arrangements of implicit surfaces and materials interfaces, are used for modeling piecewise smooth or partitioned shapes. However, accurate and numerically robust algorithms for discretizing either structure on a grid are still lacking. We present a unified approach for computing both types of surface networks for piecewise linear functions defined on a tetrahedral grid. Both algorithms are guaranteed to produce a correct combinatorial structure for any number of functions. Our main contribution is an exact and efficient method for partitioning a tetrahedron using the level sets of linear functions defined by barycentric interpolation. To further improve performance, we designed look-up tables to speed up processing of tetrahedra involving few functions and introduced an efficient algorithm for identifying nested 3D regions.

## Code

- impl_arrangement, a command line program that discretizes the arrangements of the 0-level sets of a collection of functions on a tetrahedron grid. 
- material_interface, a command line program that discretizes the material interfaces of a collection of material functions on a tetrahedron grid.

The programs have been tested on macOS 12.4.

## Build

### Mac

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

The program `impl_arrangement` and `material_interface` will be generated in the `build` subdirectory.

## Usage

### implicit_arrangement

Usage: 

    ./impl_arrangement [OPTIONS] config_file

Positionals:
- config_file (REQUIRED): Configuration file, specifying input halfspaces and samples. See `/examples/xxx/input/config.json` for examples.

Options:
- -h,--help : Print help message and exit.
- -T,--timing-only BOOLEAN    Record timing without output result
- -R,--robust-test BOOLEAN    Perform robustness test.

Example:

    ./BSH_CLI ../examples/figure16/tori/input/grid_64.json -P ../examples/figure16/tori/input/param.json -A mesh ../examples/figure16/tori/input/config.json  ../examples/figure16/tori/output/result.grid

### material_interface

Usage: 

    ./material_interface [OPTIONS] config_file

Positionals:
- config_file (REQUIRED): Configuration file, specifying input halfspaces and samples. See `/examples/xxx/input/config.json` for examples.

Options:
- -h,--help : Print help message and exit.
- -T,--timing-only BOOLEAN    Record timing without output result
- -R,--robust-test BOOLEAN    Perform robustness test.
