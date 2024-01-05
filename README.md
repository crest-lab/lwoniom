
#  lwONIOM

This repository contains a standalone implementation of a *light-weight* ONIOM interface.

The default CMake and meson builds compile a statically linked library (`liblwoniom.a`) that can be linked in other projects.

The library contains implementations for partitioning of a molecular system into the individual ONIOM fragments, including saturation of bonds with hydrogen linking atoms, and the reconstruction of respective energies, gradients and Hessian matrices.

Adhereing to the eponymous *light-weight* aspect, no potentials are implemented in the library itself as it is intended for the use in other codes.
A standandalone command line app can be built with the `-Dbuild_exe=true` option, which offers a demonstration of input file formats and capabilities to test the ONIOM partition.
`main.f90` in `app/` demonstrates the library's in-code usage.

## Building the Project

Make sure you have the following dependencies installed:

- CMake and `make`, or meson and ninja build systems
- Fortran and C compilers (e.g., `gfortran`/`gcc` or `ifort`/`icc`)

### Instructions

Follow these steps to build the project:

1. Create a build directory and navigate to it
   ```bash
   mkdir _build
   cd _build
   ```

2. Export the compilers (here for example `ifort`/`icc`) and depending on your chosen build system set up the build:
   - generate the build files using CMake:
     ```bash
     FC=ifort CC=icc cmake ..
     ```
   - generate the build files using meson:
     ```bash
     FC=ifort CC=icc meson ..
     ```
   **If you wish to build the an app binary, add `-Dbuild_exe=true` to either the `cmake` or `meson` setup command.**


3. Depending on your chosen build system, build the project. If you have multiple cores/processors, you can speed up the build process by specifying the number of cores to use with the `-j` option. For example, to use 4 cores:
   - With CMake/`make`:
     ```shell
     make -j4
     ```
   - With meson/`ninja`:
     ```shell
     ninja -j4
     ```
### Cleaning the Build

To clean the build files, simply delete the `build` directory:

```shell
rm -rf _build
```

