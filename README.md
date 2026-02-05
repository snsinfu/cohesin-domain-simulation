# Chromatin domain simulations

This codebase contains the program, scripts, and notebooks for analysis of a
computational polymer model for chromatin domains.


## Running a simulation

Ue [pixi](https://prefix.dev/) for managing reproducible environment. Inside
this directory, run the following command to install dependencies:

```
pixi install
```

Build the simulation program:

```
cd src/simulation
pixi run make -j
```

You may want to edit `overrides.mk` in the source directory to correctly specify
the CPU architecture. Default is to use the one on the machine compiling the
program. Change `OPTFLAGS` macro in the file like the following to compile and
optimize for specific hardware:

```
OPTFLAGS += \
  -march=x86-64-v3 \
  -tune=znver2 \
  -flto \
  -fuse-ld=lld
```

The simulation program (executable named `main`) takes a JSON configuration file
and writes a custom HDF5 trajectory file. Run a sample simulation:

```
./main sample.json
```

The trajectory file (output.h5) can be converted into a GSD file (output.gsd)
and visualized using [OVITO](https://www.ovito.org/).

```
pixi run python ../dump_gsd.py output.h5 output.gsd
```

For statistical analysis under different conditions, see the directories under
the `tasks` directory.
