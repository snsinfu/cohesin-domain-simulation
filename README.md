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
The simulation program (executable named `main`) takes a JSON configuration and
writes a custom HDF5 trajectory file. Run a sample simulation:

```
./main sample.json
```

The trajectory file (output.h5) can be converted into a GSD file (output.gsd)
and visualized using [OVITO](https://www.ovito.org/).

```
pixi run python ../dump_gsd.py output.h5 output.gsd
```

For statistical analysis under different conditions, see the
`tasks/P1-cohesin_nucleosome` directory.
