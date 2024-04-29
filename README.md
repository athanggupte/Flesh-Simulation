# Stable Neo-Hookean Flesh Simulation using Finite Element Method

## Build

Generate project files using CMake:

    mkdir build
    cmake -S . -B build

Then build the project:

    cmake --build build

## Run

Run the simulation:

	./build/driver [options...]

    Options:
    -n <frames>         Number of frames to write to file
    -w                  Write simulation to file
    -p <file>           Playback simulation from file
    -f <file>           Model file to load (default: cube_bar)
    -Y <youngs modulus> Young's modulus (default: 1e6)
    -P <poisson ratio>  Poisson ratio (default: 0.49)
    -m <material type>  Material type [stvk, stneo] (default: stneo)

Examples:
  - To run a new simulation with default parameters:
	```
	./build/driver
	```
  - To run a new simulation and write the simulation to file:
    ```
    ./build/driver -w -n 100
    ```
    Simulation is saved to `sim.dat`
  - To playback a simulation from file:
	```
	./build/driver -p sim.dat
	```
