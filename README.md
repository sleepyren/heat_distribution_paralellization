Heat Distribution Simulation
Overview

This program simulates the heat distribution across a 2D grid using both a sequential CPU version and a parallelized version utilizing OpenMP. The simulation calculates the temperature distribution within a square space where the edges have fixed temperatures and the interior point temperatures depend on their neighboring points.
Compilation and Execution

To compile and run this program, ensure that you have GCC with OpenMP support on your system. The following steps will guide you through compiling and executing the program:

Load the necessary compiler:


module load gcc-12.2

Compile the code:

gcc -Wall -std=c99 -fopenmp -o heatdist heatdist.c -lm

Execute the program:


./heatdist <dimension> <iterations> <version> <threads>

        dimension: Dimension of the 2D square matrix (NxN).
        iterations: Number of iterations the simulation runs.
        version: 0 for sequential CPU version, 1 for OpenMP parallel version.
        threads: Number of threads to use for the OpenMP version.

Example:


./heatdist 100 1000 1 4

Code Structure
Main Components

    main(): Orchestrates the initialization, execution of the heat distribution function (sequential or parallel), and timing of the simulation.
    seq_heat_dist(float *, unsigned int, unsigned int): Sequential version of the heat distribution calculation.
    parallel_heat_dist(float *, unsigned int, unsigned int): Parallel version utilizing OpenMP for heat distribution.
    check_result(int, unsigned int, float *): Validates that the parallel version produces the same results as the sequential version.

Parallel Heat Distribution

The parallel_heat_dist function is the core of the OpenMP implementation. It uses a double-buffering technique to manage dependencies between iterations effectively. Two arrays are used: one holds the current temperatures, and the other stores updated temperatures. At the end of each iteration, the arrays are swapped.
Utilities

    index(i, j, N): Macro to convert 2D array indices into a 1D array index, facilitating easier data management in a linear memory space.

Notes

    Ensure that the dimensions and the number of iterations are chosen considering the capabilities of the executing system, especially for large sizes, as the computational and memory requirements increase significantly.
    The OpenMP version will likely perform better on systems with multiple cores, particularly for larger grid sizes and higher iteration counts.

