# Electromagnetic Field Simulation Project

This project simulates and visualizes the electromagnetic field around a metal strip using the Method of Moments (MoM) and the Symmetric Successive Over-Relaxation (SSOR) iterative method.

## Scripts Overview

### `run_all.bat`
This is the main script to run the entire project. It orchestrates the compilation of the C code, the execution of the C programs, and the running of the Python scripts for analysis and visualization. Finally, it cleans up all generated files.

### `mom.c`
This C program implements the **Method of Moments (MoM)**. Its primary responsibilities are:
- Calculating the "exact" current distribution on a metal strip when it's subjected to an incident electromagnetic field from a transmitter.
- Setting up the linear system of equations `Z * x = V`, where `Z` is the impedance matrix, `x` is the unknown current distribution, and `V` is the incident field vector.
- It outputs `x_exact.txt` (the exact current coefficients) and `function_result.txt` (the `Z` matrix and `V` vector).

### `ssor.c`
This C program solves the complex linear system `Z * x = V` that was set up by `mom.c`. It uses the **Symmetric Successive Over-Relaxation (SSOR)** iterative method. Key features include:
- A universal 3-step omega optimization routine to find the best relaxation parameter for faster convergence.
- It takes `function_result.txt` as input and outputs `x_optimized.txt`, which contains the current distribution calculated using the optimized SSOR method.

### `magnet.py`
This Python script is responsible for the visualization of the electromagnetic field. It:
- Takes the `x_exact.txt` (exact current distribution from `mom.c`) and transmitter position as input.
- Computes and visualizes the total electric field (incident + scattered) around the metal strip.
- Saves the resulting plot as `total_fields.png`.

### `ploter.py`
This Python script compares the "exact" and "optimized" solutions for the current distribution. It:
- Reads `x_exact.txt` (from `mom.c`) and `x_optimized.txt` (from `ssor.c`).
- Plots the magnitudes of both currents for comparison.
- Calculates and prints error metrics to evaluate the accuracy of the SSOR solution.
- Saves the comparison plot as `current_comparison.png`.
