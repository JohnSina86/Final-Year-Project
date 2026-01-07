# Electromagnetic Field Simulation Project

This project simulates and visualizes the electromagnetic field around a metal strip using the Method of Moments (MoM) and the Symmetric Successive Over-Relaxation (SSOR) iterative method. The workflow has been streamlined for efficiency and clarity.

## Scripts Overview

### `run_all.bat`
This is the main script to run the entire project. It orchestrates the compilation and execution of the C code, and the running of the Python scripts for analysis and visualization. It now performs a more focused cleanup, removing generated executables and intermediate data files. The waiting mechanism for Python scripts has been optimized for a faster experience.

### `mom.c`
This C program implements the **Method of Moments (MoM)** and now integrates the **Symmetric Successive Over-Relaxation (SSOR)** iterative solver. Its primary responsibilities are:
- Calculating the current distribution on a metal strip when it's subjected to an incident electromagnetic field from a transmitter.
- Setting up and solving the complex linear system of equations `Z * x = V`, where `Z` is the impedance matrix, `x` is the unknown current distribution, and `V` is the incident field vector.
- Includes an optimized routine to find the best relaxation parameter (`omega`) for faster SSOR convergence.
- It outputs `current_distribution.txt`, which contains the current distribution calculated using the optimized SSOR method.

### `magnet.py`
This Python script is responsible for the visualization of the electromagnetic field. It:
- Takes the `current_distribution.txt` (current distribution from `mom.c`) and transmitter position as input.
- Computes and visualizes the total electric field (incident + scattered) around the metal strip.
- Saves the resulting plot as `total_fields.png`.

### `ploter.py`
This Python script visualizes the calculated current distribution. It:
- Reads `current_distribution.txt` (from `mom.c`).
- Plots the magnitude and phase of the current distribution.
- Saves the plot as `current_distribution_plot.png`.