# L-Bracket Electromagnetic Field Simulation

This project simulates and visualizes the electromagnetic field around an L-bracket antenna using the Method of Moments (MoM) and an integrated Symmetric Successive Over-Relaxation (SSOR) iterative solver.

## Geometry

The simulation analyzes an L-bracket antenna with two perpendicular arms of equal length:
- **Horizontal arm:** Extends along the positive x-axis from the origin.
- **Vertical arm:** Extends along the positive y-axis from the origin.
- **Junction:** Located at the origin (0, 0).

## Scripts Overview

### `run_all.bat`
The main script to run the entire simulation and visualization workflow. It performs the following steps:
1.  Compiles the C code (`mom.c`).
2.  Runs the compiled executable to perform the MoM/SSOR calculation.
3.  Simultaneously starts the Python scripts for visualization.
4.  Waits for the plot windows to be closed.
5.  Cleans up all generated executables, data files, and image files (`mom.exe`, `current_distribution.txt`, `geometry.txt`, `total_fields.png`, `current_distribution_plot.png`).

### `mom.c`
This C program implements the MoM and SSOR solver for the L-bracket geometry. Its primary responsibilities are:
-   **Geometry Definition:** Defines the two perpendicular strips forming the L-bracket.
-   **Matrix Assembly:** Sets up the complex impedance matrix (`Z`) and incident field vector (`V`).
-   **Solver:** Integrates an SSOR solver with an omega optimization routine to solve the system `Z * x = V` for the unknown current distribution `x`.
-   **Output:**
    -   `current_distribution.txt`: The calculated complex current on each segment of the L-bracket. The header specifies which segments belong to each arm.
    -   `geometry.txt`: The (x, y) coordinates of each segment, used for plotting.

### `magnet.py`
This Python script visualizes the total electromagnetic field around the L-bracket. It:
-   Reads the segment positions from `geometry.txt` and the current distribution from `current_distribution.txt`.
-   Calculates the total E-field (incident + scattered) on a 2D grid.
-   Overlays the L-bracket geometry on top of the field plot for context.
-   Saves the final visualization as `total_fields.png`.

### `ploter.py`
This Python script provides a detailed visualization of the current distribution on the L-bracket. It:
-   Reads `geometry.txt` to understand the structure of the antenna.
-   Reads `current_distribution.txt` for the complex current values.
-   Generates a 4-panel plot showing:
    1.  Current magnitude vs. segment index for each arm.
    2.  Current phase vs. segment index for each arm.
    3.  A 2D plot of the L-bracket with current magnitude represented by color.
    4.  Current magnitude plotted along the length of each arm.
-   Saves the final plot as `current_distribution_plot.png`.

## Input File

### `Transmitter pos`
This file configures the simulation parameters:
```
3.0           <- Strip length (applies to both arms of the L-bracket)
25.0          <- Transmitter X position
25.0          <- Transmitter Y position
```
