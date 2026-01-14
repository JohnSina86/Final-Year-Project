@echo off
echo.
echo ============================================================
echo    L-Bracket MoM Analysis - Complete Pipeline
echo ============================================================
echo.

REM ============================================================
REM Step 1: Compile with OpenMP and optimizations
REM ============================================================
echo [1/3] Compiling mom.c with optimizations...
gcc mom.c -o mom.exe -fopenmp -Ofast -march=native -lm

if errorlevel 1 (
    echo.
    echo *** Compilation FAILED! ***
    echo Trying fallback compilation...
    gcc mom.c -o mom.exe -fopenmp -O3 -lm

    if errorlevel 1 (
        echo *** Fallback also failed! Check your code. ***
        pause
        exit /b 1
    )
)

echo    Compilation successful!
echo.

REM ============================================================
REM Step 2: Run MoM solver
REM ============================================================
echo ============================================================
echo [2/3] Running MoM solver...
echo ============================================================
echo.

mom.exe

if errorlevel 1 (
    echo.
    echo *** Solver FAILED! ***
    pause
    exit /b 1
)

echo.
echo    Solver completed successfully!
echo.

REM ============================================================
REM Step 3: Launch visualization scripts
REM ============================================================
echo ============================================================
echo [3/3] Starting Python visualization...
echo ============================================================
echo.

REM Check if output files exist
if not exist "current_distribution.txt" (
    echo *** WARNING: current_distribution.txt not found! ***
)

if not exist "geometry.txt" (
    echo *** WARNING: geometry.txt not found! ***
)

echo Launching visualization scripts...
echo   - Magnetic field plot (magnet.py)
echo   - Current distribution plot (ploter.py)
echo.

start "Magnetic Field" python magnet.py
start "Current Plotter" python ploter.py

echo Visualization windows opened.
echo Close the plot windows when done viewing.
echo.

REM Wait for Python scripts to finish
:wait
timeout /t 1 /nobreak >nul
tasklist /fi "imagename eq python.exe" 2>nul | find "python.exe" >nul
if not errorlevel 1 goto wait

echo.
echo ============================================================
echo    Pipeline Complete!
echo ============================================================
echo.
echo Generated files:
if exist "current_distribution.txt" echo   [OK] current_distribution.txt
if exist "geometry.txt" echo   [OK] geometry.txt
if exist "total_fields.png" echo   [OK] total_fields.png
if exist "current_distribution_plot.png" echo   [OK] current_distribution_plot.png
echo.

REM ============================================================
REM Optional cleanup - only delete executable
REM ============================================================
echo Cleaning up...
if exist "mom.exe" (
    del mom.exe
    echo   [Deleted] mom.exe
)

echo.
echo ============================================================
echo All tasks complete! Results saved.
echo ============================================================
pause
