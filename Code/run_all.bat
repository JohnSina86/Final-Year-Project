@echo off

:: Compile C source files
echo Compiling C source files...
gcc mom.c -o mom.exe -lm
if errorlevel 1 (
    echo mom.c compilation failed.
    exit /b 1
)
gcc ssor.c -o ssor.exe -lm
if errorlevel 1 (
    echo ssor.c compilation failed.
    del mom.exe
    exit /b 1
)

:: Run C executables
echo Running C executables...
mom.exe
ssor.exe

:: Run Python scripts simultaneously
echo Starting Python scripts...
echo Close the plot windows to continue the script.
start "magnet" python magnet.py
start "ploter" python ploter.py

:wait
echo Waiting for python scripts to close...
timeout /t 5 /nobreak >nul
tasklist /fi "imagename eq python.exe" | find "python.exe" >nul
if not errorlevel 1 goto wait

:: Cleanup
echo Cleaning up generated files...
del mom.exe
del ssor.exe
del function_result.txt
del x_exact.txt
del x_estimate.txt
del x_optimized.txt

echo All tasks complete.