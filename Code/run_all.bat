@echo off
echo Compiling mom.c...
gcc mom.c -o mom.exe -lm
if errorlevel 1 (
    echo mom.c compilation failed.
    exit /b 1
)

echo Running MoM solver with integrated SSOR...
mom.exe

echo Starting Python visualization scripts...
start "Magnetic Field" python magnet.py
start "Current Plotter" python ploter.py

:wait
timeout /t 1 /nobreak >nul
tasklist /fi "imagename eq python.exe" | find "python.exe" >nul
if not errorlevel 1 goto wait

echo Cleaning up executable...
del mom.exe
del current_distribution.txt

echo All tasks complete.