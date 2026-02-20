@echo off
SETLOCAL EnableDelayedExpansion

SET SERVER_USER=sina
SET SERVER_HOST=192.168.0.226
SET SERVER_PATH=~/gpu_benchmark

echo.
echo ======================================================================
echo   FREQUENCY SWEEP: 0.1-10 GHz (0.1 GHz steps)
echo ======================================================================
echo   Total simulations: 100
echo ======================================================================
echo.

if not exist mom_freq.c (
    echo ERROR: mom_freq.c not found!
    pause
    exit /b 1
)

gcc --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: gcc not found!
    pause
    exit /b 1
)

echo Testing server...
ssh %SERVER_USER%@%SERVER_HOST% "echo OK" >nul 2>&1
if errorlevel 1 (
    echo ERROR: Cannot connect!
    pause
    exit /b 1
)
echo   OK
echo.

echo Compiling...
gcc mom_freq.c -o mom_freq.exe -fopenmp -Ofast -march=native -lm 2>nul
if errorlevel 1 (
    gcc mom_freq.c -o mom_freq.exe -fopenmp -O3 -lm 2>nul
)
echo   OK
echo.

if not exist frequency_sweep_results mkdir frequency_sweep_results

SET /A count=0
FOR /L %%f IN (1,1,100) DO (
    SET /A count+=1
    SET /A whole=%%f / 10
    SET /A dec=%%f %% 10

    echo [!count!/100] !whole!.!dec! GHz...

    SET "freq_dir=frequency_sweep_results\freq_!whole!.!dec!GHz"
    if not exist "!freq_dir!" mkdir "!freq_dir!"

    echo !whole!.!dec! > frequency.txt

    mom_freq.exe > "!freq_dir!\mom_output.txt" 2>&1

    if exist current_distribution.txt (
        move /Y current_distribution.txt "!freq_dir!" >nul 2>&1
        move /Y geometry.txt "!freq_dir!" >nul 2>&1
        move /Y omega.txt "!freq_dir!" >nul 2>&1
        copy /Y "Transmitter pos" "!freq_dir!" >nul 2>&1

        ssh %SERVER_USER%@%SERVER_HOST% "mkdir -p %SERVER_PATH%/freq_!whole!.!dec!GHz" 2>nul
        scp -q "!freq_dir!\current_distribution.txt" %SERVER_USER%@%SERVER_HOST%:%SERVER_PATH%/freq_!whole!.!dec!GHz/ 2>nul
        scp -q "!freq_dir!\geometry.txt" %SERVER_USER%@%SERVER_HOST%:%SERVER_PATH%/freq_!whole!.!dec!GHz/ 2>nul
        scp -q "!freq_dir!\Transmitter pos" %SERVER_USER%@%SERVER_HOST%:%SERVER_PATH%/freq_!whole!.!dec!GHz/ 2>nul
        scp -q "!freq_dir!\omega.txt" %SERVER_USER%@%SERVER_HOST%:%SERVER_PATH%/freq_!whole!.!dec!GHz/ 2>nul
        scp -q frequency.txt %SERVER_USER%@%SERVER_HOST%:%SERVER_PATH%/freq_!whole!.!dec!GHz/ 2>nul

        ssh %SERVER_USER%@%SERVER_HOST% "cd %SERVER_PATH% && python3 ploter_gpu.py freq_!whole!.!dec!GHz && python3 magnet_gpu.py freq_!whole!.!dec!GHz" 2>nul

        scp -q %SERVER_USER%@%SERVER_HOST%:%SERVER_PATH%/freq_!whole!.!dec!GHz/total_fields_gpu.png "!freq_dir!" 2>nul
        scp -q %SERVER_USER%@%SERVER_HOST%:%SERVER_PATH%/freq_!whole!.!dec!GHz/current_distribution_plot.png "!freq_dir!" 2>nul
    )
)

echo.
echo Generating summary...
ssh %SERVER_USER%@%SERVER_HOST% "cd %SERVER_PATH% && python3 analyze_frequency_sweep.py" 2>nul
scp -q %SERVER_USER%@%SERVER_HOST%:%SERVER_PATH%/omega_vs_frequency.png frequency_sweep_results/ 2>nul
scp -q %SERVER_USER%@%SERVER_HOST%:%SERVER_PATH%/frequency_sweep_summary.csv frequency_sweep_results/ 2>nul

echo.
echo ======================================================================
echo   DONE! Check frequency_sweep_results\
echo ======================================================================
echo.

del frequency.txt 2>nul
del mom_freq.exe 2>nul

pause
