@echo off
rem ============================================================
rem  SUPPLEMENTARY TESTS - article 1 (runs R1 and R2, 4 November 2026)
rem    R1_P0_Tb161_direto  : 161Tb photon lines versus the same lines folded on the energy grid
rem    R2_P0_piso_concreto : patient-alone run P0 with a 30 cm concrete floor under the patient
rem  CONFIGURATION: set the MCNP executable and the number of threads
rem ============================================================
set MCNP=mcnp6
set TASKS=20
rem 1 = delete the RUNTPE file after each successful run (the tests are not continued)
set APAGA_RUNTPE=1
rem ============================================================
rem  Folders (relative to this file):
rem    MCTAL (official results)  -> ..\2_MCNP_results
rem    full output files (.o)    -> ..\..\Dados\saidas_MCNP_completas
rem    RUNTPE restart files (.r) -> ..\..\Dados\continuacao\trabalho
rem  A run is skipped when its MCTAL file already exists, so the official results are never overwritten.
rem  CTME 1200 (computer time summed over the threads): about 1 h per run on 20 threads.
rem ============================================================
cd /d "%~dp0"
set RES=..\2_MCNP_results
set SAI=..\..\Dados\saidas_MCNP_completas
set TRAB=..\..\Dados\continuacao\trabalho
set LOG=%SAI%\log_testes_R1_R2.txt
if not exist %RES% mkdir %RES%
if not exist %SAI% mkdir %SAI%
if not exist %TRAB% mkdir %TRAB%

echo SUPPLEMENTARY TESTS - article 1 (runs R1 and R2)
call :roda R1_P0_Tb161_direto
call :roda R2_P0_piso_concreto
echo.
echo Done: SUPPLEMENTARY TESTS - article 1 (runs R1 and R2). See %LOG%
echo Analysis: python a1_05_supplementary_tests.py in ..\4_analysis
pause
goto :eof

:roda
set N=%~1
if exist "%RES%\%N%.m" (
  echo [skipped] %N% already done ^(%RES%\%N%.m exists^)
  goto :eof
)
echo.
echo ==================== %N% ====================
echo [%date% %time%] INICIO %N% >> %LOG%
if exist "%SAI%\%N%.o" del /q "%SAI%\%N%.o"
if exist "%TRAB%\%N%.r" del /q "%TRAB%\%N%.r"
%MCNP% inp=%N%.i outp=%SAI%\%N%.o runtpe=%TRAB%\%N%.r mctal=%RES%\%N%.m tasks %TASKS%
if exist "%RES%\%N%.m" (
  echo [%date% %time%] FIM OK %N% >> %LOG%
  if "%APAGA_RUNTPE%"=="1" del /q "%TRAB%\%N%.r"
) else (
  echo [%date% %time%] FALHOU %N% - check %SAI%\%N%.o >> %LOG%
)
goto :eof
