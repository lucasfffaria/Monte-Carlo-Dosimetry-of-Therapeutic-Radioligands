@echo off
rem ============================================================
rem  PHASE 2 - article 1, standing patient (runs P0 and T01-T08)
rem  Lu-177 PSMA project - MCNP6.3 production runs (article 1)
rem  CONFIGURATION: set the MCNP executable and the number of threads
rem ============================================================
set MCNP=mcnp6
set TASKS=20
rem 1 = delete the RUNTPE file after each successful run (saves disk space)
set APAGA_RUNTPE=0
rem ============================================================
rem  Folders (relative to this file):
rem    MCTAL (official results)  -> ..\2_MCNP_results
rem    full output files (.o)    -> ..\..\Dados\saidas_MCNP_completas
rem    RUNTPE restart files (.r) -> ..\..\Dados\continuacao\trabalho
rem  A run is skipped when its MCTAL file already exists, so the official results are never overwritten.
rem ============================================================
cd /d "%~dp0"
set RES=..\2_MCNP_results
set SAI=..\..\Dados\saidas_MCNP_completas
set TRAB=..\..\Dados\continuacao\trabalho
set LOG=%SAI%\log_execucao.txt
if not exist %RES% mkdir %RES%
if not exist %SAI% mkdir %SAI%
if not exist %TRAB% mkdir %TRAB%

echo PHASE 2 - article 1, standing patient (runs P0 and T01-T08)
call :roda P0_paciente_em_pe
call :roda T01_AM_g030
call :roda T02_AM_g050
call :roda T03_AM_g100
call :roda T04_AM_g200
call :roda T05_AF_g030
call :roda T06_AF_g050
call :roda T07_AF_g100
call :roda T08_AF_g200
echo.
echo Done: PHASE 2 - article 1, standing patient (runs P0 and T01-T08). See %LOG%
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
