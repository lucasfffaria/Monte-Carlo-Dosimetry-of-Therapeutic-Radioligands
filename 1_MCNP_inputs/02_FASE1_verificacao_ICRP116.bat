@echo off
rem ============================================================
rem  PHASE 1 - verification against ICRP 116 (runs V01-V08)
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

echo PHASE 1 - verification against ICRP 116 (runs V01-V08)
call :roda V01_AM_AP
call :roda V02_AM_PA
call :roda V03_AM_LLAT
call :roda V04_AM_RLAT
call :roda V05_AF_AP
call :roda V06_AF_PA
call :roda V07_AF_LLAT
call :roda V08_AF_RLAT
echo.
echo Done: PHASE 1 - verification against ICRP 116 (runs V01-V08). See %LOG%
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
