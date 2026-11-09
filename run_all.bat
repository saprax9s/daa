@echo off
setlocal enabledelayedexpansion
echo ==========================================================
echo Starting DAA End-to-End Build and Benchmark Pipeline
echo ==========================================================

REM Detect python interpreter with pandas/matplotlib
set "PYTHON_EXE="
where py >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    py -c "import matplotlib, pandas" >nul 2>&1
    if !ERRORLEVEL% EQU 0 set "PYTHON_EXE=py"
)

if "%PYTHON_EXE%"=="" (
    if exist "%LOCALAPPDATA%\Programs\Python\Python313\python.exe" (
        "%LOCALAPPDATA%\Programs\Python\Python313\python.exe" -c "import matplotlib, pandas" >nul 2>&1
        if !ERRORLEVEL! EQU 0 (
            set "PYTHON_EXE=%LOCALAPPDATA%\Programs\Python\Python313\python.exe"
        )
    )
)

if "%PYTHON_EXE%"=="" (
    python -c "import matplotlib, pandas" >nul 2>&1
    if !ERRORLEVEL! EQU 0 (
        set "PYTHON_EXE=python"
    ) else (
        set "PYTHON_EXE=python"
        echo [WARNING] Python with matplotlib/pandas not detected in PATH.
    )
)

echo [INFO] Using Python: %PYTHON_EXE%
echo.

cd /d "%~dp0\codon\cpp"
call build.bat
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] Codon C++ build failed.
    exit /b %ERRORLEVEL%
)
.\build\codon_tool.exe --selftest
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] Codon selftest failed.
    exit /b %ERRORLEVEL%
)
.\build\codon_tool.exe --bench
cd ..\python
%PYTHON_EXE% plots.py
%PYTHON_EXE% make_tables.py

echo.
echo ==========================================================
echo Codon Pipeline Complete. Starting Merkle Pipeline...
echo ==========================================================
cd /d "%~dp0\merkle\cpp"
call build.bat
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] Merkle C++ build failed.
    exit /b %ERRORLEVEL%
)
.\build\merkle_tool.exe --selftest
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] Merkle selftest failed.
    exit /b %ERRORLEVEL%
)
.\build\merkle_tool.exe --bench
cd ..\python
%PYTHON_EXE% plots.py
%PYTHON_EXE% make_tables.py

echo.
echo ==========================================================
echo ALL PIPELINES COMPLETED SUCCESSFULLY!
echo ==========================================================
cd /d "%~dp0"
