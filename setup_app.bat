@echo off
:: ==========================================
:: Script: setup_conda_app.bat
:: Purpose: Create Conda environment and install dependencies for Codon AI app
:: ==========================================

:: Check if conda exists
where conda >nul 2>&1
IF %ERRORLEVEL% NEQ 0 (
    echo ⚠️ Conda is not installed. Please install Anaconda or Miniconda first.
    pause
    exit /b 1
)

echo ✅ Conda detected.

:: Environment name
set ENV_NAME=codon_ai_env

:: Create environment
echo Creating Conda environment "%ENV_NAME%" with Python 3.10...
conda create -y -n %ENV_NAME% python=3.10

IF %ERRORLEVEL% NEQ 0 (
    echo ❌ Failed to create Conda environment.
    pause
    exit /b 1
)

echo Activating environment...
call conda activate %ENV_NAME%

:: Install required packages
echo Installing required packages...
conda install -y -c conda-forge streamlit biopython matplotlib

:: Optional: If you have extra pip packages
IF EXIST requirements.txt (
    pip install -r requirements.txt
)

echo ==================================
echo ✅ Setup Complete!
echo To run the app:
echo call conda activate %ENV_NAME%
echo cd streamlit_app
echo streamlit run app.py
pause