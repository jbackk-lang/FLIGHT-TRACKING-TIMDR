@echo off
title FLIGHT-TRACKING-TIMDR -- Przyklady
color 0B
cls

cd /d "%~dp0"

echo ============================================================
echo   FLIGHT-TRACKING-TIMDR: uruchamianie przykladow
echo   Katalog roboczy: %cd%
echo ============================================================
echo.

if exist "venv\Scripts\activate.bat" (
    echo [OK] Aktywacja venv...
    call venv\Scripts\activate.bat
) else if exist ".venv\Scripts\activate.bat" (
    echo [OK] Aktywacja .venv...
    call .venv\Scripts\activate.bat
) else (
    echo [INFO] Uzywanie systemowej instalacji Pythona.
)

echo.
echo [1/3] Instalacja pakietow pip...
python -m pip install --upgrade pip --disable-pip-version-check
python -m pip install numpy scipy matplotlib pytest

if %ERRORLEVEL% NEQ 0 (
    echo [BLAD] Instalacja pakietow nie powiodla sie.
    pause
    exit /b 1
)

echo.
echo [2/3] Testy (pytest -q)...
python -m pytest -q
if %ERRORLEVEL% NEQ 0 (
    echo [UWAGA] Niektore testy nie przeszly - przykłady mimo to zostana uruchomione,
    echo ale sprawdz powyzszy wynik testow przed zaufaniem wynikom.
)

echo.
echo [3/3] Generowanie mapy z przykladami (demo_conflict_map.py)...
echo   Przyklad 1: kurs kolizyjny, ta sama wysokosc         -^> KONFLIKT
echo   Przyklad 2: kurs kolizyjny, separacja pionowa 5000ft -^> bezpiecznie
echo   Przyklad 3: lot rownolegly, odstep 10 NM             -^> bezpiecznie
echo.
python demo_conflict_map.py

if %ERRORLEVEL% NEQ 0 (
    echo [BLAD] Generowanie mapy nie powiodlo sie.
    pause
    exit /b 1
)

echo.
echo Otwieranie conflict_alert_przyklady.png...
start "" "conflict_alert_przyklady.png"

echo.
echo (Chcesz zobaczyc tez podstawowe demo tekstowe timdr_flow/twist/predict?
echo  Uruchom osobno: python demo.py)
echo.
echo ============================================================
echo Gotowe.
echo ============================================================
pause
