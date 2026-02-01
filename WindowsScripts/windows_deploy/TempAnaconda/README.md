# TempAnaconda (Windows) – stabil neben Anaconda

Problem: Anaconda ist installiert, aber die Notebook-Umgebung ist “vermischt” (conda/python/jupyter/pip im PATH).
Das führt zu kaputten Imports oder Installation in der falschen Umgebung.

Lösung: Wir erstellen eine **separate venv** im Projekt und registrieren einen **eigenen Kernel**:
`WTConc2 (venv)`. Damit ist klar, welche Pakete und welcher Code verwendet wird.

## Was du bekommst
- **Separates Python** (venv) neben Anaconda (keine conda-Abhängigkeit)
- Jupyter startet **aus der venv** (empfohlen)
- Optional: Du kannst Anaconda-Jupyter weiter nutzen, aber Kernel **WTConc2 (venv)** auswählen

## Voraussetzungen
- Windows + Python von `python.org` (empfohlen: **3.11.x**) inkl. **py launcher**

## Installation
Im Repo unter `WindowsScripts\\windows_deploy\\TempAnaconda\\`:

1) Setup (online oder offline, siehe unten):
- `setup_WTConc2_venv.bat`

2) Jupyter starten:
- `start_WTConc2_Jupyter.bat`

## Offline-Variante (ohne Internet auf dem Zielrechner)
Auf einem **Online-Rechner**:
- `setup_WTConc2_venv.bat download-wheels`
  - erzeugt `TempAnaconda\\wheels\\`

Dann `wheels\\` auf den Offline-Rechner kopieren (gleicher Pfad) und dort:
- `setup_WTConc2_venv.bat`

## Falls Anaconda-Jupyter genutzt werden muss
- Anaconda Jupyter starten
- Notebook → **Kernel** → **Change Kernel** → **WTConc2 (venv)**

## Typische Fehler
- **“py launcher not found”**: Python von python.org installieren (mit py launcher)
- **Kernel fehlt**: `setup_WTConc2_venv.bat` erneut ausführen (registriert den Kernel)


