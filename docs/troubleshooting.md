# Troubleshooting Guide

This document provides solutions for common issues encountered when running or building the HPA Cancer/IHC Image Downloader.

---

## 1️⃣ Application Does Not Start (Windows Executable)

### Symptom
Double-clicking `ImageDownloaderHPA.exe` does nothing or closes immediately.

### Possible Causes
- Windows Defender blocking the executable
- Antivirus quarantine
- Missing dependency during build

### Solutions
- Right-click the file → Properties → If "Unblock" appears, enable it.
- Add the project folder to Windows Defender exclusions.
- Rebuild the executable using the debug version (without `--windowed`).

---

## 2️⃣ ModuleNotFoundError: bs4

### Symptom
Error appears when running the Python version or the executable.

### Cause
BeautifulSoup is not installed or not included in the build.

### Fix (Development Mode)
```bash
pip install beautifulsoup4