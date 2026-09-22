@echo off
rem ============================================================================
rem  PDOS · 结构分析工作台 —— 一键打包成 Windows 便携程序包
rem
rem  用法：双击本文件，或在命令行执行
rem        build_portable.bat
rem
rem  产物：dist\PDOS-Viewer\PDOS-Viewer.exe   （双击即可运行，整目录可随意拷贝）
rem ============================================================================
setlocal EnableExtensions
cd /d "%~dp0"

set "PY=python"
where %PY% >nul 2>nul
if errorlevel 1 (
    echo [ERROR] 未找到 python，请先把 Python 加入 PATH。
    pause
    exit /b 1
)

echo.
echo [1/4] 检查 PyInstaller ...
%PY% -c "import PyInstaller" >nul 2>nul
if errorlevel 1 (
    echo       未安装，正在安装 pyinstaller ...
    %PY% -m pip install --upgrade pyinstaller
    if errorlevel 1 (
        echo [ERROR] 安装 PyInstaller 失败。
        pause
        exit /b 1
    )
)

echo.
echo [2/4] 清理旧产物 ...
if exist "build" rmdir /s /q "build"
if exist "dist\PDOS-Viewer" rmdir /s /q "dist\PDOS-Viewer"

echo.
echo [3/4] PyInstaller 打包（QtWebEngine 体积较大，请耐心等待）...
%PY% -m PyInstaller --noconfirm --clean PDOS-Viewer.spec
if errorlevel 1 (
    echo [ERROR] 打包失败，请查看上面的日志。
    pause
    exit /b 1
)

echo.
echo [4/4] 补充运行期需要的文件 ...
if not exist "dist\PDOS-Viewer\3dmol" mkdir "dist\PDOS-Viewer\3dmol"
copy /y "3dmol\3Dmol-min.js" "dist\PDOS-Viewer\3dmol\" >nul
if not exist "dist\PDOS-Viewer\assets" mkdir "dist\PDOS-Viewer\assets"
copy /y "assets\app.ico" "dist\PDOS-Viewer\assets\" >nul
if exist "README.md" copy /y "README.md" "dist\PDOS-Viewer\" >nul
copy /y "使用说明-便携版.txt" "dist\PDOS-Viewer\" >nul 2>nul

echo.
echo ============================================================
echo  打包完成！
echo  程序目录： %CD%\dist\PDOS-Viewer
echo  启动文件： dist\PDOS-Viewer\PDOS-Viewer.exe
echo ============================================================
echo.
pause
exit /b 0
