# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller 打包配置 —— PDOS · 结构分析工作台（Windows 便携版）

构建：
    pyinstaller --noconfirm --clean PDOS-Viewer.spec

产物：
    dist/PDOS-Viewer/PDOS-Viewer.exe    ← 双击即可运行
    dist/PDOS-Viewer/_internal/         ← 运行时依赖（QtWebEngine / Qt / matplotlib）

打包后程序的目录约定（见 pdos_viewer.py 顶部 APP_DIR / RES_DIR）：
    APP_DIR = exe 所在目录      —— 数据文件夹、导出文件、运行期生成的 _viewer.html
    RES_DIR = sys._MEIPASS      —— 打包进来的只读资源（assets/、3dmol/）
"""

from pathlib import Path

ROOT = Path(SPECPATH)

# 只读资源：程序图标 + 3Dmol.js（运行期需要写 _viewer.html，会复制到 exe 同级 3dmol/）
datas = [
    (str(ROOT / "assets" / "app.ico"), "assets"),
    (str(ROOT / "3dmol" / "3Dmol-min.js"), "3dmol"),
]

hiddenimports = [
    # QtWebEngine 相关（动态导入，显式声明更稳）
    "PySide6.QtWebEngineWidgets",
    "PySide6.QtWebEngineCore",
    "PySide6.QtWebChannel",
    "PySide6.QtPrintSupport",
    # matplotlib 的 Qt 后端
    "matplotlib.backends.backend_qtagg",
    "matplotlib.backends.backend_agg",
]

a = Analysis(
    [str(ROOT / "pdos_viewer.py")],
    pathex=[str(ROOT)],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "tkinter",
        "PyQt5", "PyQt6", "PySide2",
        "IPython", "jupyter", "notebook", "pytest",
    ],
    noarchive=False,
    optimize=0,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="PDOS-Viewer",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,            # 无控制台窗口（GUI 程序）
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=str(ROOT / "assets" / "app.ico"),
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="PDOS-Viewer",
)
