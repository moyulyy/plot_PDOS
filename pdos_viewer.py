#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PDOS · 结构分析工作台
=====================

iOS / macOS 风格 PySide6 窗口，整体采用「左 / 中 / 右」三栏结构：

    左栏 —— 项目（结构、总态密度、元素 PDOS、原子-轨道-PDOS、元素-轨道-PDOS、高斯圆滑）
    中栏 —— 结构（3Dmol.js 交互式 3D 原子结构 / DOS 曲线）
    右栏 —— 参数（能量范围、费米能级、自旋、轨道、结构显示等）

数据目录（相对本文件）：
    Data-total/    TDOS.dat、FERMI_ENERGY、BAND_CENTER
    Data-element/  元素 PDOS、元素 d 带中心、SELECTED_ATOMS_LIST_*
    Data-atom/     每个原子的 PDOS（PDOS_EIG_UP/DW_<n>.dat）

高斯圆滑页可任选上述三个目录中的一个数据文件，对指定列做一维高斯圆滑并导出。

运行：
    python pdos_viewer.py
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

from PySide6.QtCore import (QAbstractTableModel, QBuffer, QEvent, QModelIndex,
                            QObject, QPointF, QRectF, QSize, Qt, QTimer, QUrl,
                            Signal, Slot)
from PySide6.QtGui import (QBrush, QColor, QFont, QIcon, QLinearGradient,
                           QPainter, QPainterPath, QPen, QPixmap, QPolygonF)
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWidgets import (QAbstractItemView, QAbstractSpinBox,
                               QApplication, QButtonGroup, QColorDialog,
                               QComboBox, QDoubleSpinBox, QFileDialog, QFrame,
                               QGridLayout, QHBoxLayout, QHeaderView, QLabel,
                               QLineEdit, QListWidget, QListWidgetItem,
                               QMessageBox, QPushButton, QScrollArea,
                               QSizePolicy, QStackedWidget, QTableView,
                               QVBoxLayout, QWidget)

from PySide6.QtWebEngineWidgets import QWebEngineView

import matplotlib
matplotlib.use("QtAgg")
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei UI", "Microsoft YaHei",
                                           "SimHei", "PingFang SC", "Segoe UI",
                                           "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

from ui_kit import (LIGHT, DARK, THEME, Card, ComboBox, SegmentedControl,
                    Switch, apply_theme, field_label, hint_label)

# ==========================================================================
# 路径 / 常量
# ==========================================================================
# APP_DIR —— 可写目录（源码运行时为脚本目录；打包后为 exe 所在目录），
#            数据文件夹、导出文件、运行时生成的 _viewer.html 都放这里。
# RES_DIR —— 只读资源目录（打包后为 PyInstaller 解包目录 _MEIPASS）。
if getattr(sys, "frozen", False):
    APP_DIR = Path(sys.executable).resolve().parent
    RES_DIR = Path(getattr(sys, "_MEIPASS", str(APP_DIR)))
else:
    APP_DIR = Path(__file__).resolve().parent
    RES_DIR = APP_DIR

BASE_DIR = APP_DIR
DATA_TOTAL = BASE_DIR / "Data-total"
DATA_ELEMENT = BASE_DIR / "Data-element"
DATA_ATOM = BASE_DIR / "Data-atom"
THREEDMOL_JS = RES_DIR / "3dmol" / "3Dmol-min.js"
VIEWER_DIR = APP_DIR / "3dmol"
APP_ICO = RES_DIR / "assets" / "app.ico"

_viewer_dir_cache = None


def viewer_html_dir():
    """返回可写且能访问到 3Dmol-min.js 的目录（打包后资源只读，需复制一份）。"""
    global _viewer_dir_cache
    if _viewer_dir_cache is not None:
        return _viewer_dir_cache
    cands = [VIEWER_DIR, Path(tempfile.gettempdir()) / "pdos_viewer"]
    for cand in cands:
        try:
            cand.mkdir(parents=True, exist_ok=True)
            js = cand / "3Dmol-min.js"
            if not js.exists() and THREEDMOL_JS.exists():
                shutil.copyfile(THREEDMOL_JS, js)
            probe = cand / ".wtest"
            probe.write_text("1", encoding="utf-8")
            probe.unlink()
            _viewer_dir_cache = cand
            return cand
        except OSError:
            continue
    _viewer_dir_cache = Path(tempfile.gettempdir())
    return _viewer_dir_cache

# VESTA 默认元素颜色 / 共价半径（Å）
VESTA_COLORS = {
    "H": "#FFCCCC", "C": "#804929", "N": "#3050F8", "O": "#FE0300",
    "Na": "#AB5CF2", "Mg": "#8AFF00", "Al": "#BFA6A6", "Si": "#F0C8A0",
    "P": "#FF8000", "S": "#FFFF00", "Cl": "#1FF01F", "K": "#8F40D4",
    "Ca": "#3DFF00", "Ti": "#BFC2C7", "V": "#A6A6AB", "Cr": "#8A99C7",
    "Mn": "#9C7AC4", "Fe": "#E06633", "Co": "#0000AF", "Ni": "#B7BBBD",
    "Cu": "#C88020", "Zn": "#7D80B0", "Pd": "#006985", "Ag": "#C0C0C0",
    "Pt": "#D0D0E0", "Au": "#FFD123",
}
VESTA_RADII = {
    "H": 0.46, "C": 0.77, "N": 0.75, "O": 0.74, "Na": 1.02, "Mg": 0.72,
    "Al": 0.54, "Si": 1.17, "P": 1.10, "S": 1.04, "Cl": 0.99, "K": 1.38,
    "Ca": 1.00, "Ti": 0.86, "V": 0.79, "Cr": 0.94, "Mn": 0.90, "Fe": 0.83,
    "Co": 1.25, "Ni": 1.25, "Cu": 1.17, "Zn": 1.25, "Pd": 1.20, "Ag": 1.34,
    "Pt": 1.30, "Au": 1.34,
}

# 六个标准视角
VIEW_ITEMS = [("front", "正视"), ("back", "后视"), ("top", "俯视"),
              ("bottom", "仰视"), ("right", "右视"), ("left", "左视")]
VIEW_HINTS = {
    "front": "a-c 面平行屏幕, a → 屏幕右",
    "back": "a-c 面平行屏幕, a → 屏幕左",
    "top": "a-b 面平行屏幕, a → 屏幕右",
    "bottom": "a-b 面平行屏幕, a → 屏幕左",
    "right": "b-c 面平行屏幕, b → 屏幕右",
    "left": "b-c 面平行屏幕, b → 屏幕左",
}


def _fallback_color(elem):
    """未知元素的稳定配色（按名称哈希，保证不太暗）。"""
    import hashlib
    h = int(hashlib.md5(elem.encode("utf-8")).hexdigest()[:6], 16)
    r = max((h >> 16) & 0xFF, 80)
    g = max((h >> 8) & 0xFF, 80)
    b = max(h & 0xFF, 80)
    return f"#{r:02X}{g:02X}{b:02X}"


def elem_color(elem):
    return VESTA_COLORS.get(elem, _fallback_color(elem))


def elem_radius(elem):
    return VESTA_RADII.get(elem, 1.2)


# 元素 PDOS 多条曲线的默认配色（需彼此区分）
ELEM_DOS_PALETTE = ("#1F77B4", "#D62728", "#2CA02C", "#9467BD",
                    "#FF7F0E", "#17BECF", "#8C564B", "#E377C2",
                    "#7F7F7F", "#BCBD22")

# 原子-轨道-PDOS 各条曲线的默认配色（共 13 条，需彼此区分）
ATOM_DOS_PALETTE = ("#1F77B4", "#D62728", "#2CA02C", "#9467BD", "#FF7F0E",
                    "#17BECF", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22",
                    "#4C72B0", "#DD8452", "#937860")

# Data-atom/PDOS_EIG_UP|DW_<n>.dat.bak 的列（顺序与文件表头一致）
ATOM_DOS_COLS = ("s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz",
                 "dx2-y2", "tot", "sum_s", "sum_p", "sum_d")


def make_random_app_icon(size=256):
    """随机生成一个个性化应用图标：渐变圆角底 + 随机化学主题图案。"""
    import random
    rng = random.Random()
    base_hue = rng.randint(0, 359)
    hue2 = (base_hue + rng.randint(40, 150)) % 360

    pm = QPixmap(size, size)
    pm.fill(Qt.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing)

    margin = size * 0.05
    rect = QRectF(margin, margin, size - 2 * margin, size - 2 * margin)
    radius = size * 0.22
    path = QPainterPath()
    path.addRoundedRect(rect, radius, radius)

    grad = QLinearGradient(rect.topLeft(), rect.bottomRight())
    grad.setColorAt(0.0, QColor.fromHsv(base_hue, 185, 238))
    grad.setColorAt(1.0, QColor.fromHsv(hue2, 205, 185))
    p.fillPath(path, QBrush(grad))

    hl = QLinearGradient(rect.topLeft(),
                         QPointF(rect.left(), rect.center().y()))
    hl.setColorAt(0.0, QColor(255, 255, 255, 70))
    hl.setColorAt(1.0, QColor(255, 255, 255, 0))
    p.fillPath(path, QBrush(hl))

    cx, cy = rect.center().x(), rect.center().y()
    r = size * 0.25
    white = QColor(255, 255, 255, 235)
    p.setPen(QPen(white, max(2.0, size * 0.024), Qt.SolidLine,
                  Qt.RoundCap, Qt.RoundJoin))
    p.setBrush(Qt.NoBrush)
    motif = rng.choice(["atom", "hex", "orbit", "dots"])
    if motif == "atom":
        p.setBrush(white)
        p.drawEllipse(QPointF(cx, cy), size * 0.075, size * 0.075)
        p.setBrush(Qt.NoBrush)
        for ang in (rng.uniform(0, 180), rng.uniform(0, 180)):
            p.save()
            p.translate(cx, cy)
            p.rotate(ang)
            p.drawEllipse(QPointF(0, 0), r, r * 0.42)
            p.restore()
    elif motif == "hex":
        poly = QPolygonF()
        for i in range(6):
            a = math.radians(60 * i - 90)
            poly.append(QPointF(cx + r * math.cos(a), cy + r * math.sin(a)))
        p.drawPolygon(poly)
        p.setBrush(white)
        p.drawEllipse(QPointF(cx, cy), size * 0.06, size * 0.06)
    elif motif == "orbit":
        p.drawEllipse(QPointF(cx, cy), r, r)
        p.setBrush(white)
        p.drawEllipse(QPointF(cx + r * 0.8, cy - r * 0.5),
                      size * 0.07, size * 0.07)
        p.drawEllipse(QPointF(cx - r * 0.7, cy + r * 0.55),
                      size * 0.05, size * 0.05)
    else:
        p.setBrush(white)
        for dx, dy, rr in ((0.0, 0.0, 0.09), (-0.45, -0.42, 0.06),
                           (0.48, -0.36, 0.05), (0.34, 0.5, 0.07)):
            p.drawEllipse(QPointF(cx + dx * size, cy + dy * size),
                          rr * size, rr * size)

    p.setPen(QPen(QColor(255, 255, 255, 70), size * 0.012))
    p.setBrush(Qt.NoBrush)
    p.drawRoundedRect(rect.adjusted(1, 1, -1, -1), radius, radius)
    p.end()

    icon = QIcon()
    for s in (16, 24, 32, 48, 64, 128, 256):
        icon.addPixmap(pm.scaled(s, s, Qt.KeepAspectRatio,
                                 Qt.SmoothTransformation))
    return icon


def write_ico(pm, path, sizes=(16, 24, 32, 48, 64, 128, 256)):
    """把 QPixmap 写成多分辨率 .ico（PNG 负载，Vista+）。"""
    import struct
    pngs = []
    for s in sizes:
        img = pm.scaled(s, s, Qt.KeepAspectRatio,
                        Qt.SmoothTransformation).toImage()
        buf = QBuffer()
        buf.open(QBuffer.WriteOnly)
        img.save(buf, "PNG")
        pngs.append((s, bytes(buf.data())))
    n = len(pngs)
    header = struct.pack("<HHH", 0, 1, n)
    offset = 6 + 16 * n
    entries = b""
    for s, data in pngs:
        wh = 0 if s >= 256 else s
        entries += struct.pack("<BBBBHHII", wh, wh, 0, 0, 1, 32,
                               len(data), offset)
        offset += len(data)
    with open(path, "wb") as f:
        f.write(header)
        f.write(entries)
        for _, data in pngs:
            f.write(data)


# ---- 基于晶胞矢量的六个标准视角 (四元数) ----
def _unit(v):
    v = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(v))
    return np.array([1.0, 0.0, 0.0]) if n < 1e-12 else v / n


def _quat_from_matrix(m):
    m = np.asarray(m, dtype=float)
    t = float(m[0, 0] + m[1, 1] + m[2, 2])
    if t > 0.0:
        s = math.sqrt(t + 1.0) * 2.0
        w = 0.25 * s
        x = (m[2, 1] - m[1, 2]) / s
        y = (m[0, 2] - m[2, 0]) / s
        z = (m[1, 0] - m[0, 1]) / s
    elif m[0, 0] > m[1, 1] and m[0, 0] > m[2, 2]:
        s = math.sqrt(1.0 + m[0, 0] - m[1, 1] - m[2, 2]) * 2.0
        w = (m[2, 1] - m[1, 2]) / s
        x = 0.25 * s
        y = (m[0, 1] + m[1, 0]) / s
        z = (m[0, 2] + m[2, 0]) / s
    elif m[1, 1] > m[2, 2]:
        s = math.sqrt(1.0 + m[1, 1] - m[0, 0] - m[2, 2]) * 2.0
        w = (m[0, 2] - m[2, 0]) / s
        x = (m[0, 1] + m[1, 0]) / s
        y = 0.25 * s
        z = (m[1, 2] + m[2, 1]) / s
    else:
        s = math.sqrt(1.0 + m[2, 2] - m[0, 0] - m[1, 1]) * 2.0
        w = (m[1, 0] - m[0, 1]) / s
        x = (m[0, 2] + m[2, 0]) / s
        y = (m[1, 2] + m[2, 1]) / s
        z = 0.25 * s
    q = np.array([x, y, z, w], dtype=float)
    nrm = float(np.linalg.norm(q))
    if nrm < 1e-12:
        return [0.0, 0.0, 0.0, 1.0]
    q /= nrm
    if q[3] < 0:
        q = -q
    return [float(v) for v in q]


def view_quaternion(cell, name):
    """根据晶胞矢量 (行向量 a,b,c) 计算标准视角四元数 [qx,qy,qz,qw]。"""
    cell = np.asarray(cell, dtype=float).reshape(3, 3)
    a, b, c = cell[0], cell[1], cell[2]
    if name in ("front", "back"):
        ex, ey = _unit(a), c
    elif name in ("top", "bottom"):
        ex, ey = _unit(a), b
    else:
        ex, ey = _unit(b), c
    if name in ("back", "bottom", "left"):
        ex = -ex
    ey = np.asarray(ey, dtype=float) - float(np.dot(ey, ex)) * ex
    ey = _unit(ey)
    ez = _unit(np.cross(ex, ey))
    return _quat_from_matrix(np.array([ex, ey, ez]))


def compute_fit(atoms, cell):
    """覆盖原子 + 晶胞顶点的包围球 [cx,cy,cz,r]。"""
    pts = [(a["x"], a["y"], a["z"]) for a in atoms]
    if cell:
        for i in (0, 1):
            for j in (0, 1):
                for k in (0, 1):
                    pts.append(tuple(i * cell[0][d] + j * cell[1][d] +
                                     k * cell[2][d] for d in range(3)))
    if not pts:
        return [0.0, 0.0, 0.0, 1.0]
    p = np.array(pts, dtype=float)
    lo, hi = p.min(0), p.max(0)
    c = (lo + hi) / 2.0
    r = float(np.linalg.norm(p - c, axis=1).max())
    return [float(c[0]), float(c[1]), float(c[2]), max(r, 1.0)]


# ==========================================================================
# 数据读取
# ==========================================================================
def read_dos_file(path: Path):
    """读取 VASPKIT PDOS/TDOS 文件。

    返回 (energies, cols)，cols 为 {列名: 数值列表}。
    跳过第一行注释、第二行说明行。
    """
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        return [], {}
    header = lines[0].split()
    names = header[1:]
    energies = []
    cols = {n: [] for n in names}
    for line in lines[2:]:
        parts = line.split()
        if len(parts) < 1 + len(names):
            continue
        try:
            energies.append(float(parts[0]))
        except ValueError:
            continue
        for i, n in enumerate(names):
            try:
                cols[n].append(float(parts[1 + i]))
            except (ValueError, IndexError):
                cols[n].append(0.0)
    return energies, cols


def _natural_key(text):
    """自然排序键：PDOS_EIG_UP_2 < PDOS_EIG_UP_10。"""
    return [int(t) if t.isdigit() else t.lower()
            for t in re.split(r"(\d+)", str(text))]


def _gauss_file_key(name):
    """高斯页文件排序：UP 在前、DW 在后，其余按自然序。"""
    spin = 0 if "_UP_" in name else (1 if "_DW_" in name else 2)
    return (spin, _natural_key(name))


def gaussian_smooth(x, y, sigma):
    """一维高斯圆滑。

    sigma 与 x 同单位（能量 eV）；边界采用端点延拓（edge padding），
    避免两端因补零而被拉向 0。
    """
    n = len(y)
    if n == 0:
        return []
    arr = np.asarray(y, dtype=float)
    try:
        sig_eV = float(sigma)
    except (TypeError, ValueError):
        sig_eV = 0.0
    if sig_eV <= 0:
        return list(arr)
    if n == 1:
        return list(arr)
    dx = (float(x[-1]) - float(x[0])) / (n - 1)
    if not np.isfinite(dx) or abs(dx) < 1e-12:
        dx = 1.0
    sig = sig_eV / abs(dx)                    # 换算为“采样点数”为单位
    if not np.isfinite(sig) or sig < 1e-9:
        return list(arr)
    half = int(math.ceil(4.0 * sig))
    half = max(1, min(half, n - 1))
    k = np.exp(-0.5 * (np.arange(-half, half + 1) / sig) ** 2)
    tot = k.sum()
    if tot <= 0:
        return list(arr)
    k /= tot
    pad = np.pad(arr, (half, half), mode="edge")
    return list(np.convolve(pad, k, mode="valid"))


def _nice_range(lo, hi, pad=0.0):
    """把区间取整到“好看”的刻度上。"""
    lo, hi = float(lo), float(hi)
    if not np.isfinite(lo) or not np.isfinite(hi):
        return -1.0, 1.0
    if hi <= lo:
        lo, hi = lo - 1.0, hi + 1.0
    span = hi - lo
    lo -= span * pad
    hi += span * pad
    step = 10.0 ** math.floor(math.log10(max(span, 1e-9)))
    return (math.floor(lo / step) * step, math.ceil(hi / step) * step)


def _nice_tick(span):
    """根据区间宽度选一个合适的主刻度间距。"""
    span = abs(float(span))
    if not np.isfinite(span) or span <= 0:
        return None
    raw = span / 16.0
    mag = 10.0 ** math.floor(math.log10(raw))
    for m in (1.0, 2.0, 2.5, 5.0, 10.0):
        if raw <= m * mag:
            return m * mag
    return 10.0 * mag


def parse_fermi(path: Path) -> float:
    """从 FERMI_ENERGY 文件里取第一个浮点数。"""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        for tok in line.split():
            try:
                return float(tok)
            except ValueError:
                continue
    return 0.0


def parse_band_center(path: Path):
    """读取 d 带中心文件，返回 (rows, average)。

    rows: [(atom_id, s_center, p_center, d_center), ...]
    average: (s_avg, p_avg, d_avg) 或 None
    """
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    rows = []
    avg = None
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if parts and parts[0].isdigit() and len(parts) >= 4:
            try:
                rows.append((int(parts[0]), float(parts[1]),
                             float(parts[2]), float(parts[3])))
            except ValueError:
                pass
        elif parts and parts[0].startswith("#Average"):
            j = i + 1
            while j < len(lines):
                p = lines[j].split()
                if len(p) >= 3:
                    try:
                        avg = (float(p[0]), float(p[1]), float(p[2]))
                        break
                    except ValueError:
                        pass
                j += 1
        i += 1
    return rows, avg


def band_center(energy, dos):
    """态密度带中心（Band center）：

        ε = ∫ E·g(E) dE / ∫ g(E) dE

    使用梯形法在文件给出的整个能量区间上数值积分；g(E) 取绝对值以兼容
    下自旋 DOS 以负值存储的情况。分母接近 0（该轨道无态密度）时返回 None。
    """
    n = min(len(energy), len(dos))
    num = den = 0.0
    for i in range(n - 1):
        e0, e1 = energy[i], energy[i + 1]
        g0, g1 = abs(dos[i]), abs(dos[i + 1])
        de = e1 - e0
        num += 0.5 * (e0 * g0 + e1 * g1) * de
        den += 0.5 * (g0 + g1) * de
    return (num / den) if abs(den) > 1e-12 else None


def parse_structure(path: Path):
    """解析 SELECTED_ATOMS_LIST_*，得到原子列表。

    每个原子: {id, label, elem, x, y, z, selected}
    坐标为分数坐标（0~1）。
    """
    atoms = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 6 or not parts[0].isdigit():
            continue
        try:
            idx = int(parts[0])
            label = parts[1]
            x = float(parts[2])
            y = float(parts[3])
            z = float(parts[4])
            selected = parts[5] == "T"
        except (ValueError, IndexError):
            continue
        elem = "".join(ch for ch in label if not ch.isdigit())
        atoms.append({"id": idx, "label": label, "elem": elem,
                      "x": x, "y": y, "z": z, "selected": selected})
    return atoms


def parse_contcar(path: Path):
    """解析 VASP CONTCAR / POSCAR 结构文件。

    返回: {comment, cell, atoms, direct}
      cell:  3x3 晶格矢量（已应用缩放因子，单位 Å）
      atoms: [{elem, x, y, z}]，坐标为笛卡尔坐标（Å）
      direct: 原始坐标是否为分数坐标
    """
    lines = [ln for ln in Path(path).read_text(
        encoding="utf-8", errors="replace").splitlines() if ln.strip()]
    if len(lines) < 8:
        raise ValueError("CONTCAR 文件格式不完整")

    comment = lines[0].strip()
    scale = float(lines[1].split()[0])
    cell = []
    for i in range(2, 5):
        cell.append([float(x) for x in lines[i].split()[:3]])

    # 第 6 行：元素名（VASP5+）或原子数（旧格式）
    tokens = lines[5].split()
    if tokens and all(t.replace(".", "").replace("-", "").replace("+", "").replace("e", "").replace("E", "").isdigit() for t in tokens):
        counts = [int(t) for t in tokens]
        symbols = [f"El{i + 1}" for i in range(len(counts))]
        coord_line = 6
    else:
        symbols = tokens
        counts = [int(x) for x in lines[6].split()]
        coord_line = 7

    # 可选 Selective dynamics 行
    if coord_line < len(lines) and lines[coord_line].strip().lower().startswith("s"):
        coord_line += 1
    mode = lines[coord_line].strip().lower() if coord_line < len(lines) else "direct"
    direct = not (mode.startswith("c") or mode.startswith("k"))

    # 缩放晶格矢量（scale>0 时）
    if scale > 0:
        cell = [[scale * v for v in vec] for vec in cell]

    def frac_to_cart(f):
        return [f[0] * cell[0][i] + f[1] * cell[1][i] + f[2] * cell[2][i]
                for i in range(3)]

    atoms = []
    n = 0
    for si, cnt in enumerate(counts):
        for _ in range(cnt):
            idx = coord_line + 1 + n
            if idx >= len(lines):
                break
            parts = lines[idx].split()
            if len(parts) < 3:
                break
            try:
                f = (float(parts[0]), float(parts[1]), float(parts[2]))
            except ValueError:
                break
            xyz = frac_to_cart(f) if direct else list(f)
            atoms.append({"elem": symbols[si],
                          "x": xyz[0], "y": xyz[1], "z": xyz[2]})
            n += 1
    return {"comment": comment, "cell": cell, "atoms": atoms, "direct": direct}


def discover_elements(element_dir: Path):
    """从 Data-element 目录中的 PDOS_EIG_UP_*.dat(.bak) 发现元素列表。"""
    elems = []
    try:
        for f in sorted(Path(element_dir).glob("PDOS_EIG_UP_*.dat*")):
            name = f.name.split("PDOS_EIG_UP_", 1)[-1]
            for suf in (".dat.bak", ".dat"):
                if name.endswith(suf):
                    name = name[: -len(suf)]
                    break
            if name and name not in elems:
                elems.append(name)
    except OSError:
        pass
    return elems


def discover_atom_ids(atom_dir: Path):
    """从 Data-atom 目录中的 PDOS_EIG_UP_<n>.dat 发现原子编号列表。"""
    ids = []
    try:
        for f in Path(atom_dir).glob("PDOS_EIG_UP_*.dat"):
            name = f.stem.replace("PDOS_EIG_UP_", "")
            if name.isdigit():
                ids.append(int(name))
    except OSError:
        pass
    return sorted(ids)


def orbital_values(cols, orbital):
    """从 PDOS 列字典里提取 s / p / d / tot 轨道。"""
    if orbital == "tot":
        return cols["tot"]
    if orbital == "s":
        return cols["s"]
    if orbital == "p":
        return [cols["py"][i] + cols["pz"][i] + cols["px"][i]
                for i in range(len(cols["py"]))]
    if orbital == "d":
        keys = ("dxy", "dyz", "dz2", "dxz", "dx2-y2")
        return [sum(cols[k][i] for k in keys)
                for i in range(len(cols["dxy"]))]
    return cols.get("tot", [])


# ==========================================================================
# 3D 结构视图（3Dmol.js）
# ==========================================================================
class StructureBridge(QObject):
    """JS -> Python 桥：把点选中的原子加入分析列表。"""
    addRequested = Signal(str)

    @Slot(str)
    def requestAdd(self, payload):
        self.addRequested.emit(payload)


class StructureView(QWebEngineView):
    """中栏 3D 原子结构浏览器（球棍模型 + 点选原子 + 六视角）。"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.atoms = []
        self.cell = None
        self.fit = None
        self.bg = "#FFFFFF"
        self._ready = False
        self._pending = []
        self._poll_count = 0
        self.bridge = StructureBridge()
        self._channel = QWebChannel(self.page())
        self._channel.registerObject("bridge", self.bridge)
        self.page().setWebChannel(self._channel)
        # WebEngine 背景不透明，避免与半透明主窗口合成时拖动闪烁
        self._apply_page_bg()
        self.loadFinished.connect(self._on_load)

    def _apply_page_bg(self):
        try:
            self.page().setBackgroundColor(QColor(self.bg))
        except Exception:
            pass

    # ---- 加载结构（写文件后 load，性能/稳定性优于 setHtml） ----
    def load_structure(self, atoms, cell=None, params=None):
        self.atoms = list(atoms)
        self.cell = cell
        if params and params.get("bg"):
            self.bg = params["bg"]
        self._apply_page_bg()
        self.fit = compute_fit(self.atoms, self.cell)
        self._ready = False
        self._pending.clear()
        path = viewer_html_dir() / "_viewer.html"
        path.write_text(self._build_html(), encoding="utf-8")
        self.load(QUrl.fromLocalFile(str(path.resolve())))

    # ---- 等待页面就绪后再执行 JS ----
    def _on_load(self, ok):
        if not ok:
            return
        self._poll_count = 0
        QTimer.singleShot(120, self._poll_ready)

    def _poll_ready(self):
        if self._ready:
            return
        self._poll_count += 1

        def cb(val):
            if val:
                self._ready = True
                for js in self._pending:
                    self.page().runJavaScript(js)
                self._pending.clear()
            elif self._poll_count < 100:
                QTimer.singleShot(120, self._poll_ready)

        self.page().runJavaScript("window.ready === true", cb)

    def run(self, js, callback=None):
        if self._ready:
            if callback is not None:
                self.page().runJavaScript(js, callback)
            else:
                self.page().runJavaScript(js)
        else:
            self._pending.append(js)

    # ---- 对外操作 ----
    def set_orientation(self, quat):
        self.run("window.setOrientation(%s);" %
                 json.dumps([float(x) for x in quat]))

    def reset_view(self):
        self.run("window.resetView();")

    def set_bg(self, color):
        self.bg = color
        self._apply_page_bg()
        self.run("window.setBg('%s');" % color)

    def clear_selection(self):
        self.run("window.clearSelection();")

    # ---- 晶胞 12 条棱（真实晶格矢量） ----
    @staticmethod
    def _cell_edges(cell):
        a, b, c = cell

        def corner(i, j, k):
            return (i * a[0] + j * b[0] + k * c[0],
                    i * a[1] + j * b[1] + k * c[1],
                    i * a[2] + j * b[2] + k * c[2])

        c000, c100 = corner(0, 0, 0), corner(1, 0, 0)
        c110, c010 = corner(1, 1, 0), corner(0, 1, 0)
        c001, c101 = corner(0, 0, 1), corner(1, 0, 1)
        c111, c011 = corner(1, 1, 1), corner(0, 1, 1)
        pairs = [(c000, c100), (c100, c110), (c110, c010), (c010, c000),
                 (c001, c101), (c101, c111), (c111, c011), (c011, c001),
                 (c000, c001), (c100, c101), (c110, c111), (c010, c011)]
        return [list(p) + list(q) for p, q in pairs]

    def _build_html(self):
        atoms = self.atoms
        elems = sorted({a["elem"] for a in atoms})
        elem_styles = {el: {"color": elem_color(el), "radius": elem_radius(el)}
                       for el in elems}
        counts = {}
        for a in atoms:
            counts[a["elem"]] = counts.get(a["elem"], 0) + 1
        lines = [str(len(atoms)), "structure"]
        for a in atoms:
            lines.append("%s %.6f %.6f %.6f" % (a["elem"], a["x"], a["y"], a["z"]))
        xyz = "\n".join(lines)
        edges = self._cell_edges(self.cell) if self.cell else []
        atoms_json = json.dumps([{"elem": a["elem"], "x": a["x"],
                                  "y": a["y"], "z": a["z"]} for a in atoms])
        cell_inv = None
        if self.cell:
            try:
                cell_inv = np.linalg.inv(np.asarray(self.cell,
                                                    dtype=float)).tolist()
            except Exception:
                cell_inv = None
        fit = self.fit or [0.0, 0.0, 0.0, 1.0]
        bg = self.bg
        dark = bg != "#FFFFFF"
        cell_color = "#6C6C70" if dark else "#9A9AA2"
        popup_bg = "#2C2C2E" if dark else "#FFFFFF"
        popup_fg = "#F2F2F7" if dark else "#1C1C1E"
        popup_bd = "#48484A" if dark else "#E3E3E8"
        panel_bg = "rgba(44,44,46,0.94)" if dark else "rgba(255,255,255,0.94)"
        panel_bd = "rgba(255,255,255,0.10)" if dark else "rgba(0,0,0,0.08)"
        panel_fg = "#F2F2F7" if dark else "#1C1C1E"
        panel_sub = "#98989D" if dark else "#8A8A8E"
        hover_bg = "rgba(255,255,255,0.14)" if dark else "rgba(0,0,0,0.06)"
        sep_color = "rgba(255,255,255,0.16)" if dark else "rgba(0,0,0,0.14)"

        # 六个标准视角四元数（按晶胞矢量计算，放在中间工具栏）
        cell_for_views = self.cell or [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                       [0.0, 0.0, 1.0]]
        viewquats = {}
        for key, _ in VIEW_ITEMS:
            try:
                viewquats[key] = view_quaternion(cell_for_views, key)
            except Exception:
                viewquats[key] = [0.0, 0.0, 0.0, 1.0]

        vbtns = "".join(
            "<button data-view='%s' onclick=\"applyView('%s')\">%s</button>"
            % (key, key, label) for key, label in VIEW_ITEMS)
        views_html = ("<div id='views'>" + vbtns +
                      "<span class='sep'></span>"
                      "<button class='rst' onclick=\"resetView()\">重置</button>"
                      "</div>")

        lrows = "".join(
            "<div class='lg-row'><span class='lg-dot' style='background:%s'>"
            "</span><span class='lg-name'>%s</span>"
            "<span class='lg-cnt'>%d</span></div>"
            % (elem_color(el), el, counts.get(el, 0)) for el in elems)
        legend_html = ("<div id='legend'><div class='lg-title'>元素图例</div>"
                       + lrows + "</div>")

        head = (
            "<!DOCTYPE html><html><head><meta charset='utf-8'>"
            "<title>structure</title>"
            "<style>html,body{width:100%;height:100%;margin:0;padding:0;"
            "overflow:hidden;background:" + bg + ";}"
            "*{font-family:'Microsoft YaHei UI','Segoe UI',sans-serif;}"
            "#v canvas{display:block;}"
            "#popup{position:absolute;display:none;background:" + popup_bg + ";"
            "border:1px solid " + popup_bd + ";border-radius:10px;"
            "padding:10px 12px;font-size:12px;color:" + popup_fg + ";z-index:30;"
            "box-shadow:0 8px 28px rgba(0,0,0,0.22);max-width:230px;}"
            "#popup .ttl{font-size:13px;font-weight:600;margin-bottom:4px;}"
            "#popup .row{margin-top:3px;}"
            "#popup button{margin-top:9px;width:100%;background:#007AFF;"
            "color:#fff;border:none;border-radius:7px;padding:6px 10px;"
            "font-size:12px;cursor:pointer;}"
            "#popup button:hover{background:#2A90FF;}"
            "#views{position:absolute;top:12px;left:50%;transform:translateX(-50%);"
            "white-space:nowrap;background:" + panel_bg + ";"
            "border:1px solid " + panel_bd + ";border-radius:12px;padding:6px 8px;"
            "box-shadow:0 6px 22px rgba(0,0,0,0.16);z-index:20;font-size:0;}"
            "#views button{display:inline-block;vertical-align:middle;"
            "border:none;background:transparent;color:" + panel_fg + ";"
            "font-size:12px;line-height:1;white-space:nowrap;"
            "writing-mode:horizontal-tb;padding:6px 12px;margin:0 2px;"
            "border-radius:8px;cursor:pointer;transition:background .15s;}"
            "#views button:hover{background:" + hover_bg + ";}"
            "#views button.active{background:#007AFF;color:#fff;font-weight:600;}"
            "#views .sep{display:inline-block;vertical-align:middle;width:1px;"
            "height:16px;background:" + sep_color + ";margin:0 4px;}"
            "#legend{position:absolute;left:12px;bottom:12px;background:" + panel_bg +
            ";border:1px solid " + panel_bd + ";border-radius:12px;padding:9px 11px;"
            "box-shadow:0 6px 22px rgba(0,0,0,0.16);z-index:20;font-size:12px;"
            "color:" + panel_fg + ";max-height:46%;overflow:auto;}"
            "#legend .lg-title{font-weight:600;font-size:11px;color:" + panel_sub +
            ";letter-spacing:1px;margin-bottom:6px;}"
            "#legend .lg-row{display:flex;align-items:center;gap:7px;margin-top:3px;}"
            "#legend .lg-dot{width:11px;height:11px;border-radius:6px;"
            "border:1px solid rgba(0,0,0,0.25);box-sizing:border-box;}"
            "#legend .lg-name{color:" + panel_fg + ";}"
            "#legend .lg-cnt{margin-left:auto;color:" + panel_sub + ";font-size:11px;}"
            "</style>"
            "<script src='3Dmol-min.js'></script>"
            "<script src='qrc:///qtwebchannel/qwebchannel.js'></script>"
            "</head><body>"
            "<div id='v' style='position:absolute;left:0;top:0;width:100%;"
            "height:100%;'></div>"
            + views_html + legend_html + "<div id='popup'></div>"
        )
        script = (
            "<script>"
            "var ATOMS=" + atoms_json + ";"
            "var EDGES=" + json.dumps(edges) + ";"
            "var ELEM=" + json.dumps(elem_styles) + ";"
            "var VIEWS=" + json.dumps(viewquats) + ";"
            "var CELL_INV=" + json.dumps(cell_inv) + ";"
            "var FIT=" + json.dumps([float(v) for v in fit]) + ";"
            "var ZOOM=0.92;var ORIENT=[0,0,0,1];"
            "var SEL=[];var selShapes=[];var lastClick=false;"
            "var viewer=$3Dmol.createViewer(document.getElementById('v'),"
            "{backgroundColor:'" + bg + "',alpha:false,antialias:true});"
            "viewer.addModel(" + json.dumps(xyz) + ",'xyz');"
            "for(var el in ELEM){var s=ELEM[el];"
            "viewer.setStyle({elem:el},{sphere:{radius:s.radius*0.5,color:s.color},"
            "stick:{radius:s.radius*0.18,color:s.color}});}"
            "for(var i=0;i<EDGES.length;i++){var e=EDGES[i];"
            "viewer.addLine({start:{x:e[0],y:e[1],z:e[2]},"
            "end:{x:e[3],y:e[4],z:e[5]},color:'" + cell_color + "',"
            "dashed:true,linewidth:1});}"
            "var HALF_H=1.0;var ASPECT=1.0;"
            "function refreshAspect(){var el=document.getElementById('v');"
            "var w=(el&&el.clientWidth)||viewer.WIDTH||0;"
            "var h=(el&&el.clientHeight)||viewer.HEIGHT||0;"
            "ASPECT=(w>0&&h>0)?(w/h):1.0;}"
            "refreshAspect();"
            "viewer.setSlabAndFog=function(){var cam=this.camera;"
            "var e=Math.max(1,cam.position.z-this.rotationGroup.position.z);"
            "var rr=Math.max(1,FIT[3]);cam.ortho=true;"
            "cam.near=Math.max(0.1,e-4*rr);cam.far=e+4*rr;"
            "cam.top=HALF_H;cam.bottom=-HALF_H;"
            "cam.left=-HALF_H*ASPECT;cam.right=HALF_H*ASPECT;"
            "cam.updateProjectionMatrix();};"
            "function applyCamera(q){"
            "if(!q)q=ORIENT;"
            "refreshAspect();HALF_H=FIT[3]/ZOOM;"
            "var half=Math.tan(Math.PI/180*((viewer.camera&&viewer.camera.fov)||50)/2);"
            "var z=viewer.CAMERA_Z-FIT[3]/(half*ZOOM);"
            "viewer.setView([-FIT[0],-FIT[1],-FIT[2],z,q[0],q[1],q[2],q[3]]);"
            "viewer.setSlabAndFog();viewer.render();}"
            "applyCamera(ORIENT);"
            "var __rt=0;window.addEventListener('resize',function(){"
            "clearTimeout(__rt);__rt=setTimeout(function(){"
            "refreshAspect();viewer.setSlabAndFog();viewer.render();},120);});"
            "document.getElementById('v').addEventListener('wheel',function(ev){"
            "ev.preventDefault();ev.stopPropagation();"
            "HALF_H*=(ev.deltaY>0?1.1:0.9);"
            "HALF_H=Math.max(FIT[3]*0.05,Math.min(FIT[3]*20,HALF_H));"
            "viewer.setSlabAndFog();viewer.render();"
            "},{passive:false,capture:true});"
            "function setActiveView(name){"
            "var b=document.querySelectorAll('#views button');"
            "for(var i=0;i<b.length;i++){"
            "b[i].classList.toggle('active',b[i].getAttribute('data-view')===name);}}"
            "window.applyView=function(name){if(!VIEWS[name])return;"
            "ORIENT=VIEWS[name];applyCamera(ORIENT);setActiveView(name);};"
            "function fracOf(a){if(!CELL_INV)return [a.x,a.y,a.z];"
            "var p=[a.x,a.y,a.z],f=[0,0,0];"
            "for(var j=0;j<3;j++){var s=0;"
            "for(var i=0;i<3;i++)s+=p[i]*CELL_INV[i][j];f[j]=s;}return f;}"
            "function updateSel(){"
            "for(var i=0;i<selShapes.length;i++)viewer.removeShape(selShapes[i]);"
            "selShapes=[];"
            "for(var k=0;k<SEL.length;k++){var a=ATOMS[SEL[k]];"
            "var r=(ELEM[a.elem]?ELEM[a.elem].radius:1.2)*0.62;"
            "var sh=viewer.addSphere({center:{x:a.x,y:a.y,z:a.z},radius:r,"
            "color:'black',wireframe:true,alpha:0.9});selShapes.push(sh);}"
            "viewer.render();}"
            "function hidePopup(){document.getElementById('popup').style.display='none';}"
            "function showPopup(idx,event){var a=ATOMS[idx];var f=fracOf(a);var h='';"
            "if(SEL.length>1){h='<div class=\"ttl\">'+SEL.length+' 个原子已选</div>';"
            "h+='<button onclick=\"addSelToList()\">将选中原子添加到分析列表</button>';}"
            "else{h='<div class=\"ttl\">原子 #'+(idx+1)+'</div>';"
            "h+='<div class=\"row\">元素：'+a.elem+'</div>';"
            "h+='<div class=\"row\">分数坐标：<br>('+f[0].toFixed(4)+', '"
            "+f[1].toFixed(4)+', '+f[2].toFixed(4)+')</div>';"
            "h+='<button onclick=\"addSelToList()\">添加到指定原子分析列表</button>';}"
            "var p=document.getElementById('popup');p.innerHTML=h;"
            "p.style.display='block';var x=event.clientX+14,y=event.clientY+14;"
            "var w=p.offsetWidth||220,ht=p.offsetHeight||120;"
            "if(x+w>window.innerWidth-8)x=window.innerWidth-w-8;"
            "if(y+ht>window.innerHeight-8)y=window.innerHeight-ht-8;"
            "if(x<8)x=8;if(y<8)y=8;p.style.left=x+'px';p.style.top=y+'px';}"
            "function addSelToList(){var out=[];"
            "for(var k=0;k<SEL.length;k++){var idx=SEL[k],a=ATOMS[idx],f=fracOf(a);"
            "out.push({index:idx+1,elem:a.elem,frac:[f[0],f[1],f[2]]});}"
            "if(window.bridge)window.bridge.requestAdd(JSON.stringify(out));"
            "hidePopup();}"
            "function onAtomClick(atom,event){"
            "lastClick=true;"
            "if(!atom){if(!event.shiftKey){SEL=[];updateSel();hidePopup();}return;}"
            "var idx=(typeof atom.index==='number')?atom.index:-1;"
            "if(idx<0){for(var i=0;i<ATOMS.length;i++){var a=ATOMS[i];"
            "if(Math.abs(a.x-atom.x)<1e-6&&Math.abs(a.y-atom.y)<1e-6&&"
            "Math.abs(a.z-atom.z)<1e-6){idx=i;break;}}}"
            "if(idx<0)return;"
            "if(event.shiftKey){var pos=SEL.indexOf(idx);"
            "if(pos>=0)SEL.splice(pos,1);else SEL.push(idx);}else{SEL=[idx];}"
            "updateSel();showPopup(idx,event);}"
            "viewer.setClickable({},true,function(atom,v,event,container){"
            "onAtomClick(atom,event);});"
            "window.onAtomClick=onAtomClick;"
            "document.getElementById('v').addEventListener('click',function(e){"
            "if(!lastClick){if(!e.shiftKey){SEL=[];updateSel();hidePopup();}}"
            "lastClick=false;},false);"
            "window.setOrientation=function(q){ORIENT=q;applyCamera(q);"
            "setActiveView(null);};"
            "window.resetView=function(){ORIENT=[0,0,0,1];applyCamera(ORIENT);"
            "setActiveView(null);SEL=[];updateSel();hidePopup();};"
            "window.clearSelection=function(){SEL=[];updateSel();hidePopup();};"
            "window.setBg=function(c){viewer.setBackgroundColor(c,1.0);viewer.render();};"
            "if(typeof qt!=='undefined'&&qt.webChannelTransport){"
            "new QWebChannel(qt.webChannelTransport,function(ch){"
            "window.bridge=ch.objects.bridge;});}"
            "window.ready=true;"
            "</script></body></html>"
        )
        return head + script


# ==========================================================================
# DOS 曲线视图（matplotlib）
# ==========================================================================
class DosTableModel(QAbstractTableModel):
    """只读表格模型：展示 DOS 文件的格式化内容。"""

    def __init__(self, headers, rows, parent=None):
        super().__init__(parent)
        self._headers = list(headers)
        self._rows = rows
        self._hi = -1

    def set_highlight(self, col):
        old = self._hi
        self._hi = -1 if col is None else int(col)
        if old != self._hi:
            if old >= 0:
                self.headerDataChanged.emit(Qt.Horizontal, old, old)
            if self._hi >= 0:
                self.headerDataChanged.emit(Qt.Horizontal, self._hi, self._hi)

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self._rows)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self._headers)

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        r, c = index.row(), index.column()
        if r >= len(self._rows) or c >= len(self._headers):
            return None
        val = self._rows[r][c]
        if role == Qt.DisplayRole:
            return val if isinstance(val, str) else f"{val:.6f}"
        if role == Qt.TextAlignmentRole:
            return int(Qt.AlignRight | Qt.AlignVCenter)
        return None

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if orientation == Qt.Horizontal:
            if role == Qt.DisplayRole:
                return (self._headers[section]
                        if section < len(self._headers) else "")
            if section != self._hi:
                return None
            if role == Qt.ForegroundRole:
                return QBrush(QColor("#0A84FF"))
            if role == Qt.FontRole:
                f = QFont()
                f.setBold(True)
                return f
            return None
        if role == Qt.DisplayRole:
            return str(section + 1)
        return None


class DosView(FigureCanvasQTAgg):
    def __init__(self, parent=None):
        self.figure = Figure(figsize=(6, 4), dpi=100)
        super().__init__(self.figure)
        self.setParent(parent)
        self.ax = self.figure.add_subplot(111)
        self._dark = False

    def set_theme(self, dark):
        self._dark = dark
        self.figure.patch.set_facecolor("#1B1B1D" if dark else "#FFFFFF")
        self.ax.set_facecolor("#1B1B1D" if dark else "#FFFFFF")
        for spine in self.ax.spines.values():
            spine.set_color("#48484A" if dark else "#E3E3E8")
        self.ax.tick_params(colors="#98989D" if dark else "#6E6E73")
        self.ax.xaxis.label.set_color("#98989D" if dark else "#6E6E73")
        self.ax.yaxis.label.set_color("#98989D" if dark else "#6E6E73")
        self.ax.title.set_color("#F2F2F7" if dark else "#1C1C1E")

    def plot(self, series, *, fermi=0.0, show_fermi=True, xlim=None, ylim=None,
             xtick=None, ytick=None, show_legend=False, show_zero=True,
             title="", xlabel="Energy (eV)", ylabel="DOS (states/eV)",
             gap=0.0, vlines=None):
        self.ax.clear()
        for s in series:
            self.ax.plot(s["x"], s["y"], color=s["color"],
                         lw=s.get("lw", 1.4), ls=s.get("ls", "-"),
                         label=s["label"])
            if s.get("fill"):
                self.ax.fill_between(s["x"], s["y"], 0,
                                     color=s.get("fill_color", s["color"]),
                                     alpha=s.get("alpha", 0.18), lw=0)
        if show_zero:
            self.ax.axhline(0, color="#000000" if not self._dark else "#FFFFFF",
                            lw=1.0, zorder=3)
        if show_fermi:
            self.ax.axvline(fermi, color="#FF9500" if self._dark else "#8A8A8E",
                            ls="--", lw=1.2, zorder=2)
        if gap and gap > 0:
            for gx in (-gap / 2.0, gap / 2.0):
                self.ax.axvline(gx,
                                color="#C7C7CC" if not self._dark else "#48484A",
                                ls=":", lw=1.0, zorder=1)
        if xlim and xlim[0] is not None and xlim[1] is not None:
            self.ax.set_xlim(float(xlim[0]), float(xlim[1]))
        if ylim and ylim[0] is not None and ylim[1] is not None:
            self.ax.set_ylim(float(ylim[0]), float(ylim[1]))
        if xtick and xtick > 0:
            self.ax.xaxis.set_major_locator(MultipleLocator(xtick))
        if ytick and ytick > 0:
            self.ax.yaxis.set_major_locator(MultipleLocator(ytick))
        for vl in (vlines or []):
            self.ax.axvline(vl["x"], color=vl.get("color", "#FF3B30"),
                            lw=vl.get("lw", 1.5), ls=vl.get("ls", "--"),
                            alpha=vl.get("alpha", 0.9), zorder=4)
            if vl.get("label"):
                # 标签放在竖线右侧一点
                self.ax.annotate(vl["label"],
                                 xy=(vl["x"], vl.get("y", 0.99)),
                                 xycoords=self.ax.get_xaxis_transform(),
                                 xytext=(4, 0), textcoords="offset points",
                                 color=vl.get("color", "#FF3B30"),
                                 fontsize=vl.get("fs", 9), fontweight="bold",
                                 ha="left", va="top", zorder=5)
        self.ax.set_title(title, fontsize=11, pad=8)
        self.ax.set_xlabel(xlabel, fontsize=9)
        self.ax.set_ylabel(ylabel, fontsize=9)
        self.ax.tick_params(labelsize=8)
        if show_legend and series:
            leg = self.ax.legend(fontsize=8, frameon=False)
            if leg:
                for t in leg.get_texts():
                    t.set_color("#F2F2F7" if self._dark else "#3A3A3C")
        self.set_theme(self._dark)
        self.figure.tight_layout(pad=0.8)
        self.draw()


# ==========================================================================
# 主窗口
# ==========================================================================
class MainWindow(QWidget):
    RESIZE_MARGIN = 6

    def __init__(self):
        super().__init__()
        self.setWindowTitle("PDOS")
        self.setWindowFlags(Qt.Window | Qt.FramelessWindowHint)
        self.setAttribute(Qt.WA_TranslucentBackground)
        self.setMouseTracking(True)
        self.setMinimumSize(1180, 720)
        self.app_icon = make_random_app_icon()
        self.setWindowIcon(self.app_icon)
        self._ensure_ico_file(self.app_icon)
        self._dark = False
        _app = QApplication.instance()
        if _app is not None:
            _app.installEventFilter(self)
        self._switches = []
        self._current_key = None
        self._atom_sel = 1
        self._elements = []
        self._struct_atoms = []
        self._struct_cell = None
        self._struct_sig = None
        self._shadow_pix = None
        self._shadow_key = None
        self.analysis_atoms = []          # 指定原子分析列表
        self.tdos_style = None            # 总态密度曲线参数/样式
        self.elem_style = None            # 元素 PDOS 曲线参数/样式
        self.atom_style = None            # 原子-轨道-PDOS 曲线参数/样式
        self._atom_band_cache = {}         # 原子 -> 各曲线带中心
        self._atom_warned = False
        self.eorbit_style = None           # 元素-轨道-PDOS 曲线参数/样式
        self._elem_orbit_cache = {}
        self._eorbit_band_cache = {}
        self._eorbit_warned = False
        self.gauss_style = None            # 高斯圆滑页面参数
        self._gauss_data = None            # (x, {列名: 列表})
        self._gauss_path = None            # 当前加载的文件 Path
        self._gauss_smoothed = None        # 圆滑后的 y
        self.data_root = None
        self._img_ready = False            # 当前页面是否已有可保存的曲线
        self.paths = {"contcar": None, "atom": None,
                      "element": None, "total": None}

        self._data = self._load_data()
        self._atom_cache = {}

        self._build()
        self._select_nav("structure")

    def eventFilter(self, obj, event):
        # 禁止鼠标滚轮改数值/改下拉值：拦截后转发给所在滚动面板，保持可滚动
        if event.type() == QEvent.Wheel:
            wdg = obj
            if (isinstance(wdg, QLineEdit)
                    and isinstance(wdg.parentWidget(), QAbstractSpinBox)):
                wdg = wdg.parentWidget()
            if isinstance(wdg, (QAbstractSpinBox, QComboBox)):
                parent = obj.parentWidget()
                while parent is not None and not isinstance(parent, QScrollArea):
                    parent = parent.parentWidget()
                if parent is not None:
                    QApplication.sendEvent(parent.viewport(), event)
                return True
        return super().eventFilter(obj, event)

    def _install_wheel_guard(self):
        """给右栏所有输入框/下拉框（含 spinbox 内部 QLineEdit）装滚轮拦截。"""
        for sp in self.right_area.findChildren(QAbstractSpinBox):
            sp.installEventFilter(self)
            le = sp.lineEdit()
            if le is not None:
                le.installEventFilter(self)
        for cb in self.right_area.findChildren(QComboBox):
            cb.installEventFilter(self)

    @staticmethod
    def _ensure_ico_file(icon):
        """首次运行生成 assets/app.ico，供打包/任务栏使用。"""
        ico = APP_ICO
        if ico.exists():
            return
        ico = APP_DIR / "assets" / "app.ico"
        if ico.exists():
            return
        try:
            ico.parent.mkdir(parents=True, exist_ok=True)
            write_ico(icon.pixmap(256, 256), ico)
        except Exception:
            pass

    # ------------------------------------------------------------------
    # data 根目录与自动发现
    # ------------------------------------------------------------------
    def _set_data_root(self, root):
        """在给定根目录下自动查找 Data-atom / Data-element / Data-total
        三个子目录以及 CONTCAR / POSCAR 文件。"""
        root = Path(root)
        self.data_root = root
        try:
            entries = list(root.iterdir())
        except OSError:
            entries = []

        def find_dir(name):
            for p in entries:
                if p.is_dir() and p.name.lower() == name.lower():
                    return p
            return None

        self.paths["atom"] = find_dir("Data-atom")
        self.paths["element"] = find_dir("Data-element")
        self.paths["total"] = find_dir("Data-total")

        contcar = None
        for nm in ("CONTCAR", "POSCAR"):
            p = root / nm
            if p.is_file():
                contcar = p
                break
        self.paths["contcar"] = contcar
        return self.paths

    # ------------------------------------------------------------------
    # 数据装载
    # ------------------------------------------------------------------
    def _load_data(self):
        data = {}
        total = self.paths.get("total")
        element = self.paths.get("element")
        atom = self.paths.get("atom")

        fermi_file = (total / "FERMI_ENERGY") if total else None
        data["fermi"] = (parse_fermi(fermi_file)
                         if (fermi_file and fermi_file.exists()) else 0.0)

        # 1) 初步读 SELECTED_ATOMS_LIST（用于元素兜底）
        prelim_atoms = []
        if element:
            for cand in sorted(element.glob("SELECTED_ATOMS_LIST_*")):
                if "p-dbandcenter" in cand.name:
                    continue
                prelim_atoms = parse_structure(cand)
                break

        # 2) 从 Data-element 自动发现元素
        self._elements = discover_elements(element) if element else []
        if not self._elements and prelim_atoms:
            self._elements = sorted({a["elem"] for a in prelim_atoms})

        # 3) 解析元素 PDOS / d 带中心，并识别“金属”元素（d 带中心非零）
        data["elem"] = {}
        metal_elems = []
        if element:
            for elem in self._elements:
                up_p = element / f"PDOS_EIG_UP_{elem}.dat.bak"
                dw_p = element / f"PDOS_EIG_DW_{elem}.dat.bak"
                if not up_p.exists():
                    up_p = element / f"PDOS_EIG_UP_{elem}.dat"
                if not dw_p.exists():
                    dw_p = element / f"PDOS_EIG_DW_{elem}.dat"
                entry = {}
                if up_p.exists():
                    entry["up"] = read_dos_file(up_p)
                if dw_p.exists():
                    entry["down"] = read_dos_file(dw_p)
                bc = element / f"BAND_CENTER_p-dbandcenter_{elem}"
                if bc.exists():
                    rows, avg = parse_band_center(bc)
                    entry["band"] = (rows, avg)
                    if avg and abs(avg[2]) > 1e-6:
                        metal_elems.append(elem)
                data["elem"][elem] = entry

        # 4) SELECTED_ATOMS_LIST：优先金属元素对应的文件（否则取第一个）
        atoms_file = None
        if element:
            for elem in metal_elems:
                cand = element / f"SELECTED_ATOMS_LIST_{elem}"
                if cand.exists():
                    atoms_file = cand
                    break
            if atoms_file is None:
                for cand in sorted(element.glob("SELECTED_ATOMS_LIST_*")):
                    if "p-dbandcenter" in cand.name:
                        continue
                    atoms_file = cand
                    break
        data["atoms"] = parse_structure(atoms_file) if atoms_file else []
        data["selected_elems"] = set(metal_elems) if metal_elems else {
            a["elem"] for a in data["atoms"] if a["selected"]}
        if not data["selected_elems"]:
            data["selected_elems"] = set(self._elements)
        data["atom_ids"] = discover_atom_ids(atom) if atom else []

        tdos_file = (total / "TDOS.dat") if total else None
        if tdos_file and tdos_file.exists():
            data["tdos"] = read_dos_file(tdos_file)

        bc_total = (total / "BAND_CENTER") if total else None
        if bc_total and bc_total.exists():
            data["band_total"] = parse_band_center(bc_total)

        # CONTCAR 结构
        data["contcar"] = None
        contcar = self.paths.get("contcar")
        if contcar and Path(contcar).exists():
            try:
                data["contcar"] = parse_contcar(Path(contcar))
            except Exception:
                data["contcar"] = None
        return data

    def _get_atom_dos(self, n):
        if n in self._atom_cache:
            return self._atom_cache[n]
        entry = {}
        atom_dir = self.paths.get("atom")
        if atom_dir:
            for spin, tag in (("up", "UP"), ("down", "DW")):
                # 优先 .dat.bak（含 sum_s / sum_p / sum_d），否则退回 .dat
                p = Path(atom_dir) / f"PDOS_EIG_{tag}_{n}.dat.bak"
                if not p.exists():
                    p = Path(atom_dir) / f"PDOS_EIG_{tag}_{n}.dat"
                if p.exists():
                    entry[spin] = read_dos_file(p)
        self._atom_cache[n] = entry
        return entry

    def _atom_label(self, n):
        for a in self._data.get("atoms", []):
            if a["id"] == n:
                return a["label"]
        contcar = self._data.get("contcar")
        if contcar and 1 <= n <= len(contcar["atoms"]):
            return f"{contcar['atoms'][n - 1]['elem']}{n}"
        return str(n)

    def _atom_entries(self):
        """原子选择下拉框条目 [(id, label), ...]。"""
        entries = [(a["id"], a["label"]) for a in self._data.get("atoms", [])]
        if entries:
            return entries
        ids = self._data.get("atom_ids", [])
        contcar = self._data.get("contcar")
        out = []
        for i, n in enumerate(ids):
            if contcar and i < len(contcar["atoms"]):
                label = f"{contcar['atoms'][i]['elem']}{n}"
            else:
                label = str(n)
            out.append((n, label))
        return out

    def _atom_detail(self, n):
        for a in self._data.get("atoms", []):
            if a["id"] == n:
                return {"elem": a["elem"], "coord": (a["x"], a["y"], a["z"]),
                        "selected": a["selected"], "frac": True}
        contcar = self._data.get("contcar")
        if contcar and 1 <= n <= len(contcar["atoms"]):
            a = contcar["atoms"][n - 1]
            return {"elem": a["elem"], "coord": (a["x"], a["y"], a["z"]),
                    "selected": a["elem"] in self._data.get("selected_elems", set()),
                    "frac": False}
        return None

    # ==================================================================
    # 界面搭建
    # ==================================================================
    def _build(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(18, 18, 18, 18)
        self.outer = outer

        self.container = QFrame()
        self.container.setObjectName("Window")
        outer.addWidget(self.container)

        root = QVBoxLayout(self.container)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addWidget(self._build_titlebar())

        body = QHBoxLayout()
        body.setContentsMargins(0, 0, 0, 0)
        body.setSpacing(0)
        body.addWidget(self._build_sidebar())

        # 右栏参数（定宽 + 独立滚动），与中栏页面处于同一行，
        # 使其高度与「中间可视化」区域完全一致（底部对齐）。
        self.right_area, self.right_lay = self._side_panel()
        # 右栏底部固定一个「保存当前图片」按钮（仅 DOS 曲线类页面显示），
        # 固定在右下角，不随右栏滚动。
        self.right_footer = self._build_right_footer()
        right_col = QWidget()
        rcol = QVBoxLayout(right_col)
        rcol.setContentsMargins(0, 0, 0, 0)
        rcol.setSpacing(0)
        rcol.addWidget(self.right_area, 1)
        rcol.addWidget(self.right_footer, 0)

        self.pages = QStackedWidget()
        self.pages.addWidget(self._page_structure())
        self.pages.addWidget(self._page_dos())

        mid = QHBoxLayout()
        mid.setContentsMargins(0, 0, 0, 0)
        mid.setSpacing(0)
        mid.addWidget(self.pages, 1)
        mid.addWidget(right_col, 0)

        right = QVBoxLayout()
        right.setContentsMargins(0, 0, 0, 0)
        right.setSpacing(0)
        right.addWidget(self._build_actionbar())
        right.addLayout(mid, 1)
        right.addWidget(self._build_footer())

        wrap = QWidget()
        wrap.setLayout(right)
        body.addWidget(wrap, 1)
        root.addLayout(body, 1)

    def _render_shadow(self, w, h, dark):
        """将柔和阴影预渲染到 pixmap（避免每帧重画 16 层圆角矩形）。"""
        pm = QPixmap(max(1, w), max(1, h))
        pm.fill(Qt.transparent)
        p = QPainter(pm)
        p.setRenderHint(QPainter.Antialiasing)
        m = 18
        rect = QRectF(m - 4, m - 2, w - 2 * m + 8, h - 2 * m + 8)
        base = QColor(0, 0, 0)
        layers = 16
        for i in range(layers, 0, -1):
            alpha = int(2.6 * (layers - i) / layers * (2.2 if dark else 1.6))
            grow = i * 0.9
            p.setBrush(QColor(base.red(), base.green(), base.blue(), max(0, alpha)))
            p.setPen(Qt.NoPen)
            p.drawRoundedRect(rect.adjusted(-grow, -grow + 2, grow, grow + 2),
                              16 + grow / 2, 16 + grow / 2)
        p.end()
        return pm

    def paintEvent(self, event):  # noqa: N802
        if self.isMaximized() or self.isFullScreen():
            super().paintEvent(event)
            return
        key = (self.width(), self.height(), self._dark)
        if self._shadow_pix is None or self._shadow_key != key:
            self._shadow_pix = self._render_shadow(self.width(), self.height(),
                                                   self._dark)
            self._shadow_key = key
        p = QPainter(self)
        p.drawPixmap(0, 0, self._shadow_pix)
        p.end()
        super().paintEvent(event)

    # ---- 标题栏 ----
    def _build_titlebar(self):
        bar = QWidget()
        bar.setObjectName("TitleBar")
        bar.setFixedHeight(48)
        bar.mousePressEvent = self._titlebar_press
        bar.mouseDoubleClickEvent = self._titlebar_double
        lay = QHBoxLayout(bar)
        lay.setContentsMargins(16, 0, 16, 0)
        lay.setSpacing(8)

        for name, slot in (("TrafficClose", self.close),
                           ("TrafficMin", self.showMinimized),
                           ("TrafficMax", self._toggle_max)):
            btn = QPushButton()
            btn.setObjectName(name)
            btn.setFixedSize(12, 12)
            btn.setCursor(Qt.PointingHandCursor)
            btn.clicked.connect(slot)
            lay.addWidget(btn)

        lay.addStretch(1)
        self.icon_label = QLabel()
        self.icon_label.setFixedSize(22, 22)
        self.icon_label.setPixmap(self.app_icon.pixmap(22, 22))
        lay.addWidget(self.icon_label)
        lay.addStretch(1)

        self.theme_btn = QPushButton("🌙")
        self.theme_btn.setObjectName("Ghost")
        self.theme_btn.setFixedWidth(38)
        self.theme_btn.setCursor(Qt.PointingHandCursor)
        self.theme_btn.setToolTip("切换深色 / 浅色主题")
        self.theme_btn.clicked.connect(self._toggle_theme)
        lay.addWidget(self.theme_btn)
        return bar

    def _titlebar_press(self, event):
        if event.button() == Qt.LeftButton and self.windowHandle():
            self.windowHandle().startSystemMove()
            event.accept()

    def _titlebar_double(self, event):
        if event.button() == Qt.LeftButton:
            self._toggle_max()

    def _toggle_max(self):
        self.showNormal() if self.isMaximized() else self.showMaximized()

    def changeEvent(self, event):  # noqa: N802
        super().changeEvent(event)
        if event.type() == QEvent.WindowStateChange:
            m = 0 if (self.isMaximized() or self.isFullScreen()) else 18
            self.outer.setContentsMargins(m, m, m, m)
            self.update()

    # ---- 侧边栏（项目） ----
    def _nav_items(self):
        return [("structure", "🧬 晶体结构"),
                ("tdos", "📊 总态密度 TDOS"),
                ("element", "🧪 元素 PDOS"),
                ("atom", "⚛️ 原子-轨道-PDOS"),
                ("eorbit", "🔬 元素-轨道-PDOS"),
                ("gauss", "🧮 高斯圆滑处理")]

    def _nav_title(self, key):
        for k, t in self._nav_items():
            if k == key:
                return t
        return key

    def _dos_title(self, key):
        if key == "tdos":
            return "总态密度 TDOS"
        if key == "element":
            return "元素 PDOS"
        if key == "atom":
            return "原子-轨道-PDOS"
        if key == "eorbit":
            return "元素-轨道-PDOS"
        if key == "gauss":
            return "高斯圆滑处理"
        return key

    def _populate_sidebar(self, lay):
        sub = QLabel("PDOS · TDOS · d 带中心")
        sub.setObjectName("SidebarSub")
        lay.addWidget(sub)
        lay.addSpacing(16)

        self.nav_group = QButtonGroup(self)
        self.nav_group.setExclusive(True)
        self.nav_btns = {}
        for key, text in self._nav_items():
            btn = QPushButton(text)
            btn.setObjectName("NavItem")
            btn.setCheckable(True)
            btn.setCursor(Qt.PointingHandCursor)
            btn.clicked.connect(lambda _=False, k=key: self._select_nav(k))
            self.nav_group.addButton(btn)
            lay.addWidget(btn)
            self.nav_btns[key] = btn

        lay.addStretch(1)
        self.lbl_sidebar_info = QLabel("在结构页选择 data 所在路径即可自动加载")
        self.lbl_sidebar_info.setObjectName("SidebarSub")
        self.lbl_sidebar_info.setWordWrap(True)
        lay.addWidget(self.lbl_sidebar_info)

    def _build_sidebar(self):
        side = QFrame()
        side.setObjectName("Sidebar")
        side.setFixedWidth(224)
        self.sidebar_lay = QVBoxLayout(side)
        self.sidebar_lay.setContentsMargins(16, 18, 16, 16)
        self.sidebar_lay.setSpacing(4)
        self._populate_sidebar(self.sidebar_lay)
        return side

    def _rebuild_sidebar(self):
        while self.sidebar_lay.count():
            item = self.sidebar_lay.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()
            else:
                sub = item.layout()
                if sub is not None:
                    while sub.count():
                        si = sub.takeAt(0)
                        sw = si.widget()
                        if sw is not None:
                            sw.deleteLater()
        self._populate_sidebar(self.sidebar_lay)

    # ---- 动作栏 ----
    def _build_actionbar(self):
        bar = QWidget()
        lay = QHBoxLayout(bar)
        lay.setContentsMargins(24, 16, 24, 8)
        lay.setSpacing(10)

        self.page_title = QLabel("晶体结构")
        self.page_title.setObjectName("H1")
        lay.addWidget(self.page_title)
        lay.addStretch(1)
        return bar

    # ---- 页脚 ----
    def _build_footer(self):
        foot = QWidget()
        lay = QVBoxLayout(foot)
        lay.setContentsMargins(24, 6, 24, 16)
        lay.setSpacing(6)
        self.status = QLabel("请点击右侧「选择 data 所在路径」加载数据")
        self.status.setObjectName("Footer")
        self.status.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self.status.setMinimumWidth(0)
        lay.addWidget(self.status)
        return foot

    # ---- 右栏参数容器 ----
    def _build_right_footer(self):
        """右栏底部固定按钮：保存当前图片（仅 DOS 曲线类页面显示）。"""
        box = QWidget()
        box.setFixedWidth(336)
        lay = QVBoxLayout(box)
        lay.setContentsMargins(0, 6, 12, 12)
        self.btn_save_img = QPushButton("保存当前图片")
        self.btn_save_img.setObjectName("Secondary")
        self.btn_save_img.setCursor(Qt.PointingHandCursor)
        self.btn_save_img.setToolTip("把当前页面的曲线图保存为 300 dpi 图片")
        self.btn_save_img.clicked.connect(self._save_figure)
        lay.addWidget(self.btn_save_img)
        box.setVisible(False)
        return box

    def _side_panel(self, width=336):
        area = QScrollArea()
        area.setWidgetResizable(True)
        area.setFrameShape(QFrame.NoFrame)
        area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        area.setFixedWidth(width)
        inner = QWidget()
        inner.setObjectName("PageInner")
        inner.setAttribute(Qt.WA_StyledBackground, True)
        area.setWidget(inner)
        lay = QVBoxLayout(inner)
        lay.setContentsMargins(0, 4, 12, 12)
        lay.setSpacing(12)
        return area, lay

    def _clear_layout(self, lay):
        while lay.count():
            item = lay.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()
            else:
                sub = item.layout()
                if sub is not None:
                    self._clear_layout(sub)

    # ==================================================================
    # 页面：结构 / DOS
    # ==================================================================
    def _page_structure(self):
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(20, 2, 8, 12)
        v.setSpacing(10)

        holder = QFrame()
        holder.setObjectName("Card")
        hv = QVBoxLayout(holder)
        hv.setContentsMargins(8, 8, 8, 8)
        hv.setSpacing(8)
        self.viewer = StructureView()
        self.viewer.bridge.addRequested.connect(self._on_add_atoms)
        hv.addWidget(self.viewer, 1)
        v.addWidget(holder, 1)

        hint = QLabel("左键拖动旋转 · 滚轮缩放 · 右键平移 · "
                      "点击原子选中（Shift 多选，点空白取消）")
        hint.setObjectName("Hint")
        hint.setAlignment(Qt.AlignCenter)
        v.addWidget(hint)
        return page

    def _page_dos(self):
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(20, 2, 8, 12)
        v.setSpacing(10)
        self._dos_page_lay = v

        holder = QFrame()
        holder.setObjectName("Card")
        hv = QVBoxLayout(holder)
        hv.setContentsMargins(12, 12, 12, 12)
        self.dos_view = DosView()
        hv.addWidget(self.dos_view, 1)
        self._dos_holder = holder
        v.addWidget(holder, 2)

        # 元素 PDOS：底部 1/3 曲线样式面板（其他页面隐藏）
        self.dos_bottom = QFrame()
        self.dos_bottom.setObjectName("Card")
        self.dos_bottom_lay = QVBoxLayout(self.dos_bottom)
        self.dos_bottom_lay.setContentsMargins(12, 10, 12, 12)
        self.dos_bottom_lay.setSpacing(6)
        self.dos_bottom.hide()
        v.addWidget(self.dos_bottom, 1)
        return page

    def _set_dos_split(self, bottom_mode):
        """元素/原子-轨道-PDOS：底栏按内容高度，作图占剩余空间。"""
        if bottom_mode:
            self.dos_bottom.show()
            self._dos_page_lay.setStretchFactor(self._dos_holder, 1)
            self._dos_page_lay.setStretchFactor(self.dos_bottom, 0)
        else:
            self.dos_bottom.hide()
            self._dos_page_lay.setStretchFactor(self._dos_holder, 1)
            self._dos_page_lay.setStretchFactor(self.dos_bottom, 0)

    # ==================================================================
    # 导航切换
    # ==================================================================
    def _select_nav(self, key):
        if key not in {k for k, _ in self._nav_items()}:
            key = "structure"
        self._current_key = key
        for k, btn in self.nav_btns.items():
            btn.setChecked(k == key)
        # 先重建右栏参数（会创建读取参数所需的控件），再载入页面
        self._rebuild_right(key)
        if key == "structure":
            self.pages.setCurrentIndex(0)
            self.page_title.setText("晶体结构")
            self._load_structure_page()
        else:
            self.pages.setCurrentIndex(1)
            self._set_dos_split(key in ("element", "atom", "eorbit", "gauss"))
            if key == "element":
                self._build_elem_style_panel()
            elif key == "atom":
                self._atom_warned = False
                self._build_atom_style_panel()
            elif key == "eorbit":
                self._eorbit_warned = False
                self._build_eorbit_style_panel()
            elif key == "gauss":
                self._gauss_ensure_loaded()
                self._build_gauss_table_panel()
            self.page_title.setText(self._nav_title(key))
            self._load_dos_page(key)

    def _build_structure(self):
        """组装当前结构（笛卡尔坐标原子列表 + 晶胞）。

        优先使用 CONTCAR；没有 CONTCAR 时退回到 SELECTED_ATOMS_LIST
        （分数坐标 ×10 作为演示尺度）。
        """
        contcar = self._data.get("contcar")
        if contcar:
            atoms = [{"elem": a["elem"], "x": a["x"], "y": a["y"],
                      "z": a["z"]} for a in contcar["atoms"]]
            return atoms, contcar["cell"]
        atoms = [{"elem": a["elem"], "x": a["x"] * 10,
                  "y": a["y"] * 10, "z": a["z"] * 10}
                 for a in self._data.get("atoms", [])]
        return atoms, None

    def _load_structure_page(self, force=False):
        atoms, cell = self._build_structure()
        self._struct_atoms = atoms
        self._struct_cell = cell
        if not atoms:
            self.status.setText("请先点击右侧「选择 data 所在路径」加载数据")
            return
        bg = "#1B1B1D" if self._dark else "#FFFFFF"
        h = hashlib.md5()
        h.update(("%d|%s|%s" % (len(atoms), cell, bg)).encode("utf-8"))
        for a in atoms:
            h.update(("%s%.4f%.4f%.4f"
                      % (a["elem"], a["x"], a["y"], a["z"])).encode("utf-8"))
        sig = h.hexdigest()
        if force or sig != self._struct_sig:
            self._struct_sig = sig
            self.viewer.load_structure(atoms, cell=cell, params={"bg": bg})
        src = "CONTCAR" if self._data.get("contcar") else "SELECTED_ATOMS_LIST"
        self.status.setText(f"结构已载入（{src}）：共 {len(atoms)} 个原子")

    # ---- DOS 数据组装 ----
    def _load_dos_page(self, key):
        if key == "tdos":
            self._plot_tdos()
            return
        if key == "element":
            self._plot_elem()
            return
        if key == "atom":
            self._plot_atom()
            return
        if key == "eorbit":
            self._plot_eorbit()
            return
        if key == "gauss":
            self._plot_gauss()
            return
        params = self._read_dos_params(key)
        series = self._build_series(key, params)
        fermi = params.get("fermi", self._data["fermi"])
        e_range = (params.get("emin", self._energy_range()[0]),
                   params.get("emax", self._energy_range()[1]))
        title = self._dos_title(key)
        self.dos_view.set_theme(self._dark)
        self.dos_view.plot(series, fermi=fermi,
                           show_fermi=params.get("show_fermi", True),
                           xlim=e_range, show_legend=True, show_zero=False,
                           title=title, xlabel="Energy (eV)",
                           ylabel="DOS (states/eV)")
        self.status.setText(f"{title}：{len(series)} 条曲线")

    # ---- 总态密度（TDOS）专用 ----
    def _tdos_style(self):
        if getattr(self, "tdos_style", None) is None:
            self.tdos_style = self._default_tdos_style()
        return self.tdos_style

    def _default_tdos_style(self):
        return {"xmin": -8.0, "xmax": 8.0, "xtick": 1.0,
                "ymin": -150.0, "ymax": 150.0, "ytick": 30.0,
                "fill": True, "line_color": "#000000",
                "fill_color": "#B0B0B0",
                "show_s": False, "show_p": False, "show_d": False,
                "bc_lw": 1.5, "bc_ls": "--", "bc_color": "#FF3B30",
                "bc_alpha": 0.9, "bc_label": True, "bc_fs": 9.0}

    def _plot_tdos(self):
        style = self._tdos_style()
        e, cols = self._data.get("tdos", ([], {}))
        line = style.get("line_color", "#000000")
        fillc = style.get("fill_color", "#B0B0B0")
        fill = style.get("fill", True)
        series = []
        for key in ("TDOS-UP", "TDOS-DOWN"):
            if e and key in cols:
                series.append({"x": list(e), "y": list(cols[key]),
                               "color": line, "fill_color": fillc,
                               "label": key, "fill": fill, "alpha": 0.45})
        # 带中心参考线（数据来自 Data-total/BAND_CENTER 的 #Average）
        vlines = []
        bc = self._data.get("band_total", (None, None))[1]
        if bc:
            lw = style.get("bc_lw", 1.5)
            ls = style.get("bc_ls", "--")
            col = style.get("bc_color", "#FF3B30")
            alpha = style.get("bc_alpha", 0.9)
            show_lbl = style.get("bc_label", True)
            fs = style.get("bc_fs", 9)
            for on, val, lbl in ((style.get("show_s"), bc[0], "s"),
                                 (style.get("show_p"), bc[1], "p"),
                                 (style.get("show_d"), bc[2], "d")):
                if on and val is not None:
                    txt = "\u03b5%s=%.2f" % (lbl, val) if show_lbl else ""
                    vlines.append({"x": val, "color": col, "lw": lw,
                                   "ls": ls, "alpha": alpha, "fs": fs,
                                   "label": txt})
        self._img_ready = bool(series)
        self.dos_view.set_theme(self._dark)
        self.dos_view.plot(
            series, fermi=0.0, show_fermi=True,
            xlim=(style["xmin"], style["xmax"]),
            ylim=(style["ymin"], style["ymax"]),
            xtick=style.get("xtick"), ytick=style.get("ytick"),
            show_legend=False, show_zero=True,
            title="Total Density of States",
            xlabel="E-Ef (eV)", ylabel="DOS of states (a.u.)",
            vlines=vlines)
        self.status.setText(f"总态密度：{len(series)} 条曲线")

    def _redraw_tdos(self):
        s = self._tdos_style()
        s["xmin"] = self.sp_xmin.value()
        s["xmax"] = self.sp_xmax.value()
        s["xtick"] = self.sp_xtick.value()
        s["ymin"] = self.sp_ymin.value()
        s["ymax"] = self.sp_ymax.value()
        s["ytick"] = self.sp_ytick.value()
        s["fill"] = self.sw_fill.isChecked()
        s["show_s"] = self.sw_s.isChecked()
        s["show_p"] = self.sw_p.isChecked()
        s["show_d"] = self.sw_d.isChecked()
        s["bc_lw"] = self.sp_bc_lw.value()
        s["bc_ls"] = self.cb_bc_ls.currentData() or "--"
        s["bc_alpha"] = self.sp_bc_alpha.value()
        s["bc_label"] = self.sw_bc_label.isChecked()
        s["bc_fs"] = self.sp_bc_fs.value()
        self._plot_tdos()

    # ---- 元素 PDOS 专用 ----
    def _elem_style(self):
        if getattr(self, "elem_style", None) is None:
            self.elem_style = {
                "xmin": -8.0, "xmax": 8.0, "xtick": 1.0,
                "ymin": -150.0, "ymax": 150.0, "ytick": 30.0,
                "fill": True,
                "elems": {},      # elem -> 每条曲线的样式（含带中心样式）
            }
        return self.elem_style

    @staticmethod
    def _spin_arrow(spin):
        return "\u2191" if spin == "up" else "\u2193"

    def _elem_keys(self):
        """返回 [(key, elem, spin), ...]，仅含有 tot 数据的曲线。"""
        keys = []
        for elem in self._elements:
            entry = self._data.get("elem", {}).get(elem, {})
            for spin in ("up", "down"):
                pair = entry.get(spin)
                if pair and pair[0] and pair[1].get("tot"):
                    keys.append((f"{elem}|{spin}", elem, spin))
        return keys

    def _elem_line(self, elem):
        """每个元素一份样式，上下自旋共用（颜色/线型/线宽一致）。"""
        st = self._elem_style()
        line = st["elems"].get(elem)
        if line is None:
            try:
                idx = self._elements.index(elem)
            except ValueError:
                idx = 0
            col = ELEM_DOS_PALETTE[idx % len(ELEM_DOS_PALETTE)]
            line = {"show": True, "color": col,
                    "ls": "-", "lw": 1.8, "alpha": 0.22,
                    "pick": "",
                    # 带中心颜色/标签颜色默认与曲线一致
                    "bc_color": col, "bc_ls": "--", "bc_lw": 1.5,
                    "bc_alpha": 0.9, "bc_label": True, "bc_fs": 9.0,
                    # 标签高度（坐标轴比例 0~1）：每个元素错开，避免重叠
                    "bc_y": round(0.97 - 0.075 * (idx % 6), 3)}
            st["elems"][elem] = line
        return line

    def _plot_elem(self):
        st = self._elem_style()
        series = []
        for key, elem, spin in self._elem_keys():
            line = self._elem_line(elem)
            if not line.get("show", True):
                continue
            e, cols = self._data["elem"][elem]["up" if spin == "up" else "down"]
            series.append({
                "x": list(e), "y": list(cols["tot"]),
                "color": line["color"], "ls": line["ls"], "lw": line["lw"],
                "label": f"{elem} {self._spin_arrow(spin)}",
                "fill": st.get("fill", True),
                "alpha": line.get("alpha", 0.22)})
        # 带中心参考线（每个元素可选 s/p/d/不绘制，样式各自独立）
        vlines = []
        idx_map = {"s": 0, "p": 1, "d": 2}
        for elem in self._elements:
            line = self._elem_line(elem)
            pick = line.get("pick", "")
            if pick not in idx_map:
                continue
            avg = self._data.get("elem", {}).get(elem, {}).get(
                "band", (None, None))[1]
            if not avg or avg[idx_map[pick]] is None:
                continue
            val = avg[idx_map[pick]]
            txt = f"\u03b5{pick}={val:.2f}" if line.get("bc_label", True) else ""
            # 每个元素使用各自的标签高度，避免重叠
            vlines.append({"x": val, "color": line["bc_color"],
                           "lw": line["bc_lw"], "ls": line["bc_ls"],
                           "alpha": line["bc_alpha"], "fs": line["bc_fs"],
                           "label": txt,
                           "y": line.get("bc_y", 0.97)})
        self._img_ready = bool(series)
        self.dos_view.set_theme(self._dark)
        self.dos_view.plot(
            series, fermi=0.0, show_fermi=True,
            xlim=(st["xmin"], st["xmax"]), ylim=(st["ymin"], st["ymax"]),
            xtick=st.get("xtick"), ytick=st.get("ytick"),
            show_legend=False, show_zero=True,
            title="Element Projected Density of States",
            xlabel="E-Ef (eV)", ylabel="DOS of states (a.u.)",
            vlines=vlines)
        self.status.setText(f"元素 PDOS：{len(series)} 条曲线")

    def _redraw_elem(self):
        st = self._elem_style()
        st["xmin"] = self.sp_el_xmin.value()
        st["xmax"] = self.sp_el_xmax.value()
        st["xtick"] = self.sp_el_xtick.value()
        st["ymin"] = self.sp_el_ymin.value()
        st["ymax"] = self.sp_el_ymax.value()
        st["ytick"] = self.sp_el_ytick.value()
        st["fill"] = self.sw_el_fill.isChecked()
        for elem, ctl in getattr(self, "_elem_line_ctrls", {}).items():
            line = st["elems"][elem]
            line["show"] = ctl["show"].isChecked()
            line["ls"] = ctl["ls"].currentData() or "-"
            line["lw"] = ctl["lw"].value()
            line["pick"] = ctl["bc"].currentData() or ""
            line["bc_ls"] = ctl["bc_ls"].currentData() or "--"
            line["bc_lw"] = ctl["bc_lw"].value()
            line["bc_alpha"] = ctl["bc_alpha"].value()
            line["bc_label"] = ctl["bc_label"].isChecked()
            line["bc_fs"] = ctl["bc_fs"].value()
            line["bc_y"] = ctl["bc_y"].value()
        self._plot_elem()

    def _elem_band_table(self):
        """返回 [(elem, (s, p, d) 或 None), ...]，用于信息表格。"""
        rows = []
        for elem in self._elements:
            entry = self._data.get("elem", {}).get(elem, {})
            rows.append((elem, entry.get("band", (None, None))[1]))
        return rows

    # ---- 原子-轨道-PDOS 专用 ----
    def _atom_style(self):
        if getattr(self, "atom_style", None) is None:
            self.atom_style = {
                "xmin": -8.0, "xmax": 8.0, "xtick": 1.0,
                "ymin": -10.0, "ymax": 10.0, "ytick": 2.0,
                "fill": True,
                "atom": None,     # 当前选中的原子编号
                "cols": {},       # 列名 -> 每条曲线的样式（含带中心样式）
            }
        return self.atom_style

    def _atom_index(self):
        st = self._atom_style()
        if st.get("atom") is None and self.analysis_atoms:
            st["atom"] = self.analysis_atoms[0]["index"]
        return st.get("atom")

    def _atom_line(self, col):
        """每列（s / py / ... / tot / sum_*）一份样式，上下自旋共用。"""
        st = self._atom_style()
        line = st["cols"].get(col)
        if line is None:
            idx = ATOM_DOS_COLS.index(col) if col in ATOM_DOS_COLS else 0
            c = ATOM_DOS_PALETTE[idx % len(ATOM_DOS_PALETTE)]
            line = {"show": True, "color": c, "ls": "-", "lw": 1.6,
                    "alpha": 0.20,
                    # 带中心颜色默认与曲线一致
                    "bc": False, "bc_color": c, "bc_ls": "--", "bc_lw": 1.2,
                    "bc_alpha": 0.9, "bc_label": True, "bc_fs": 9.0,
                    "bc_y": round(0.97 - 0.06 * (idx % 8), 3)}
            st["cols"][col] = line
        return line

    @staticmethod
    def _atom_keys(entry):
        """该原子文件中实际存在的列（按 ATOM_DOS_COLS 顺序）。"""
        have = set()
        for spin in ("up", "down"):
            pair = entry.get(spin)
            if pair:
                have.update(k for k, v in pair[1].items() if v)
        return [c for c in ATOM_DOS_COLS if c in have]

    def _atom_band_table(self, n=None):
        """返回 [(列名, up 带中心, down 带中心), ...]（按原子缓存）。"""
        n = n if n is not None else self._atom_index()
        if n is None:
            return []
        cache = getattr(self, "_atom_band_cache", None)
        if cache is None:
            cache = self._atom_band_cache = {}
        if n in cache:
            return cache[n]
        entry = self._get_atom_dos(n)
        rows = []
        for col in self._atom_keys(entry):
            vals = []
            for spin in ("up", "down"):
                pair = entry.get(spin)
                if pair and col in pair[1]:
                    vals.append(band_center(pair[0], pair[1][col]))
                else:
                    vals.append(None)
            rows.append((col, vals[0], vals[1]))
        cache[n] = rows
        return rows

    def _warn_no_atom(self):
        if getattr(self, "_atom_warned", False):
            return
        self._atom_warned = True
        QMessageBox.information(
            self, "请选择原子",
            "指定原子分析列表为空。\n\n"
            "请先在「晶体结构」页点选原子并加入「指定原子分析列表」，"
            "再回到「原子-轨道-PDOS」页绘图。")

    def _plot_atom(self):
        st = self._atom_style()
        n = self._atom_index()
        self.dos_view.set_theme(self._dark)
        if n is None:
            self._img_ready = False
            self.dos_view.plot(
                [], fermi=0.0, show_fermi=True,
                xlim=(st["xmin"], st["xmax"]),
                ylim=(st["ymin"], st["ymax"]),
                xtick=st.get("xtick"), ytick=st.get("ytick"),
                show_legend=False, show_zero=True,
                title="Atomic Projected Density of States",
                xlabel="E-Ef (eV)", ylabel="DOS of states (a.u.)")
            self.status.setText("原子-轨道-PDOS：请先在指定原子分析列表中选择原子")
            self._warn_no_atom()
            return
        entry = self._get_atom_dos(n)
        keys = self._atom_keys(entry)
        if not entry or not keys:
            self.status.setText(f"原子-轨道-PDOS：未找到原子 {n} 的 PDOS 文件")
            return
        series = []
        for col in keys:
            line = self._atom_line(col)
            if not line.get("show", True):
                continue
            for spin in ("up", "down"):
                pair = entry.get(spin)
                if not pair or col not in pair[1]:
                    continue
                series.append({
                    "x": list(pair[0]), "y": list(pair[1][col]),
                    "color": line["color"], "ls": line["ls"],
                    "lw": line["lw"],
                    "label": f"{col} {self._spin_arrow(spin)}",
                    "fill": st.get("fill", True),
                    "alpha": line.get("alpha", 0.2)})
        # 带中心参考线（每条曲线只画一条：取上/下自旋中带中心较高者）
        vlines = []
        table = {c: (u, d) for c, u, d in self._atom_band_table(n)}
        for col in keys:
            line = self._atom_line(col)
            if not line.get("bc"):
                continue
            u, d = table.get(col, (None, None))
            cand = [(v, s) for s, v in (("up", u), ("down", d))
                    if v is not None]
            if not cand:
                continue
            val, _spin = max(cand, key=lambda t: t[0])
            txt = (f"\u03b5{col}={val:.2f}"
                   if line.get("bc_label", True) else "")
            vlines.append({
                "x": val, "color": line["bc_color"],
                "lw": line["bc_lw"], "ls": line["bc_ls"],
                "alpha": line["bc_alpha"], "fs": line["bc_fs"],
                "label": txt, "y": line.get("bc_y", 0.97)})
        self._img_ready = bool(series)
        self.dos_view.plot(
            series, fermi=0.0, show_fermi=True,
            xlim=(st["xmin"], st["xmax"]), ylim=(st["ymin"], st["ymax"]),
            xtick=st.get("xtick"), ytick=st.get("ytick"),
            show_legend=False, show_zero=True,
            title="Atomic Projected Density of States",
            xlabel="E-Ef (eV)", ylabel="DOS of states (a.u.)",
            vlines=vlines)
        self.status.setText(
            f"原子-轨道-PDOS（{self._atom_label(n)}）：{len(series)} 条曲线")

    def _redraw_atom(self):
        st = self._atom_style()
        st["xmin"] = self.sp_at_xmin.value()
        st["xmax"] = self.sp_at_xmax.value()
        st["xtick"] = self.sp_at_xtick.value()
        st["ymin"] = self.sp_at_ymin.value()
        st["ymax"] = self.sp_at_ymax.value()
        st["ytick"] = self.sp_at_ytick.value()
        st["fill"] = self.sw_at_fill.isChecked()
        for col, ctl in getattr(self, "_atom_line_ctrls", {}).items():
            line = st["cols"][col]
            line["show"] = ctl["show"].isChecked()
            line["ls"] = ctl["ls"].currentData() or "-"
            line["lw"] = ctl["lw"].value()
            line["bc"] = ctl["bc"].isChecked()
            line["bc_ls"] = ctl["bc_ls"].currentData() or "--"
            line["bc_lw"] = ctl["bc_lw"].value()
            line["bc_alpha"] = ctl["bc_alpha"].value()
            line["bc_label"] = ctl["bc_label"].isChecked()
            line["bc_fs"] = ctl["bc_fs"].value()
            line["bc_y"] = ctl["bc_y"].value()
        self._plot_atom()
        self._refresh_atom_info()

    # ---- 元素-轨道-PDOS 专用（布局同原子-轨道-PDOS，只是按元素取数据） ----
    def _get_elem_orbit_dos(self, elem):
        """元素-轨道-PDOS：Data-element/PDOS_EIG_UP|DW_<元素>.dat.bak 的 13 列。"""
        cache = getattr(self, "_elem_orbit_cache", None)
        if cache is None:
            cache = self._elem_orbit_cache = {}
        if elem in cache:
            return cache[elem]
        entry = {}
        element_dir = self.paths.get("element")
        if element_dir:
            for spin, tag in (("up", "UP"), ("down", "DW")):
                p = Path(element_dir) / f"PDOS_EIG_{tag}_{elem}.dat.bak"
                if not p.exists():
                    p = Path(element_dir) / f"PDOS_EIG_{tag}_{elem}.dat"
                if p.exists():
                    entry[spin] = read_dos_file(p)
        cache[elem] = entry
        return entry

    def _eorbit_style(self):
        if getattr(self, "eorbit_style", None) is None:
            self.eorbit_style = {
                "xmin": -8.0, "xmax": 8.0, "xtick": 1.0,
                "ymin": -60.0, "ymax": 60.0, "ytick": 20.0,
                "fill": True,
                "elem": None,     # 当前选中的元素
                "cols": {},       # 列名 -> 每条曲线的样式
            }
        return self.eorbit_style

    def _eorbit_options(self):
        return [(e, e) for e in self._elements]

    def _eorbit_index(self):
        st = self._eorbit_style()
        if st.get("elem") is None and self._elements:
            st["elem"] = self._elements[0]
        return st.get("elem")

    def _eorbit_line(self, col):
        st = self._eorbit_style()
        line = st["cols"].get(col)
        if line is None:
            idx = ATOM_DOS_COLS.index(col) if col in ATOM_DOS_COLS else 0
            c = ATOM_DOS_PALETTE[idx % len(ATOM_DOS_PALETTE)]
            line = {"show": True, "color": c, "ls": "-", "lw": 1.6,
                    "alpha": 0.20,
                    # 带中心颜色默认与曲线一致
                    "bc": False, "bc_color": c, "bc_ls": "--", "bc_lw": 1.2,
                    "bc_alpha": 0.9, "bc_label": True, "bc_fs": 9.0,
                    "bc_y": round(0.97 - 0.06 * (idx % 8), 3)}
            st["cols"][col] = line
        return line

    def _eorbit_band_table(self, elem=None):
        """返回 [(列名, up 带中心, down 带中心), ...]（按元素缓存）。"""
        elem = elem if elem is not None else self._eorbit_index()
        if elem is None:
            return []
        cache = getattr(self, "_eorbit_band_cache", None)
        if cache is None:
            cache = self._eorbit_band_cache = {}
        if elem in cache:
            return cache[elem]
        entry = self._get_elem_orbit_dos(elem)
        rows = []
        for col in self._atom_keys(entry):
            vals = []
            for spin in ("up", "down"):
                pair = entry.get(spin)
                if pair and col in pair[1]:
                    vals.append(band_center(pair[0], pair[1][col]))
                else:
                    vals.append(None)
            rows.append((col, vals[0], vals[1]))
        cache[elem] = rows
        return rows

    def _warn_no_elem(self):
        if getattr(self, "_eorbit_warned", False):
            return
        self._eorbit_warned = True
        QMessageBox.information(
            self, "请选择元素",
            "未发现体系中的元素。\n\n"
            "请先在「晶体结构」页选择 data 路径加载 Data-element 数据。")

    def _plot_eorbit(self):
        st = self._eorbit_style()
        elem = self._eorbit_index()
        self.dos_view.set_theme(self._dark)
        if elem is None:
            self._img_ready = False
            self.dos_view.plot(
                [], fermi=0.0, show_fermi=True,
                xlim=(st["xmin"], st["xmax"]),
                ylim=(st["ymin"], st["ymax"]),
                xtick=st.get("xtick"), ytick=st.get("ytick"),
                show_legend=False, show_zero=True,
                title="Element Orbital Projected Density of States",
                xlabel="E-Ef (eV)", ylabel="DOS of states (a.u.)")
            self.status.setText("元素-轨道-PDOS：请先选择元素")
            self._warn_no_elem()
            return
        entry = self._get_elem_orbit_dos(elem)
        keys = self._atom_keys(entry)
        if not entry or not keys:
            self.status.setText(f"元素-轨道-PDOS：未找到元素 {elem} 的 PDOS 文件")
            return
        series = []
        for col in keys:
            line = self._eorbit_line(col)
            if not line.get("show", True):
                continue
            for spin in ("up", "down"):
                pair = entry.get(spin)
                if not pair or col not in pair[1]:
                    continue
                series.append({
                    "x": list(pair[0]), "y": list(pair[1][col]),
                    "color": line["color"], "ls": line["ls"],
                    "lw": line["lw"],
                    "label": f"{col} {self._spin_arrow(spin)}",
                    "fill": st.get("fill", True),
                    "alpha": line.get("alpha", 0.2)})
        # 带中心参考线（每条曲线只画一条：取上/下自旋中带中心较高者）
        vlines = []
        table = {c: (u, d) for c, u, d in self._eorbit_band_table(elem)}
        for col in keys:
            line = self._eorbit_line(col)
            if not line.get("bc"):
                continue
            u, d = table.get(col, (None, None))
            cand = [(v, s) for s, v in (("up", u), ("down", d))
                    if v is not None]
            if not cand:
                continue
            val, _spin = max(cand, key=lambda t: t[0])
            txt = (f"\u03b5{col}={val:.2f}"
                   if line.get("bc_label", True) else "")
            vlines.append({
                "x": val, "color": line["bc_color"],
                "lw": line["bc_lw"], "ls": line["bc_ls"],
                "alpha": line["bc_alpha"], "fs": line["bc_fs"],
                "label": txt, "y": line.get("bc_y", 0.97)})
        self._img_ready = bool(series)
        self.dos_view.plot(
            series, fermi=0.0, show_fermi=True,
            xlim=(st["xmin"], st["xmax"]), ylim=(st["ymin"], st["ymax"]),
            xtick=st.get("xtick"), ytick=st.get("ytick"),
            show_legend=False, show_zero=True,
            title="Element Orbital Projected Density of States",
            xlabel="E-Ef (eV)", ylabel="DOS of states (a.u.)",
            vlines=vlines)
        self.status.setText(f"元素-轨道-PDOS（{elem}）：{len(series)} 条曲线")

    def _redraw_eorbit(self):
        st = self._eorbit_style()
        st["xmin"] = self.sp_eo_xmin.value()
        st["xmax"] = self.sp_eo_xmax.value()
        st["xtick"] = self.sp_eo_xtick.value()
        st["ymin"] = self.sp_eo_ymin.value()
        st["ymax"] = self.sp_eo_ymax.value()
        st["ytick"] = self.sp_eo_ytick.value()
        st["fill"] = self.sw_eo_fill.isChecked()
        for col, ctl in getattr(self, "_eorbit_line_ctrls", {}).items():
            line = st["cols"][col]
            line["show"] = ctl["show"].isChecked()
            line["ls"] = ctl["ls"].currentData() or "-"
            line["lw"] = ctl["lw"].value()
            line["bc"] = ctl["bc"].isChecked()
            line["bc_ls"] = ctl["bc_ls"].currentData() or "--"
            line["bc_lw"] = ctl["bc_lw"].value()
            line["bc_alpha"] = ctl["bc_alpha"].value()
            line["bc_label"] = ctl["bc_label"].isChecked()
            line["bc_fs"] = ctl["bc_fs"].value()
            line["bc_y"] = ctl["bc_y"].value()
        self._plot_eorbit()
        self._refresh_eorbit_info()

    # ==================================================================
    # 页面：高斯圆滑处理
    # ==================================================================
    def _gauss_style(self):
        if getattr(self, "gauss_style", None) is None:
            self.gauss_style = {
                "group": "total",        # total / element / atom
                "file": None,             # 选中文件（完整路径）
                "col": None,              # 被圆滑的列名
                "sigma": 0.30,            # 一维高斯 sigma（eV）
                "xmin": -8.0, "xmax": 8.0, "xtick": 1.0,
                "ymin": -150.0, "ymax": 150.0, "ytick": 30.0,
                "orig_color": "#B0B0B0",
                "smooth_color": "#0A84FF",
                "show_orig": True,
                "show_smooth": True,
            }
        return self.gauss_style

    @staticmethod
    def _gauss_groups():
        return (("Data-total", "total"),
                ("Data-element", "element"),
                ("Data-atom", "atom"))

    def _gauss_files(self, group):
        """目录下的 DOS 数据文件 [(显示名, Path), ...]；有 .dat.bak 时优先取 .bak。"""
        root = self.paths.get(group)
        if not root:
            return []
        root = Path(root)
        if not root.is_dir():
            return []
        best = {}
        for p in root.glob("*.dat*"):
            if not p.is_file():
                continue
            is_bak = p.name.endswith(".dat.bak")
            base = p.name[:-4] if is_bak else p.name
            prev = best.get(base)
            if prev is None or (is_bak and not prev.name.endswith(".dat.bak")):
                best[base] = p
        return [(n, best[n]) for n in sorted(best, key=_gauss_file_key)]

    def _gauss_ensure_loaded(self):
        """进入页面时若尚未加载任何文件，默认加载当前目录的第一个文件。"""
        if getattr(self, "_gauss_data", None):
            return
        st = self._gauss_style()
        opts = self._gauss_files(st.get("group") or "total")
        if not opts:
            return
        want = str(st.get("file") or "")
        pick = next((p for n, p in opts if str(p) == want), opts[0][1])
        self._gauss_load(pick)

    def _gauss_load(self, path=None):
        st = self._gauss_style()
        raw = str(path if path is not None else (st.get("file") or ""))
        if not raw:
            return
        p = Path(raw)
        if not p.is_file():
            self.status.setText(f"高斯圆滑：文件不存在 {p.name}")
            return
        try:
            x, cols = read_dos_file(p)
        except Exception:
            x, cols = [], {}
        self._gauss_path = p
        st["file"] = str(p)
        if x and cols:
            self._gauss_data = (x, cols)
            if st.get("col") not in cols:
                st["col"] = next(iter(cols))
            self._gauss_autoscale()
            self._gauss_apply()
            self.status.setText(
                f"高斯圆滑：已加载 {p.name}（{len(x)} 行 × {len(cols)} 列）")
        else:
            self._gauss_data = None
            self._gauss_smoothed = None
            st["col"] = None
            self.status.setText(f"高斯圆滑：无法解析 {p.name}")
        self._refresh_gauss_controls()
        self._build_gauss_table_panel()
        self._plot_gauss()

    def _gauss_autoscale(self):
        """根据新载入的数据自动设置横/纵轴范围与刻度。"""
        st = self._gauss_style()
        data = getattr(self, "_gauss_data", None)
        if not data:
            return
        x, cols = data
        col = st.get("col")
        if not x or col not in cols:
            return
        st["xmin"], st["xmax"] = _nice_range(x[0], x[-1])
        st["xtick"] = _nice_tick(st["xmax"] - st["xmin"])
        y = cols[col]
        st["ymin"], st["ymax"] = _nice_range(min(min(y), 0.0),
                                              max(max(y), 0.0), pad=0.08)
        st["ytick"] = _nice_tick(st["ymax"] - st["ymin"])

    def _gauss_apply(self):
        st = self._gauss_style()
        data = getattr(self, "_gauss_data", None)
        self._gauss_smoothed = None
        if not data:
            return
        x, cols = data
        col = st.get("col")
        if col not in cols:
            return
        self._gauss_smoothed = gaussian_smooth(
            x, cols[col], st.get("sigma", 0.3))

    def _plot_gauss(self):
        st = self._gauss_style()
        self.dos_view.set_theme(self._dark)
        data = getattr(self, "_gauss_data", None)
        col = st.get("col")
        kw = dict(show_legend=False, show_fermi=False, show_zero=True,
                  xlim=(st["xmin"], st["xmax"]),
                  ylim=(st["ymin"], st["ymax"]),
                  xtick=st.get("xtick"), ytick=st.get("ytick"),
                  title="Gaussian Smoothing",
                  xlabel="Energy (eV)", ylabel="DOS (a.u.)")
        if not data or col not in data[1]:
            self.dos_view.plot([], **kw)
            return
        x, cols = data
        sig = float(st.get("sigma", 0.3))
        series = []
        if st.get("show_orig", True):
            series.append({"x": x, "y": cols[col],
                           "color": st.get("orig_color", "#B0B0B0"),
                           "lw": 1.0, "ls": "-", "label": "original"})
        if st.get("show_smooth", True) and self._gauss_smoothed:
            series.append({"x": x, "y": self._gauss_smoothed,
                           "color": st.get("smooth_color", "#0A84FF"),
                           "lw": 1.8, "ls": "-",
                           "label": f"quasi_dos_sigma={sig:g}"})
        kw["show_legend"] = True
        self.dos_view.plot(series, **kw)

    def _redraw_gauss(self):
        st = self._gauss_style()
        st["xmin"] = self.sp_g_xmin.value()
        st["xmax"] = self.sp_g_xmax.value()
        st["xtick"] = self.sp_g_xtick.value()
        st["ymin"] = self.sp_g_ymin.value()
        st["ymax"] = self.sp_g_ymax.value()
        st["ytick"] = self.sp_g_ytick.value()
        st["sigma"] = self.sp_g_sigma.value()
        st["show_orig"] = self.sw_g_orig.isChecked()
        st["show_smooth"] = self.sw_g_smooth.isChecked()
        self._gauss_apply()
        self._plot_gauss()
        if self._gauss_smoothed:
            self.status.setText(
                f"高斯圆滑：{st.get('col')}，sigma = {st['sigma']:g} eV")

    def _gauss_export(self):
        st = self._gauss_style()
        data = getattr(self, "_gauss_data", None)
        y = getattr(self, "_gauss_smoothed", None)
        col = st.get("col")
        if not data or not y or not col:
            QMessageBox.information(self, "无可导出数据",
                                    "请先选择文件，并点击「绘制」生成圆滑曲线。")
            return
        x = data[0]
        sig = float(st.get("sigma", 0.3))
        p = getattr(self, "_gauss_path", None)
        stem = p.name[:-4] if (p and p.name.endswith(".dat.bak")) else \
            (p.stem if p else "dos")
        where = p.parent if p else BASE_DIR
        default = str(Path(where) /
                      f"{stem}_{col}_gauss_sigma{sig:g}.dat")
        path, _ = QFileDialog.getSaveFileName(
            self, "导出圆滑曲线数据", default,
            "DAT 文件 (*.dat);;所有文件 (*)")
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("# Gaussian smoothed DOS\n")
                fh.write(f"# source : {p.name if p else '-'}\n")
                fh.write(f"# column : {col}\n")
                fh.write(f"# sigma  : {sig:.6f} eV\n")
                fh.write(f"# Energy(eV)   quasi_dos_sigma={sig:.6f}\n")
                for xi, yi in zip(x, y):
                    fh.write(f"{xi:16.8f}  {yi:16.8f}\n")
        except OSError as exc:
            QMessageBox.warning(self, "导出失败", str(exc))
            return
        self.status.setText(f"已导出圆滑曲线：{Path(path).name}")

    def _refresh_gauss_controls(self):
        """把状态回写到右栏控件（不动信号）。"""
        st = self._gauss_style()
        data = getattr(self, "_gauss_data", None)
        cb = getattr(self, "cb_gauss_col", None)
        if cb is not None:
            names = list(data[1].keys()) if data else []
            cb.blockSignals(True)
            cb.clear()
            for n in names:
                cb.addItem(n, n)
            i = cb.findData(st.get("col"))
            cb.setCurrentIndex(i if i >= 0 else 0)
            cb.blockSignals(False)
            if names:
                st["col"] = cb.currentData()
        for attr, key in (("sp_g_xmin", "xmin"), ("sp_g_xmax", "xmax"),
                          ("sp_g_xtick", "xtick"), ("sp_g_ymin", "ymin"),
                          ("sp_g_ymax", "ymax"), ("sp_g_ytick", "ytick"),
                          ("sp_g_sigma", "sigma")):
            w = getattr(self, attr, None)
            if w is not None:
                w.blockSignals(True)
                w.setValue(float(st.get(key) or 0.0))
                w.blockSignals(False)
        model = getattr(self, "_gauss_model", None)
        if model is not None:
            headers = model._headers
            try:
                model.set_highlight(headers.index(st["col"])
                                    if st.get("col") else -1)
            except ValueError:
                model.set_highlight(-1)

    def _on_gauss_group_changed(self, *_):
        st = self._gauss_style()
        st["group"] = self.cb_gauss_group.currentData()
        st["file"] = None
        st["col"] = None
        self._gauss_data = None
        self._gauss_smoothed = None
        self._gauss_path = None
        self._refresh_gauss_file_combo()
        if st.get("file"):
            self._gauss_load()
        else:
            self._refresh_gauss_controls()
            self._build_gauss_table_panel()
            self._plot_gauss()

    def _on_gauss_file_changed(self, *_):
        st = self._gauss_style()
        path = self.cb_gauss_file.currentData()
        if not path or path == str(st.get("file") or ""):
            return
        st["col"] = None
        self._gauss_load(path)

    def _on_gauss_col_changed(self, *_):
        st = self._gauss_style()
        col = self.cb_gauss_col.currentData()
        if not col or col == st.get("col"):
            return
        st["col"] = col
        self._gauss_autoscale()
        self._refresh_gauss_controls()
        self._gauss_apply()
        self._plot_gauss()

    def _refresh_gauss_file_combo(self):
        st = self._gauss_style()
        cb = getattr(self, "cb_gauss_file", None)
        if cb is None:
            return
        opts = self._gauss_files(st.get("group") or "total")
        cb.blockSignals(True)
        cb.clear()
        for name, p in opts:
            cb.addItem(name, str(p))
        i = cb.findData(str(st.get("file") or ""))
        if i < 0 and opts:
            i = 0
        if i >= 0:
            cb.setCurrentIndex(i)
        cb.blockSignals(False)
        st["file"] = cb.currentData() if opts else None
        lbl = getattr(self, "lbl_gauss_file", None)
        if lbl is not None:
            lbl.setText(f"共 {len(opts)} 个文件（优先 .dat.bak）" if opts
                        else "该目录下未找到 *.dat 数据文件")

    def _build_gauss_table_panel(self):
        """底栏：表格化展示文件内容，高亮当前被圆滑的列。"""
        self._clear_layout(self.dos_bottom_lay)
        st = self._gauss_style()
        data = getattr(self, "_gauss_data", None)
        headers, rows = ["#Energy"], []
        if data:
            x, cols = data
            headers = ["#Energy"] + list(cols.keys())
            keys = list(cols.keys())
            rows = [[x[i]] + [cols[k][i] for k in keys]
                    for i in range(min(len(x), max((len(cols[k]) for k in keys),
                                                   default=0)))]
        p = getattr(self, "_gauss_path", None)
        cap = QLabel()
        cap.setObjectName("Hint")
        if p is not None and data:
            cap.setText(f"文件内容：{p.name}　（{len(rows)} 行 × "
                        f"{max(len(headers) - 1, 0)} 列）"
                        f"　　平滑列：{st.get('col') or '—'}")
        else:
            cap.setText("文件内容：（尚未加载数据文件）")
        self.dos_bottom_lay.addWidget(cap)

        self._gauss_model = DosTableModel(headers, rows)
        tv = QTableView()
        tv.setModel(self._gauss_model)
        tv.setEditTriggers(QAbstractItemView.NoEditTriggers)
        tv.setSelectionBehavior(QAbstractItemView.SelectRows)
        tv.setSelectionMode(QAbstractItemView.ExtendedSelection)
        tv.setAlternatingRowColors(True)
        tv.setShowGrid(False)
        tv.setWordWrap(False)
        tv.verticalHeader().setVisible(False)
        tv.verticalHeader().setDefaultSectionSize(20)
        hh = tv.horizontalHeader()
        # 用前 200 行采样估算列宽（ResizeToContents 会扫描全部行，2000×14 会卡住）
        fm = tv.fontMetrics()
        sample = min(len(rows), 200)
        hh.setSectionResizeMode(QHeaderView.Interactive)
        hh.setStretchLastSection(False)
        for c, htext in enumerate(headers):
            wd = fm.horizontalAdvance(str(htext))
            for r in range(sample):
                wd = max(wd, fm.horizontalAdvance(f"{rows[r][c]:.6f}"))
            tv.setColumnWidth(c, wd + 20)
        tv.setMinimumHeight(110)
        tv.setMaximumHeight(210)
        # 表格样式跟随主题（表格不进 QSS 全局规则，避免与表头模型角色冲突）
        if self._dark:
            tv.setStyleSheet(
                "QTableView{background:#1F1F22;alternate-background-color:#26262A;"
                "color:#E5E5EA;border:1px solid #3A3A3C;border-radius:8px;"
                "selection-background-color:#0A84FF;selection-color:#FFFFFF;}"
                "QTableView::item{padding:0 4px;}"
                "QHeaderView::section{background:#2C2C2E;color:#E5E5EA;border:0;"
                "border-right:1px solid #3A3A3C;border-bottom:1px solid #3A3A3C;"
                "padding:3px 6px;}")
        else:
            tv.setStyleSheet(
                "QTableView{background:#FFFFFF;alternate-background-color:#F7F7FA;"
                "color:#1C1C1E;border:1px solid #E3E3E8;border-radius:8px;"
                "selection-background-color:#0A84FF;selection-color:#FFFFFF;}"
                "QTableView::item{padding:0 4px;}"
                "QHeaderView::section{background:#F2F2F7;color:#1C1C1E;border:0;"
                "border-right:1px solid #E3E3E8;border-bottom:1px solid #E3E3E8;"
                "padding:3px 6px;}")
        self.dos_bottom_lay.addWidget(tv, 1)
        self._gauss_table = tv
        try:
            self._gauss_model.set_highlight(
                headers.index(st["col"]) if st.get("col") else -1)
        except ValueError:
            self._gauss_model.set_highlight(-1)

    # ---- 高斯圆滑：右栏 ----
    def _right_gauss_file_card(self):
        st = self._gauss_style()
        card = Card("数据文件")
        grid = self._make_info_grid()
        grid.setVerticalSpacing(7)

        self.cb_gauss_group = ComboBox()
        for label, key in self._gauss_groups():
            self.cb_gauss_group.addItem(label, key)
        i = self.cb_gauss_group.findData(st.get("group"))
        self.cb_gauss_group.setCurrentIndex(i if i >= 0 else 0)
        grid.addWidget(field_label("目录", 76), 0, 0)
        grid.addWidget(self.cb_gauss_group, 0, 1)

        self.cb_gauss_file = ComboBox()
        self.cb_gauss_file.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        grid.addWidget(field_label("文件", 76), 1, 0)
        grid.addWidget(self.cb_gauss_file, 1, 1)

        card.add_layout(grid)
        self.lbl_gauss_file = QLabel()
        self.lbl_gauss_file.setObjectName("Hint")
        self.lbl_gauss_file.setSizePolicy(QSizePolicy.Ignored,
                                          QSizePolicy.Preferred)
        card.add(self.lbl_gauss_file)
        self.right_lay.addWidget(card)

        self._refresh_gauss_file_combo()
        self.cb_gauss_group.currentIndexChanged.connect(
            self._on_gauss_group_changed)
        self.cb_gauss_file.currentIndexChanged.connect(
            self._on_gauss_file_changed)

    def _right_gauss_param_card(self):
        st = self._gauss_style()
        card = Card("曲线参数")
        card.vbox.setContentsMargins(12, 14, 12, 16)
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(9)
        grid.setColumnStretch(1, 1)
        r = 0

        self.cb_gauss_col = ComboBox()
        self.cb_gauss_col.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        grid.addWidget(field_label("平滑列", 76), r, 0)
        grid.addWidget(self.cb_gauss_col, r, 1)
        r += 1

        self.sp_g_xmin = self._dos_spin(st["xmin"], -10000, 10000)
        self.sp_g_xmax = self._dos_spin(st["xmax"], -10000, 10000)
        grid.addWidget(field_label("横轴区间", 76), r, 0)
        grid.addLayout(self._pair(self.sp_g_xmin, self.sp_g_xmax), r, 1)
        r += 1

        self.sp_g_xtick = self._dos_spin(st["xtick"], 0.0, 5000.0, step=0.5)
        grid.addWidget(field_label("横轴刻度", 76), r, 0)
        grid.addWidget(self.sp_g_xtick, r, 1)
        r += 1

        self.sp_g_ymin = self._dos_spin(st["ymin"], -1e5, 1e5, step=5.0)
        self.sp_g_ymax = self._dos_spin(st["ymax"], -1e5, 1e5, step=5.0)
        grid.addWidget(field_label("纵轴区间", 76), r, 0)
        grid.addLayout(self._pair(self.sp_g_ymin, self.sp_g_ymax), r, 1)
        r += 1

        self.sp_g_ytick = self._dos_spin(st["ytick"], 0.0, 1e5, step=1.0)
        grid.addWidget(field_label("纵轴刻度", 76), r, 0)
        grid.addWidget(self.sp_g_ytick, r, 1)
        r += 1

        card.add_layout(grid)
        self.right_lay.addWidget(card)
        self.cb_gauss_col.currentIndexChanged.connect(self._on_gauss_col_changed)
        self._refresh_gauss_controls()

    def _right_gauss_smooth_card(self):
        st = self._gauss_style()
        card = Card("高斯圆滑")
        grid = self._make_info_grid()
        self.sp_g_sigma = self._dos_spin(st["sigma"], 0.0, 50.0, step=0.05)
        grid.addWidget(field_label("sigma (eV)", 76), 0, 0)
        grid.addWidget(self.sp_g_sigma, 0, 1)
        self.sw_g_orig = self._grid_switch(grid, 1, "显示原曲线",
                                          st.get("show_orig", True))
        self.sw_g_smooth = self._grid_switch(grid, 2, "显示圆滑曲线",
                                            st.get("show_smooth", True))
        card.add_layout(grid)
        card.add(hint_label("一维高斯核 w ∝ exp(-ΔE²/2σ²)，边界按端点延拓"))
        self.right_lay.addWidget(card)

        btn = QPushButton("绘制")
        btn.setObjectName("Primary")
        btn.setCursor(Qt.PointingHandCursor)
        btn.clicked.connect(self._redraw_gauss)
        self.right_lay.addWidget(btn)

    def _right_gauss_export_card(self):
        card = Card("导出")
        btn = QPushButton("导出圆滑曲线 (.dat)")
        btn.setObjectName("Secondary")
        btn.setCursor(Qt.PointingHandCursor)
        btn.clicked.connect(self._gauss_export)
        card.add(btn)
        card.add(hint_label("两列：Energy(eV) 与 quasi_dos_sigma=<σ>"))
        self.right_lay.addWidget(card)

    # ---- 保存当前图片（总态密度 / 元素 / 原子 / 元素-轨道） ----
    def _img_default_name(self):
        """按页面 + 当前选择生成默认文件名。"""
        key = self._current_key
        if key == "tdos":
            stem = "TDOS"
        elif key == "element":
            stem = "Element-PDOS"
        elif key == "atom":
            stem = f"Atom-Orbital-PDOS_{self._atom_index() or 'none'}"
        elif key == "eorbit":
            stem = (f"Element-Orbital-PDOS_"
                    f"{self._eorbit_index() or 'none'}")
        else:
            stem = (self._dos_title(key) or "figure")
            stem = "_".join(stem.split())
        where = self.data_root if self.data_root is not None else BASE_DIR
        return Path(where) / f"{stem}.png"

    def _save_figure(self):
        """把当前页面的曲线图保存为图片（所见即所得）。"""
        if not getattr(self, "_img_ready", False):
            QMessageBox.information(self, "无可保存的图片",
                                    "当前页面还没有绘制曲线，请先加载数据。")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "保存当前图片", str(self._img_default_name()),
            "PNG 图片 (*.png);;PDF 文档 (*.pdf);;SVG 矢量图 (*.svg);;"
            "JPEG 图片 (*.jpg *.jpeg);;所有文件 (*)")
        if not path:
            return
        # 没写扩展名时按所选过滤器的默认补 .png
        if not Path(path).suffix:
            path += ".png"
        fig = self.dos_view.figure
        try:
            fig.savefig(path, dpi=300, facecolor=fig.get_facecolor())
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self, "保存失败", str(exc))
            return
        self.status.setText(f"已保存图片：{Path(path).name}")

    # ---- 自动特征识别（总态密度信息） ----
    @staticmethod
    def _detect_gap(e, tot, ymax):
        n = len(e)
        if n == 0 or ymax <= 0:
            return "—", 0.0
        thr = 0.02 * ymax
        i0 = min(range(n), key=lambda i: abs(e[i]))
        if tot[i0] >= thr:
            return "导体", 0.0
        lo = i0
        while lo > 0 and tot[lo - 1] < thr:
            lo -= 1
        hi = i0
        while hi < n - 1 and tot[hi + 1] < thr:
            hi += 1
        gap = e[hi] - e[lo]
        return ("半导体" if gap >= 0.05 else "导体"), gap

    @staticmethod
    def _detect_magnetism(up, down):
        if not up or not down:
            return "无"
        num = sum(abs(abs(a) - abs(b)) for a, b in zip(up, down))
        den = sum(abs(a) + abs(b) for a, b in zip(up, down)) or 1.0
        return "有" if num / den > 0.05 else "无"

    def _tdos_features(self):
        e, cols = self._data.get("tdos", ([], {}))
        if not e:
            return [("状态", "未加载数据（请先选择 data 路径）")]
        up = cols.get("TDOS-UP")
        down = cols.get("TDOS-DOWN")
        n = len(e)
        tot = [abs(up[i]) + abs(down[i] if down else 0.0)
               for i in range(n)] if up else [0.0] * n
        ymax = max(tot) or 1.0
        i0 = min(range(n), key=lambda i: abs(e[i]))
        d_ef = tot[i0]
        char, gap = self._detect_gap(e, tot, ymax)
        mag = self._detect_magnetism(up, down)
        feats = [("带隙特征", char)]
        if char == "半导体":
            feats.append(("带隙大小", f"{gap:.3f} eV"))
        feats.append(("体系磁性", mag))
        feats.append(("费米能级处态密度", f"{d_ef:.3f}"))
        # 带中心（来自 Data-total/BAND_CENTER 的 #Average）
        bc = self._data.get("band_total", (None, None))[1]
        if bc:
            feats.append(("s 带中心", f"{bc[0]:.3f} eV"))
            feats.append(("p 带中心", f"{bc[1]:.3f} eV"))
            feats.append(("d 带中心", f"{bc[2]:.3f} eV"))
        return feats

    def _read_dos_params(self, key):
        e_min, e_max = self._energy_range()
        p = {"spin": "both", "orbital": "tot", "atom": self._atom_sel,
             "show_fermi": True, "emin": e_min, "emax": e_max,
             "fermi": self._data["fermi"]}
        if hasattr(self, "sp_emin"):
            p["emin"] = self.sp_emin.value()
            p["emax"] = self.sp_emax.value()
        if hasattr(self, "ed_fermi"):
            try:
                p["fermi"] = float(self.ed_fermi.text())
            except ValueError:
                pass
        if hasattr(self, "seg_spin"):
            v = self.seg_spin.value()
            if v:
                p["spin"] = v
        if hasattr(self, "cb_orbital"):
            p["orbital"] = self.cb_orbital.currentData() or "tot"
        if hasattr(self, "cb_atom"):
            p["atom"] = int(self.cb_atom.currentData() or self._atom_sel)
        if hasattr(self, "sw_fermi"):
            p["show_fermi"] = self.sw_fermi.isChecked()
        return p

    def _energy_range(self):
        if "tdos" in self._data:
            e, _ = self._data["tdos"]
            if e:
                return (min(e), max(e))
        return (-20.0, 20.0)

    def _build_series(self, key, params):
        series = []
        spin = params.get("spin", "both")
        orbital = params.get("orbital", "tot")
        up_c = "#007AFF"
        dw_c = "#FF3B30"

        def add(spin_name, e, y, color, label):
            if not e or len(e) != len(y):
                return
            series.append({"x": list(e), "y": list(y), "color": color,
                           "label": label, "fill": True, "alpha": 0.20})

        if key == "tdos":
            e, cols = self._data.get("tdos", ([], {}))
            if spin in ("both", "up") and "TDOS-UP" in cols:
                add("up", e, cols["TDOS-UP"], up_c, "TDOS ↑")
            if spin in ("both", "down") and "TDOS-DOWN" in cols:
                y = cols["TDOS-DOWN"]
                if spin == "down":
                    y = [-v for v in y]
                add("down", e, y, dw_c, "TDOS ↓")
            return series

        if key.startswith("elem_"):
            elem = key.split("_", 1)[1]
            entry = self._data["elem"].get(elem, {})
            if spin in ("both", "up") and "up" in entry:
                e, cols = entry["up"]
                add("up", e, orbital_values(cols, orbital), up_c,
                    f"{elem} {orbital} ↑")
            if spin in ("both", "down") and "down" in entry:
                e, cols = entry["down"]
                y = orbital_values(cols, orbital)
                if spin == "down":
                    y = [-v for v in y]
                add("down", e, y, dw_c, f"{elem} {orbital} ↓")
            return series

        if key == "atom":
            n = int(params.get("atom", 1))
            entry = self._get_atom_dos(n)
            label = self._atom_label(n)
            if spin in ("both", "up") and "up" in entry:
                e, cols = entry["up"]
                add("up", e, orbital_values(cols, orbital), up_c,
                    f"{label} {orbital} ↑")
            if spin in ("both", "down") and "down" in entry:
                e, cols = entry["down"]
                y = orbital_values(cols, orbital)
                if spin == "down":
                    y = [-v for v in y]
                add("down", e, y, dw_c, f"{label} {orbital} ↓")
            return series

        return series

    # ==================================================================
    # 右栏参数（按项目重建）
    # ==================================================================
    def _rebuild_right(self, key):
        self._clear_layout(self.right_lay)
        self._switches = []
        # 右栏控件已销毁，置空避免旧 C++ 对象被误用
        self.analysis_list = None
        self.lbl_analysis_count = None
        if key == "structure":
            self._right_structure()
        elif key == "tdos":
            self._right_tdos_card()
            self._right_tdos_info()
        elif key == "element":
            self._right_elem_card()
            self._right_elem_info_card()
        elif key == "atom":
            self._right_atom_pick_card()
            self._right_atom_card()
            self._right_atom_info_card()
        elif key == "eorbit":
            self._right_eorbit_pick_card()
            self._right_eorbit_card()
            self._right_eorbit_info_card()
        elif key == "gauss":
            self._right_gauss_file_card()
            self._right_gauss_param_card()
            self._right_gauss_smooth_card()
            self._right_gauss_export_card()
        # 底部固定的「保存当前图片」按钮：只在有曲线图的页面出现
        if getattr(self, "right_footer", None) is not None:
            self.right_footer.setVisible(
                key in ("tdos", "element", "atom", "eorbit"))
        self._install_wheel_guard()

    def _right_structure(self):
        self._right_paths_card()
        self._right_analysis_card()

    def _right_paths_card(self):
        card = Card("数据路径")
        btn = QPushButton("选择 data 所在路径")
        btn.setObjectName("Secondary")
        btn.setCursor(Qt.PointingHandCursor)
        btn.clicked.connect(self._pick_data_root)
        card.add(btn)
        self.lbl_paths = QLabel()
        self.lbl_paths.setObjectName("Hint")
        self.lbl_paths.setWordWrap(True)
        self.lbl_paths.setTextFormat(Qt.PlainText)
        # 长路径不得撑宽定宽右栏：忽略宽度尺寸提示，允许任意处换行
        self.lbl_paths.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self.lbl_paths.setMinimumWidth(0)
        self.lbl_paths.setTextInteractionFlags(Qt.TextSelectableByMouse)
        card.add(self.lbl_paths)
        self._update_paths_label()
        self.right_lay.addWidget(card)

    def _right_analysis_card(self):
        card = Card("指定原子分析列表")
        card.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        head = QHBoxLayout()
        head.setSpacing(8)
        self.lbl_analysis_count = QLabel("0 个原子")
        self.lbl_analysis_count.setObjectName("Chip")
        head.addWidget(self.lbl_analysis_count)
        head.addStretch(1)
        btn_clear = QPushButton("清空")
        btn_clear.setObjectName("Ghost")
        btn_clear.setCursor(Qt.PointingHandCursor)
        btn_clear.clicked.connect(self._clear_analysis)
        head.addWidget(btn_clear)
        card.add_layout(head)

        self.analysis_list = QListWidget()
        self.analysis_list.setObjectName("AnalysisList")
        self.analysis_list.setFrameShape(QFrame.NoFrame)
        self.analysis_list.setSelectionMode(QAbstractItemView.NoSelection)
        self.analysis_list.setVerticalScrollMode(
            QAbstractItemView.ScrollPerPixel)
        self.analysis_list.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.analysis_list.setSpacing(5)
        self.analysis_list.setMinimumHeight(120)
        self.analysis_list.setSizePolicy(QSizePolicy.Preferred,
                                         QSizePolicy.Expanding)
        card.body.addWidget(self.analysis_list, 1)
        self.right_lay.addWidget(card, 1)
        self._refresh_analysis_list()

    def _analysis_row(self, item):
        row = QFrame()
        row.setObjectName("AnalysisRow")
        row.setAttribute(Qt.WA_StyledBackground, True)
        row.setMinimumHeight(34)
        h = QHBoxLayout(row)
        h.setContentsMargins(9, 4, 5, 4)
        h.setSpacing(7)

        dot = QLabel()
        dot.setFixedSize(11, 11)
        dot.setStyleSheet(
            "background:%s;border-radius:5px;border:1px solid #00000033;"
            % elem_color(item["elem"]))
        h.addWidget(dot)

        name = QLabel(item["elem"])
        name.setObjectName("AnalysisName")
        h.addWidget(name)

        idx = QLabel("#%d" % item["index"])
        idx.setObjectName("AnalysisIndex")
        h.addWidget(idx)

        f = item["frac"]
        coords = QLabel("%.3f %.3f %.3f" % (f[0], f[1], f[2]))
        coords.setObjectName("AnalysisCoords")
        coords.setToolTip("分数坐标 (x, y, z)")
        h.addStretch(1)
        h.addWidget(coords)

        btn = QPushButton("\u2715")
        btn.setObjectName("AnalysisDel")
        btn.setToolTip("删除该原子")
        btn.setCursor(Qt.PointingHandCursor)
        btn.setFixedSize(22, 22)
        btn.clicked.connect(
            lambda _=False, i=item["index"]: self._remove_analysis_atom(i))
        h.addWidget(btn)
        return row

    def _right_dos_common(self):
        card = Card("曲线参数")

        # 能量范围
        e_min, e_max = self._energy_range()
        row = QHBoxLayout()
        row.setSpacing(10)
        row.addWidget(field_label("能量范围", 88))
        self.sp_emin = QDoubleSpinBox()
        self.sp_emin.setRange(-50, 50)
        self.sp_emin.setValue(round(e_min, 2))
        self.sp_emin.setDecimals(2)
        self.sp_emin.setSingleStep(0.5)
        self.sp_emax = QDoubleSpinBox()
        self.sp_emax.setRange(-50, 50)
        self.sp_emax.setValue(round(e_max, 2))
        self.sp_emax.setDecimals(2)
        self.sp_emax.setSingleStep(0.5)
        self.sp_emin.valueChanged.connect(lambda _: self._on_curve_param())
        self.sp_emax.valueChanged.connect(lambda _: self._on_curve_param())
        row.addWidget(self.sp_emin)
        row.addWidget(QLabel("~"))
        row.addWidget(self.sp_emax)
        row.addStretch(1)
        card.add_layout(row)

        # 费米能级
        row2 = QHBoxLayout()
        row2.setSpacing(10)
        row2.addWidget(field_label("费米能级", 88))
        self.ed_fermi = QLineEdit(f"{self._data['fermi']:.4f}")
        self.ed_fermi.setFixedWidth(100)
        self.ed_fermi.editingFinished.connect(self._on_curve_param)
        row2.addWidget(self.ed_fermi)
        row2.addStretch(1)
        card.add_layout(row2)

        self.sw_fermi, boxf = self._switch_box("显示费米线", True,
                                               "在曲线中绘制 E_F 虚线")
        self.sw_fermi.toggled.connect(self._on_curve_param)
        card.add_layout(boxf)

        # 自旋
        self.seg_spin = SegmentedControl([("both", "↑↓"), ("up", "↑"),
                                          ("down", "↓")], value="both")
        self.seg_spin.changed.connect(lambda _: self._on_curve_param())
        row3 = QHBoxLayout()
        row3.setSpacing(12)
        row3.addWidget(field_label("自旋", 88))
        row3.addWidget(self.seg_spin)
        row3.addStretch(1)
        card.add_layout(row3)

        # 轨道
        if self._current_key != "tdos":
            row4 = QHBoxLayout()
            row4.setSpacing(12)
            row4.addWidget(field_label("轨道", 88))
            self.cb_orbital = ComboBox()
            for label, val in (("tot", "tot"), ("s", "s"),
                               ("p", "p"), ("d", "d")):
                self.cb_orbital.addItem(label, val)
            self.cb_orbital.currentIndexChanged.connect(
                lambda _: self._on_curve_param())
            row4.addWidget(self.cb_orbital)
            row4.addStretch(1)
            card.add_layout(row4)

        self.right_lay.addWidget(card)

    # ---- 总态密度：曲线参数卡片 ----
    @staticmethod
    def _dos_spin(value, lo, hi, step=0.5):
        sp = QDoubleSpinBox()
        sp.setRange(lo, hi)
        sp.setDecimals(2)
        sp.setSingleStep(step)
        sp.setValue(float(value))
        sp.setMinimumWidth(56)
        sp.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        return sp

    @staticmethod
    def _pair(w1, w2):
        lay = QHBoxLayout()
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(6)
        lay.addWidget(w1)
        lay.addWidget(QLabel("~"), 0, Qt.AlignCenter)
        lay.addWidget(w2)
        return lay

    @staticmethod
    def _style_color_btn(btn, color):
        col = QColor(color)
        lum = col.red() * 0.299 + col.green() * 0.587 + col.blue() * 0.114
        btn.setText(str(color).upper())
        btn.setStyleSheet(
            "background:%s;color:%s;border:1px solid #00000033;"
            "border-radius:8px;padding:3px 8px;font-weight:600;"
            % (color, "#FFFFFF" if lum < 140 else "#1C1C1E"))

    def _color_button(self, color, key):
        btn = QPushButton()
        btn.setObjectName("Secondary")
        btn.setCursor(Qt.PointingHandCursor)
        btn.setFixedHeight(28)
        self._style_color_btn(btn, color)
        btn.clicked.connect(
            lambda _=False, k=key, b=btn: self._pick_color(k, b))
        return btn

    def _pick_color(self, key, btn):
        cur = self._tdos_style().get(key, "#000000")
        dlg = QColorDialog(QColor(cur), self)
        dlg.setWindowTitle("选择颜色")
        dlg.setWindowModality(Qt.ApplicationModal)
        dlg.setWindowFlag(Qt.WindowStaysOnTopHint, True)
        if dlg.exec():
            c = dlg.currentColor()
            if c.isValid():
                self._tdos_style()[key] = c.name()
                self._style_color_btn(btn, c.name())

    @staticmethod
    def _swatch_btn(btn, color):
        btn.setText("")
        btn.setStyleSheet(
            "background:%s;border:1px solid #00000033;border-radius:8px;"
            % color)

    def _line_color_button(self, key, color, which="color"):
        btn = QPushButton()
        btn.setObjectName("Secondary")
        btn.setCursor(Qt.PointingHandCursor)
        btn.setFixedHeight(26)
        self._swatch_btn(btn, color)
        btn.clicked.connect(
            lambda _=False, k=key, b=btn, w=which:
            self._pick_line_color(k, b, w))
        return btn

    def _pick_line_color(self, elem, btn, which="color"):
        line = self._elem_style()["elems"][elem]
        dlg = QColorDialog(QColor(line[which]), self)
        dlg.setWindowTitle("选择颜色")
        dlg.setWindowModality(Qt.ApplicationModal)
        dlg.setWindowFlag(Qt.WindowStaysOnTopHint, True)
        if dlg.exec():
            c = dlg.currentColor()
            if c.isValid():
                old = line.get(which)
                line[which] = c.name()
                self._swatch_btn(btn, c.name())
                # 曲线颜色改变时，带中心颜色若仍与旧曲线一致则自动跟随
                if which == "color" and line.get("bc_color") == old:
                    line["bc_color"] = c.name()
                    ctl = getattr(self, "_elem_line_ctrls", {}).get(elem)
                    bcb = (ctl or {}).get("bc_color_btn")
                    if bcb is not None:
                        self._swatch_btn(bcb, c.name())

    def _elem_bc_color_button(self, color):
        btn = QPushButton()
        btn.setObjectName("Secondary")
        btn.setCursor(Qt.PointingHandCursor)
        btn.setFixedHeight(28)
        self._style_color_btn(btn, color)
        btn.clicked.connect(lambda _=False, b=btn: self._pick_el_bc_color(b))
        return btn

    def _pick_el_bc_color(self, btn):
        st = self._elem_style()
        dlg = QColorDialog(QColor(st.get("bc_color", "#FF3B30")), self)
        dlg.setWindowTitle("选择颜色")
        dlg.setWindowModality(Qt.ApplicationModal)
        dlg.setWindowFlag(Qt.WindowStaysOnTopHint, True)
        if dlg.exec():
            c = dlg.currentColor()
            if c.isValid():
                st["bc_color"] = c.name()
                self._style_color_btn(btn, c.name())

    def _grid_switch(self, grid, row, text, value):
        sw = Switch()
        sw.setChecked(value, animate=False)
        sw.apply_theme(THEME)
        self._switches.append(sw)
        grid.addWidget(field_label(text, 76), row, 0)
        grid.addWidget(sw, row, 1, Qt.AlignLeft)
        return sw

    def _compact_switch(self, value):
        sw = Switch()
        sw.setFixedSize(QSize(38, 22))
        sw.setChecked(value, animate=False)
        sw.apply_theme(THEME)
        self._switches.append(sw)
        return sw

    @staticmethod
    def _table_sep():
        ln = QFrame()
        ln.setFixedHeight(1)
        ln.setStyleSheet("background: rgba(142,142,147,90); border: none;")
        return ln

    @staticmethod
    def _slim(combo):
        """去掉下拉框的固定最小宽度（QSS 默认 min-width 90px），可自适应列宽。"""
        combo.setStyleSheet(
            "QComboBox{min-width:0px;padding:1px 0px 1px 6px;}")
        combo.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        combo.setFixedHeight(24)
        return combo

    def _right_tdos_card(self):
        style = self._tdos_style()
        card = Card("曲线参数")
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(9)
        grid.setColumnStretch(1, 1)
        r = 0

        self.sp_xmin = self._dos_spin(style["xmin"], -10000, 10000)
        self.sp_xmax = self._dos_spin(style["xmax"], -10000, 10000)
        grid.addWidget(field_label("横轴区间", 76), r, 0)
        grid.addLayout(self._pair(self.sp_xmin, self.sp_xmax), r, 1)
        r += 1

        self.sp_xtick = self._dos_spin(style["xtick"], 0.0, 5000.0, step=0.5)
        grid.addWidget(field_label("横轴刻度", 76), r, 0)
        grid.addWidget(self.sp_xtick, r, 1)
        r += 1

        self.sp_ymin = self._dos_spin(style["ymin"], -1e5, 1e5, step=5.0)
        self.sp_ymax = self._dos_spin(style["ymax"], -1e5, 1e5, step=5.0)
        grid.addWidget(field_label("纵轴区间", 76), r, 0)
        grid.addLayout(self._pair(self.sp_ymin, self.sp_ymax), r, 1)
        r += 1

        self.sp_ytick = self._dos_spin(style["ytick"], 0.0, 1e5, step=1.0)
        grid.addWidget(field_label("纵轴刻度", 76), r, 0)
        grid.addWidget(self.sp_ytick, r, 1)
        r += 1

        self.sw_fill = Switch()
        self.sw_fill.setChecked(style.get("fill", True), animate=False)
        self.sw_fill.apply_theme(THEME)
        self._switches.append(self.sw_fill)
        grid.addWidget(field_label("是否填充", 76), r, 0)
        grid.addWidget(self.sw_fill, r, 1, Qt.AlignLeft)
        r += 1

        grid.addWidget(field_label("曲线颜色", 76), r, 0)
        self.btn_line = self._color_button(
            style.get("line_color", "#000000"), "line_color")
        grid.addWidget(self.btn_line, r, 1)
        r += 1

        grid.addWidget(field_label("填充颜色", 76), r, 0)
        self.btn_fill = self._color_button(
            style.get("fill_color", "#B0B0B0"), "fill_color")
        grid.addWidget(self.btn_fill, r, 1)
        r += 1

        # ---- 带中心参考线 ----
        sep = QLabel("带中心参考线")
        sep.setObjectName("CardTitle")
        grid.addWidget(sep, r, 0, 1, 2)
        r += 1

        self.sw_s = self._grid_switch(grid, r, "s 带中心",
                                      style.get("show_s", False))
        r += 1
        self.sw_p = self._grid_switch(grid, r, "p 带中心",
                                      style.get("show_p", False))
        r += 1
        self.sw_d = self._grid_switch(grid, r, "d 带中心",
                                      style.get("show_d", False))
        r += 1

        self.sp_bc_lw = self._dos_spin(style.get("bc_lw", 1.5), 0.5, 10.0,
                                       step=0.5)
        grid.addWidget(field_label("线宽", 76), r, 0)
        grid.addWidget(self.sp_bc_lw, r, 1)
        r += 1

        self.cb_bc_ls = ComboBox()
        for label, val in (("实线", "-"), ("虚线", "--"),
                           ("点线", ":"), ("点划线", "-.")):
            self.cb_bc_ls.addItem(label, val)
        self.cb_bc_ls.setCurrentIndex(
            max(0, self.cb_bc_ls.findData(style.get("bc_ls", "--"))))
        grid.addWidget(field_label("线型", 76), r, 0)
        grid.addWidget(self.cb_bc_ls, r, 1)
        r += 1

        grid.addWidget(field_label("颜色", 76), r, 0)
        self.btn_bc = self._color_button(
            style.get("bc_color", "#FF3B30"), "bc_color")
        grid.addWidget(self.btn_bc, r, 1)
        r += 1

        self.sp_bc_alpha = self._dos_spin(style.get("bc_alpha", 0.9),
                                          0.0, 1.0, step=0.05)
        grid.addWidget(field_label("透明度", 76), r, 0)
        grid.addWidget(self.sp_bc_alpha, r, 1)
        r += 1

        self.sw_bc_label = self._grid_switch(grid, r, "显示标签",
                                             style.get("bc_label", True))
        r += 1

        self.sp_bc_fs = self._dos_spin(style.get("bc_fs", 9), 6.0, 40.0,
                                       step=1.0)
        grid.addWidget(field_label("标签字号", 76), r, 0)
        grid.addWidget(self.sp_bc_fs, r, 1)
        r += 1

        card.add_layout(grid)
        btn = QPushButton("重新绘制")
        btn.setObjectName("Primary")
        btn.setCursor(Qt.PointingHandCursor)
        btn.clicked.connect(self._redraw_tdos)
        card.add(btn)
        self.right_lay.addWidget(card)

    def _right_tdos_info(self):
        card = Card("总态密度信息")
        card.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        grid = self._make_info_grid()
        for r, (k, v) in enumerate(self._tdos_features()):
            grid.addWidget(field_label(k, 96), r, 0)
            lab = QLabel(v)
            lab.setObjectName("Value")
            lab.setWordWrap(True)
            grid.addWidget(lab, r, 1)
        card.add_layout(grid)
        self.right_lay.addWidget(card, 1)

    # ---- 元素 PDOS：曲线参数卡片 ----
    def _right_elem_card(self):
        st = self._elem_style()
        card = Card("曲线参数")
        card.vbox.setContentsMargins(12, 14, 12, 16)
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(9)
        grid.setColumnStretch(1, 1)
        r = 0

        self.sp_el_xmin = self._dos_spin(st["xmin"], -10000, 10000)
        self.sp_el_xmax = self._dos_spin(st["xmax"], -10000, 10000)
        grid.addWidget(field_label("横轴区间", 76), r, 0)
        grid.addLayout(self._pair(self.sp_el_xmin, self.sp_el_xmax), r, 1)
        r += 1

        self.sp_el_xtick = self._dos_spin(st["xtick"], 0.0, 5000.0, step=0.5)
        grid.addWidget(field_label("横轴刻度", 76), r, 0)
        grid.addWidget(self.sp_el_xtick, r, 1)
        r += 1

        self.sp_el_ymin = self._dos_spin(st["ymin"], -1e5, 1e5, step=5.0)
        self.sp_el_ymax = self._dos_spin(st["ymax"], -1e5, 1e5, step=5.0)
        grid.addWidget(field_label("纵轴区间", 76), r, 0)
        grid.addLayout(self._pair(self.sp_el_ymin, self.sp_el_ymax), r, 1)
        r += 1

        self.sp_el_ytick = self._dos_spin(st["ytick"], 0.0, 1e5, step=1.0)
        grid.addWidget(field_label("纵轴刻度", 76), r, 0)
        grid.addWidget(self.sp_el_ytick, r, 1)
        r += 1

        self.sw_el_fill = self._grid_switch(grid, r, "是否填充",
                                            st.get("fill", True))
        r += 1

        card.add_layout(grid)
        self.right_lay.addWidget(card)

        btn = QPushButton("重新绘制")
        btn.setObjectName("Primary")
        btn.setCursor(Qt.PointingHandCursor)
        btn.clicked.connect(self._redraw_elem)
        self.right_lay.addWidget(btn)

    # ---- 元素 PDOS：底部曲线样式面板（含带中心样式） ----
    def _build_elem_style_panel(self):
        self._clear_layout(self.dos_bottom_lay)

        inner = QWidget()
        grid = QGridLayout(inner)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(4)
        grid.setVerticalSpacing(6)
        # 一行一个元素：元素 | 曲线(显示/色/型/宽) | 带中心(选择/色/型/宽/透/签/号/高)
        heads = ("元素", "显示", "颜色", "线型", "线宽",
                 "带中心", "颜色", "线型", "线宽",
                 "透明度", "标签", "字号", "标签高度")
        mins = (26, 30, 30, 42, 34, 42, 30, 42, 34, 36, 30, 34, 46)
        for c, txt in enumerate(heads):
            grid.setColumnStretch(c, 1)
            lab = QLabel(txt)
            lab.setObjectName("FieldLabel")
            lab.setMinimumWidth(mins[c])
            if c == 12:
                lab.setToolTip("标签高度：坐标轴比例 0~1（0=底部，1=顶部）\n"
                               "每个元素独立设置，避免标签彼此重叠以及与 y=0 轴重叠")
            grid.addWidget(lab, 0, c)
        grid.addWidget(self._table_sep(), 1, 0, 1, 13)

        bc_choices = (("不绘制", ""), ("s 带", "s"),
                      ("p 带", "p"), ("d 带", "d"))
        self._elem_line_ctrls = {}
        r = 2
        for i, elem in enumerate(self._elements):
            line = self._elem_line(elem)

            name = QLabel(elem)
            name.setObjectName("FieldLabel")
            name.setMinimumWidth(mins[0])
            grid.addWidget(name, r, 0)

            sw = self._compact_switch(line.get("show", True))
            grid.addWidget(sw, r, 1, Qt.AlignLeft)

            btn_c = self._line_color_button(elem, line["color"])
            btn_c.setMinimumWidth(mins[2])
            btn_c.setMaximumWidth(120)
            btn_c.setFixedHeight(24)
            btn_c.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(btn_c, r, 2)

            cb = ComboBox()
            for lab, val in (("实线", "-"), ("虚线", "--"),
                             ("点线", ":"), ("点划线", "-.")):
                cb.addItem(lab, val)
            cb.setCurrentIndex(max(0, cb.findData(line["ls"])))
            self._slim(cb)
            grid.addWidget(cb, r, 3)

            lw = self._dos_spin(line["lw"], 0.5, 10.0, step=0.5)
            lw.setMinimumWidth(mins[4])
            lw.setMaximumWidth(120)
            lw.setFixedHeight(24)
            lw.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(lw, r, 4)

            bc = ComboBox()
            for lab, val in bc_choices:
                bc.addItem(lab, val)
            bc.setCurrentIndex(max(0, bc.findData(line.get("pick", ""))))
            self._slim(bc)
            grid.addWidget(bc, r, 5)

            bcb = self._line_color_button(elem, line["bc_color"],
                                          which="bc_color")
            bcb.setMinimumWidth(mins[6])
            bcb.setMaximumWidth(120)
            bcb.setFixedHeight(24)
            bcb.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(bcb, r, 6)

            bcls = ComboBox()
            for lab, val in (("实线", "-"), ("虚线", "--"),
                             ("点线", ":"), ("点划线", "-.")):
                bcls.addItem(lab, val)
            bcls.setCurrentIndex(max(0, bcls.findData(line["bc_ls"])))
            self._slim(bcls)
            grid.addWidget(bcls, r, 7)

            bclw = self._dos_spin(line["bc_lw"], 0.5, 10.0, step=0.5)
            bclw.setMinimumWidth(mins[8])
            bclw.setMaximumWidth(120)
            bclw.setFixedHeight(24)
            bclw.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(bclw, r, 8)

            bcalpha = self._dos_spin(line["bc_alpha"], 0.0, 1.0, step=0.05)
            bcalpha.setMinimumWidth(mins[9])
            bcalpha.setMaximumWidth(120)
            bcalpha.setFixedHeight(24)
            bcalpha.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(bcalpha, r, 9)

            bclabel = self._compact_switch(line.get("bc_label", True))
            grid.addWidget(bclabel, r, 10, Qt.AlignLeft)

            bcfs = self._dos_spin(line["bc_fs"], 6.0, 40.0, step=1.0)
            bcfs.setMinimumWidth(mins[11])
            bcfs.setMaximumWidth(120)
            bcfs.setFixedHeight(24)
            bcfs.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(bcfs, r, 11)

            bcy = self._dos_spin(line.get("bc_y", 0.97), 0.0, 1.0, step=0.02)
            bcy.setButtonSymbols(QAbstractSpinBox.NoButtons)
            bcy.setMinimumWidth(mins[12])
            bcy.setMaximumWidth(120)
            bcy.setFixedHeight(24)
            bcy.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(bcy, r, 12)

            self._elem_line_ctrls[elem] = {
                "show": sw, "ls": cb, "lw": lw, "bc": bc,
                "bc_ls": bcls, "bc_lw": bclw, "bc_alpha": bcalpha,
                "bc_label": bclabel, "bc_fs": bcfs, "bc_y": bcy,
                "bc_color_btn": bcb}

            r += 1
            if i < len(self._elements) - 1:
                grid.addWidget(self._table_sep(), r, 0, 1, 13)
                grid.setRowMinimumHeight(r, 7)
                r += 1

        grid.setRowStretch(r, 1)

        if not self._elem_line_ctrls:
            hint = QLabel("未找到 PDOS_EIG_*_*.dat.bak")
            hint.setObjectName("Hint")
            hint.setWordWrap(True)
            grid.addWidget(hint, r, 0, 1, 13)

        # 最多显示约 3 个元素的高度，更多元素时纵向滚动，
        # 不再进一步挤占中间 DOS 图的高度
        scroll = QScrollArea()
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setWidget(inner)
        if len(self._elements) > 4:
            # 只保留前 4 行的高度（第 4 个元素行的底部），超出的部分滚动查看
            grid.activate()
            it = grid.itemAtPosition(2 + 2 * 3, 0)
            limit = 176
            if (it is not None and it.widget() is not None
                    and it.widget().geometry().bottom() > 0):
                limit = it.widget().geometry().bottom() + 4
            scroll.setMaximumHeight(limit)
        self.dos_bottom_lay.addWidget(scroll)

    def _right_elem_info_card(self):
        card = Card("元素 PDOS 信息")
        card.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        rows = self._elem_band_table()
        table = QGridLayout()
        table.setContentsMargins(0, 0, 0, 0)
        table.setHorizontalSpacing(6)
        table.setVerticalSpacing(0)
        if not rows:
            table.addWidget(field_label("信息", 60), 0, 0)
            lab = QLabel("未加载数据（请先选择 data 路径）")
            lab.setObjectName("Value")
            lab.setWordWrap(True)
            table.addWidget(lab, 0, 1)
            card.add_layout(table)
            card.body.addStretch(1)
            self.right_lay.addWidget(card, 1)
            return
        # 横向 s/p/d，纵向每个元素
        for c, h in enumerate(("元素", "s 带中心", "p 带中心", "d 带中心")):
            lab = QLabel(h)
            lab.setObjectName("CardTitle")
            lab.setAlignment(Qt.AlignCenter)
            lab.setContentsMargins(0, 5, 0, 5)
            table.addWidget(lab, 0, c)
        table.addWidget(self._table_sep(), 1, 0, 1, 4)
        for i, (elem, avg) in enumerate(rows):
            rr = 2 + 2 * i
            name = QLabel(elem)
            name.setObjectName("Value")
            name.setAlignment(Qt.AlignCenter)
            name.setContentsMargins(0, 5, 0, 5)
            table.addWidget(name, rr, 0)
            for c in range(3):
                v = avg[c] if avg else None
                lab = QLabel(f"{v:.2f}" if v is not None else "—")
                lab.setObjectName("Value")
                lab.setAlignment(Qt.AlignCenter)
                lab.setContentsMargins(0, 5, 0, 5)
                table.addWidget(lab, rr, c + 1)
            if i < len(rows) - 1:
                table.addWidget(self._table_sep(), rr + 1, 0, 1, 4)
        for c in range(4):
            table.setColumnStretch(c, 1)
        card.add_layout(table)
        card.body.addStretch(1)
        self.right_lay.addWidget(card, 1)

    def _right_elem_info(self, elem):
        card = Card(f"{elem} · d 带中心")
        card.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        grid = self._make_info_grid()
        entry = self._data["elem"].get(elem, {})
        avg = entry.get("band", (None, None))[1]
        rows = []
        if avg:
            rows = [("s 带中心", f"{avg[0]:.3f} eV"),
                    ("p 带中心", f"{avg[1]:.3f} eV"),
                    ("d 带中心", f"{avg[2]:.3f} eV")]
        else:
            rows = [("信息", "未找到 BAND_CENTER")]
        for r, (k, v) in enumerate(rows):
            grid.addWidget(field_label(k, 88), r, 0)
            lab = QLabel(v)
            lab.setObjectName("Value")
            grid.addWidget(lab, r, 1)
        card.add_layout(grid)
        self.right_lay.addWidget(card, 1)

    # ---- 原子-轨道-PDOS：右栏 ----
    def _right_atom_pick_card(self):
        card = Card("指定原子")
        grid = self._make_info_grid()
        grid.addWidget(field_label("原子", 76), 0, 0)
        self.cb_atom = ComboBox()
        self.cb_atom.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        for it in self.analysis_atoms:
            self.cb_atom.addItem(
                f"{it['elem']}{it['index']}  (#{it['index']})", it["index"])
        st = self._atom_style()
        if self.analysis_atoms:
            idx = self.cb_atom.findData(st.get("atom"))
            if idx < 0:
                idx = 0
            self.cb_atom.setCurrentIndex(idx)
            st["atom"] = self.cb_atom.itemData(idx)
        else:
            self.cb_atom.addItem("（指定原子分析列表为空）", None)
        self.cb_atom.currentIndexChanged.connect(self._on_atom_changed)
        grid.addWidget(self.cb_atom, 0, 1)
        card.add_layout(grid)
        card.add(hint_label("列表来自「晶体结构」页的指定原子分析列表"))
        self.right_lay.addWidget(card)

    def _right_atom_card(self):
        st = self._atom_style()
        card = Card("曲线参数")
        card.vbox.setContentsMargins(12, 14, 12, 16)
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(9)
        grid.setColumnStretch(1, 1)
        r = 0

        self.sp_at_xmin = self._dos_spin(st["xmin"], -10000, 10000)
        self.sp_at_xmax = self._dos_spin(st["xmax"], -10000, 10000)
        grid.addWidget(field_label("横轴区间", 76), r, 0)
        grid.addLayout(self._pair(self.sp_at_xmin, self.sp_at_xmax), r, 1)
        r += 1

        self.sp_at_xtick = self._dos_spin(st["xtick"], 0.0, 5000.0, step=0.5)
        grid.addWidget(field_label("横轴刻度", 76), r, 0)
        grid.addWidget(self.sp_at_xtick, r, 1)
        r += 1

        self.sp_at_ymin = self._dos_spin(st["ymin"], -1e5, 1e5, step=1.0)
        self.sp_at_ymax = self._dos_spin(st["ymax"], -1e5, 1e5, step=1.0)
        grid.addWidget(field_label("纵轴区间", 76), r, 0)
        grid.addLayout(self._pair(self.sp_at_ymin, self.sp_at_ymax), r, 1)
        r += 1

        self.sp_at_ytick = self._dos_spin(st["ytick"], 0.0, 1e5, step=0.5)
        grid.addWidget(field_label("纵轴刻度", 76), r, 0)
        grid.addWidget(self.sp_at_ytick, r, 1)
        r += 1

        self.sw_at_fill = self._grid_switch(grid, r, "是否填充",
                                            st.get("fill", True))
        r += 1

        card.add_layout(grid)
        self.right_lay.addWidget(card)

        btn = QPushButton("重新绘制")
        btn.setObjectName("Primary")
        btn.setCursor(Qt.PointingHandCursor)
        btn.clicked.connect(self._redraw_atom)
        self.right_lay.addWidget(btn)

    def _right_atom_info_card(self):
        card = Card("原子-轨道-PDOS 信息")
        card.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        self._atom_info_cells = {}
        table = QGridLayout()
        table.setContentsMargins(0, 0, 0, 0)
        table.setHorizontalSpacing(6)
        table.setVerticalSpacing(0)
        for c, h in enumerate(("曲线", "\u2191 带中心", "\u2193 带中心")):
            lab = QLabel(h)
            lab.setObjectName("CardTitle")
            lab.setAlignment(Qt.AlignCenter)
            lab.setContentsMargins(0, 5, 0, 5)
            table.addWidget(lab, 0, c)
        table.addWidget(self._table_sep(), 1, 0, 1, 3)
        for i, col in enumerate(ATOM_DOS_COLS):
            rr = 2 + 2 * i
            name = QLabel(col)
            name.setObjectName("Value")
            name.setAlignment(Qt.AlignCenter)
            name.setContentsMargins(0, 4, 0, 4)
            table.addWidget(name, rr, 0)
            for c, spin in ((1, "up"), (2, "down")):
                lab = QLabel("—")
                lab.setObjectName("Value")
                lab.setAlignment(Qt.AlignCenter)
                lab.setContentsMargins(0, 4, 0, 4)
                table.addWidget(lab, rr, c)
                self._atom_info_cells[(col, spin)] = lab
            if i < len(ATOM_DOS_COLS) - 1:
                table.addWidget(self._table_sep(), rr + 1, 0, 1, 3)
        for c in range(3):
            table.setColumnStretch(c, 1)
        card.add_layout(table)
        card.body.addStretch(1)
        self.right_lay.addWidget(card, 1)
        self._refresh_atom_info()

    def _refresh_atom_info(self):
        cells = getattr(self, "_atom_info_cells", None)
        if not cells:
            return
        rows = {c: (u, d) for c, u, d in self._atom_band_table()}
        for (col, spin), lab in cells.items():
            u, d = rows.get(col, (None, None))
            val = u if spin == "up" else d
            lab.setText(f"{val:.2f}" if val is not None else "—")

    def _atom_color_button(self, col, color, which="color"):
        btn = QPushButton()
        btn.setObjectName("Secondary")
        btn.setCursor(Qt.PointingHandCursor)
        btn.setFixedHeight(26)
        self._swatch_btn(btn, color)
        btn.clicked.connect(
            lambda _=False, c=col, b=btn, w=which:
            self._pick_atom_color(c, b, w))
        return btn

    def _pick_atom_color(self, col, btn, which="color"):
        line = self._atom_line(col)
        dlg = QColorDialog(QColor(line[which]), self)
        dlg.setWindowTitle("选择颜色")
        dlg.setWindowModality(Qt.ApplicationModal)
        dlg.setWindowFlag(Qt.WindowStaysOnTopHint, True)
        if dlg.exec():
            c = dlg.currentColor()
            if c.isValid():
                line[which] = c.name()
                self._swatch_btn(btn, c.name())
                # 带中心颜色始终与曲线一致
                if which == "color":
                    line["bc_color"] = c.name()

    # ---- 原子-轨道-PDOS：底部曲线样式面板 ----
    def _build_atom_style_panel(self):
        self._clear_layout(self.dos_bottom_lay)

        inner = QWidget()
        grid = QGridLayout(inner)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(4)
        grid.setVerticalSpacing(6)
        # 一行一条曲线：曲线 | 曲线(显示/色/型/宽) | 带中心(开关/型/宽/透/签/号/高)
        heads = ("曲线", "显示", "颜色", "线型", "线宽",
                 "带中心", "线型", "线宽", "透明度",
                 "标签", "字号", "标签高度")
        mins = (42, 30, 30, 42, 34, 34, 42, 34, 36, 30, 34, 44)
        for c, txt in enumerate(heads):
            grid.setColumnStretch(c, 1)
            lab = QLabel(txt)
            lab.setObjectName("FieldLabel")
            lab.setMinimumWidth(mins[c])
            if c == 5:
                lab.setToolTip("勾选后绘制该曲线的带中心参考线")
            if c == 11:
                lab.setToolTip("标签高度：坐标轴比例 0~1（0=底部，1=顶部）")
            grid.addWidget(lab, 0, c)
        grid.addWidget(self._table_sep(), 1, 0, 1, 12)

        self._atom_line_ctrls = {}
        target = self._atom_index()
        entry = self._get_atom_dos(target) if target is not None else {}
        cols = self._atom_keys(entry)
        r = 2
        for i, col in enumerate(cols):
            line = self._atom_line(col)

            name = QLabel(col)
            name.setObjectName("FieldLabel")
            name.setMinimumWidth(mins[0])
            name.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
            name.setToolTip(col)
            grid.addWidget(name, r, 0)

            sw = self._compact_switch(line.get("show", True))
            grid.addWidget(sw, r, 1, Qt.AlignLeft)

            btn_c = self._atom_color_button(col, line["color"])
            btn_c.setMinimumWidth(mins[2])
            btn_c.setMaximumWidth(120)
            btn_c.setFixedHeight(24)
            btn_c.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(btn_c, r, 2)

            cb = ComboBox()
            for lab, val in (("实线", "-"), ("虚线", "--"),
                             ("点线", ":"), ("点划线", "-.")):
                cb.addItem(lab, val)
            cb.setCurrentIndex(max(0, cb.findData(line["ls"])))
            self._slim(cb)
            grid.addWidget(cb, r, 3)

            lw = self._dos_spin(line["lw"], 0.5, 10.0, step=0.5)
            lw.setMinimumWidth(mins[4])
            lw.setMaximumWidth(120)
            lw.setFixedHeight(24)
            lw.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(lw, r, 4)

            bc = self._compact_switch(line.get("bc", False))
            bc.setToolTip("绘制该曲线的带中心参考线")
            grid.addWidget(bc, r, 5, Qt.AlignLeft)

            bcls = ComboBox()
            for lab, val in (("实线", "-"), ("虚线", "--"),
                             ("点线", ":"), ("点划线", "-.")):
                bcls.addItem(lab, val)
            bcls.setCurrentIndex(max(0, bcls.findData(line["bc_ls"])))
            self._slim(bcls)
            grid.addWidget(bcls, r, 6)

            bclw = self._dos_spin(line["bc_lw"], 0.5, 10.0, step=0.5)
            bclw.setMinimumWidth(mins[7])
            bclw.setMaximumWidth(120)
            bclw.setFixedHeight(24)
            bclw.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(bclw, r, 7)

            bcalpha = self._dos_spin(line["bc_alpha"], 0.0, 1.0, step=0.05)
            bcalpha.setMinimumWidth(mins[8])
            bcalpha.setMaximumWidth(120)
            bcalpha.setFixedHeight(24)
            bcalpha.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(bcalpha, r, 8)

            bclabel = self._compact_switch(line.get("bc_label", True))
            grid.addWidget(bclabel, r, 9, Qt.AlignLeft)

            bcfs = self._dos_spin(line["bc_fs"], 6.0, 40.0, step=1.0)
            bcfs.setMinimumWidth(mins[10])
            bcfs.setMaximumWidth(120)
            bcfs.setFixedHeight(24)
            bcfs.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(bcfs, r, 10)

            bcy = self._dos_spin(line.get("bc_y", 0.97), 0.0, 1.0, step=0.02)
            bcy.setButtonSymbols(QAbstractSpinBox.NoButtons)
            bcy.setMinimumWidth(mins[11])
            bcy.setMaximumWidth(120)
            bcy.setFixedHeight(24)
            bcy.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(bcy, r, 11)

            self._atom_line_ctrls[col] = {
                "show": sw, "ls": cb, "lw": lw, "bc": bc,
                "bc_ls": bcls, "bc_lw": bclw, "bc_alpha": bcalpha,
                "bc_label": bclabel, "bc_fs": bcfs, "bc_y": bcy}

            r += 1
            if i < len(cols) - 1:
                grid.addWidget(self._table_sep(), r, 0, 1, 12)
                grid.setRowMinimumHeight(r, 7)
                r += 1

        grid.setRowStretch(r, 1)

        if not self._atom_line_ctrls:
            hint = QLabel("未找到该原子的 PDOS_EIG_*.dat(.bak)")
            hint.setObjectName("Hint")
            hint.setWordWrap(True)
            grid.addWidget(hint, r, 0, 1, 12)

        # 最多显示前 4 行，更多曲线时纵向滚动，不挤占作图高度
        scroll = QScrollArea()
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setWidget(inner)
        if len(cols) > 4:
            grid.activate()
            it = grid.itemAtPosition(2 + 2 * 3, 0)
            limit = 176
            if (it is not None and it.widget() is not None
                    and it.widget().geometry().bottom() > 0):
                limit = it.widget().geometry().bottom() + 4
            scroll.setMaximumHeight(limit)
        self.dos_bottom_lay.addWidget(scroll)

    # ---- 元素-轨道-PDOS：右栏 ----
    def _right_eorbit_pick_card(self):
        card = Card("指定元素")
        grid = self._make_info_grid()
        grid.addWidget(field_label("元素", 76), 0, 0)
        self.cb_eorbit = ComboBox()
        self.cb_eorbit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        for name, val in self._eorbit_options():
            self.cb_eorbit.addItem(name, val)
        st = self._eorbit_style()
        if self._elements:
            idx = self.cb_eorbit.findData(st.get("elem"))
            if idx < 0:
                idx = 0
            self.cb_eorbit.setCurrentIndex(idx)
            st["elem"] = self.cb_eorbit.itemData(idx)
        else:
            self.cb_eorbit.addItem("（未发现元素）", None)
        self.cb_eorbit.currentIndexChanged.connect(self._on_eorbit_changed)
        grid.addWidget(self.cb_eorbit, 0, 1)
        card.add_layout(grid)
        card.add(hint_label("列表来自体系中含有的元素（Data-element）"))
        self.right_lay.addWidget(card)

    def _right_eorbit_card(self):
        st = self._eorbit_style()
        card = Card("曲线参数")
        card.vbox.setContentsMargins(12, 14, 12, 16)
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(9)
        grid.setColumnStretch(1, 1)
        r = 0

        self.sp_eo_xmin = self._dos_spin(st["xmin"], -10000, 10000)
        self.sp_eo_xmax = self._dos_spin(st["xmax"], -10000, 10000)
        grid.addWidget(field_label("横轴区间", 76), r, 0)
        grid.addLayout(self._pair(self.sp_eo_xmin, self.sp_eo_xmax), r, 1)
        r += 1

        self.sp_eo_xtick = self._dos_spin(st["xtick"], 0.0, 5000.0, step=0.5)
        grid.addWidget(field_label("横轴刻度", 76), r, 0)
        grid.addWidget(self.sp_eo_xtick, r, 1)
        r += 1

        self.sp_eo_ymin = self._dos_spin(st["ymin"], -1e5, 1e5, step=5.0)
        self.sp_eo_ymax = self._dos_spin(st["ymax"], -1e5, 1e5, step=5.0)
        grid.addWidget(field_label("纵轴区间", 76), r, 0)
        grid.addLayout(self._pair(self.sp_eo_ymin, self.sp_eo_ymax), r, 1)
        r += 1

        self.sp_eo_ytick = self._dos_spin(st["ytick"], 0.0, 1e5, step=1.0)
        grid.addWidget(field_label("纵轴刻度", 76), r, 0)
        grid.addWidget(self.sp_eo_ytick, r, 1)
        r += 1

        self.sw_eo_fill = self._grid_switch(grid, r, "是否填充",
                                            st.get("fill", True))
        r += 1

        card.add_layout(grid)
        self.right_lay.addWidget(card)

        btn = QPushButton("重新绘制")
        btn.setObjectName("Primary")
        btn.setCursor(Qt.PointingHandCursor)
        btn.clicked.connect(self._redraw_eorbit)
        self.right_lay.addWidget(btn)

    def _right_eorbit_info_card(self):
        card = Card("元素-轨道-PDOS 信息")
        card.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        self._eorbit_info_cells = {}
        table = QGridLayout()
        table.setContentsMargins(0, 0, 0, 0)
        table.setHorizontalSpacing(6)
        table.setVerticalSpacing(0)
        for c, h in enumerate(("曲线", "\u2191 带中心", "\u2193 带中心")):
            lab = QLabel(h)
            lab.setObjectName("CardTitle")
            lab.setAlignment(Qt.AlignCenter)
            lab.setContentsMargins(0, 5, 0, 5)
            table.addWidget(lab, 0, c)
        table.addWidget(self._table_sep(), 1, 0, 1, 3)
        for i, col in enumerate(ATOM_DOS_COLS):
            rr = 2 + 2 * i
            name = QLabel(col)
            name.setObjectName("Value")
            name.setAlignment(Qt.AlignCenter)
            name.setContentsMargins(0, 4, 0, 4)
            table.addWidget(name, rr, 0)
            for c, spin in ((1, "up"), (2, "down")):
                lab = QLabel("—")
                lab.setObjectName("Value")
                lab.setAlignment(Qt.AlignCenter)
                lab.setContentsMargins(0, 4, 0, 4)
                table.addWidget(lab, rr, c)
                self._eorbit_info_cells[(col, spin)] = lab
            if i < len(ATOM_DOS_COLS) - 1:
                table.addWidget(self._table_sep(), rr + 1, 0, 1, 3)
        for c in range(3):
            table.setColumnStretch(c, 1)
        card.add_layout(table)
        card.body.addStretch(1)
        self.right_lay.addWidget(card, 1)
        self._refresh_eorbit_info()

    def _refresh_eorbit_info(self):
        cells = getattr(self, "_eorbit_info_cells", None)
        if not cells:
            return
        rows = {c: (u, d) for c, u, d in self._eorbit_band_table()}
        for (col, spin), lab in cells.items():
            u, d = rows.get(col, (None, None))
            val = u if spin == "up" else d
            lab.setText(f"{val:.2f}" if val is not None else "—")

    def _eorbit_color_button(self, col, color, which="color"):
        btn = QPushButton()
        btn.setObjectName("Secondary")
        btn.setCursor(Qt.PointingHandCursor)
        btn.setFixedHeight(26)
        self._swatch_btn(btn, color)
        btn.clicked.connect(
            lambda _=False, c=col, b=btn, w=which:
            self._pick_eorbit_color(c, b, w))
        return btn

    def _pick_eorbit_color(self, col, btn, which="color"):
        line = self._eorbit_line(col)
        dlg = QColorDialog(QColor(line[which]), self)
        dlg.setWindowTitle("选择颜色")
        dlg.setWindowModality(Qt.ApplicationModal)
        dlg.setWindowFlag(Qt.WindowStaysOnTopHint, True)
        if dlg.exec():
            c = dlg.currentColor()
            if c.isValid():
                line[which] = c.name()
                self._swatch_btn(btn, c.name())
                if which == "color":
                    line["bc_color"] = c.name()

    # ---- 元素-轨道-PDOS：底部曲线样式面板 ----
    def _build_eorbit_style_panel(self):
        self._clear_layout(self.dos_bottom_lay)

        inner = QWidget()
        grid = QGridLayout(inner)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(4)
        grid.setVerticalSpacing(6)
        heads = ("曲线", "显示", "颜色", "线型", "线宽",
                 "带中心", "线型", "线宽", "透明度",
                 "标签", "字号", "标签高度")
        mins = (42, 30, 30, 42, 34, 34, 42, 34, 36, 30, 34, 44)
        for c, txt in enumerate(heads):
            grid.setColumnStretch(c, 1)
            lab = QLabel(txt)
            lab.setObjectName("FieldLabel")
            lab.setMinimumWidth(mins[c])
            if c == 5:
                lab.setToolTip("勾选后绘制该曲线的带中心参考线")
            if c == 11:
                lab.setToolTip("标签高度：坐标轴比例 0~1（0=底部，1=顶部）")
            grid.addWidget(lab, 0, c)
        grid.addWidget(self._table_sep(), 1, 0, 1, 12)

        self._eorbit_line_ctrls = {}
        target = self._eorbit_index()
        entry = self._get_elem_orbit_dos(target) if target is not None else {}
        cols = self._atom_keys(entry)
        r = 2
        for i, col in enumerate(cols):
            line = self._eorbit_line(col)

            name = QLabel(col)
            name.setObjectName("FieldLabel")
            name.setMinimumWidth(mins[0])
            name.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
            name.setToolTip(col)
            grid.addWidget(name, r, 0)

            sw = self._compact_switch(line.get("show", True))
            grid.addWidget(sw, r, 1, Qt.AlignLeft)

            btn_c = self._eorbit_color_button(col, line["color"])
            btn_c.setMinimumWidth(mins[2])
            btn_c.setMaximumWidth(120)
            btn_c.setFixedHeight(24)
            btn_c.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(btn_c, r, 2)

            cb = ComboBox()
            for lab, val in (("实线", "-"), ("虚线", "--"),
                             ("点线", ":"), ("点划线", "-.")):
                cb.addItem(lab, val)
            cb.setCurrentIndex(max(0, cb.findData(line["ls"])))
            self._slim(cb)
            grid.addWidget(cb, r, 3)

            lw = self._dos_spin(line["lw"], 0.5, 10.0, step=0.5)
            lw.setMinimumWidth(mins[4])
            lw.setMaximumWidth(120)
            lw.setFixedHeight(24)
            lw.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(lw, r, 4)

            bc = self._compact_switch(line.get("bc", False))
            bc.setToolTip("绘制该曲线的带中心参考线")
            grid.addWidget(bc, r, 5, Qt.AlignLeft)

            bcls = ComboBox()
            for lab, val in (("实线", "-"), ("虚线", "--"),
                             ("点线", ":"), ("点划线", "-.")):
                bcls.addItem(lab, val)
            bcls.setCurrentIndex(max(0, bcls.findData(line["bc_ls"])))
            self._slim(bcls)
            grid.addWidget(bcls, r, 6)

            bclw = self._dos_spin(line["bc_lw"], 0.5, 10.0, step=0.5)
            bclw.setMinimumWidth(mins[7])
            bclw.setMaximumWidth(120)
            bclw.setFixedHeight(24)
            bclw.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(bclw, r, 7)

            bcalpha = self._dos_spin(line["bc_alpha"], 0.0, 1.0, step=0.05)
            bcalpha.setMinimumWidth(mins[8])
            bcalpha.setMaximumWidth(120)
            bcalpha.setFixedHeight(24)
            bcalpha.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(bcalpha, r, 8)

            bclabel = self._compact_switch(line.get("bc_label", True))
            grid.addWidget(bclabel, r, 9, Qt.AlignLeft)

            bcfs = self._dos_spin(line["bc_fs"], 6.0, 40.0, step=1.0)
            bcfs.setMinimumWidth(mins[10])
            bcfs.setMaximumWidth(120)
            bcfs.setFixedHeight(24)
            bcfs.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(bcfs, r, 10)

            bcy = self._dos_spin(line.get("bc_y", 0.97), 0.0, 1.0, step=0.02)
            bcy.setButtonSymbols(QAbstractSpinBox.NoButtons)
            bcy.setMinimumWidth(mins[11])
            bcy.setMaximumWidth(120)
            bcy.setFixedHeight(24)
            bcy.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            grid.addWidget(bcy, r, 11)

            self._eorbit_line_ctrls[col] = {
                "show": sw, "ls": cb, "lw": lw, "bc": bc,
                "bc_ls": bcls, "bc_lw": bclw, "bc_alpha": bcalpha,
                "bc_label": bclabel, "bc_fs": bcfs, "bc_y": bcy}

            r += 1
            if i < len(cols) - 1:
                grid.addWidget(self._table_sep(), r, 0, 1, 12)
                grid.setRowMinimumHeight(r, 7)
                r += 1

        grid.setRowStretch(r, 1)

        if not self._eorbit_line_ctrls:
            hint = QLabel("未找到该元素的 PDOS_EIG_*.dat(.bak)")
            hint.setObjectName("Hint")
            hint.setWordWrap(True)
            grid.addWidget(hint, r, 0, 1, 12)

        scroll = QScrollArea()
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setWidget(inner)
        if len(cols) > 4:
            grid.activate()
            it = grid.itemAtPosition(2 + 2 * 3, 0)
            limit = 176
            if (it is not None and it.widget() is not None
                    and it.widget().geometry().bottom() > 0):
                limit = it.widget().geometry().bottom() + 4
            scroll.setMaximumHeight(limit)
        self.dos_bottom_lay.addWidget(scroll)

    def _analysis_changed(self):
        """指定原子分析列表变化后同步界面。"""
        self._refresh_analysis_list()
        self._atom_warned = False
        self._atom_band_cache = {}
        if self._current_key == "atom":
            self._rebuild_right("atom")
            self._build_atom_style_panel()
            self._plot_atom()

    # ---- 小工具 ----
    def _make_info_grid(self):
        from PySide6.QtWidgets import QGridLayout
        grid = QGridLayout()
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(8)
        grid.setColumnStretch(1, 1)
        return grid

    def _switch_box(self, text, value, hint_text=None):
        box = QVBoxLayout()
        box.setSpacing(2)
        row = QHBoxLayout()
        row.setSpacing(10)
        lab = QLabel(text)
        lab.setObjectName("FieldLabel")
        row.addWidget(lab)
        row.addStretch(1)
        sw = Switch()
        sw.setChecked(value, animate=False)
        sw.apply_theme(THEME)
        self._switches.append(sw)
        row.addWidget(sw)
        box.addLayout(row)
        if hint_text:
            box.addWidget(hint_label(hint_text))
        return sw, box

    # ---- 参数回调 ----
    def _on_add_atoms(self, payload):
        try:
            items = json.loads(payload)
        except Exception:
            return
        added = 0
        for it in items:
            try:
                idx = int(it.get("index", 0))
            except (TypeError, ValueError):
                continue
            if idx <= 0 or any(a["index"] == idx for a in self.analysis_atoms):
                continue
            frac = it.get("frac") or [0.0, 0.0, 0.0]
            self.analysis_atoms.append({
                "index": idx, "elem": it.get("elem", "?"),
                "frac": [float(x) for x in frac]})
            added += 1
        self._analysis_changed()
        if added:
            self.status.setText(f"已添加 {added} 个原子到指定原子分析列表")

    def _refresh_analysis_list(self):
        lbl = getattr(self, "lbl_analysis_count", None)
        if lbl is not None:
            try:
                lbl.setText(f"{len(self.analysis_atoms)} 个原子")
            except RuntimeError:
                self.lbl_analysis_count = None
        lst = getattr(self, "analysis_list", None)
        if lst is None:
            return
        try:
            lst.clear()
            for item in self.analysis_atoms:
                row = self._analysis_row(item)
                li = QListWidgetItem()
                li.setSizeHint(QSize(0, 36))
                lst.addItem(li)
                lst.setItemWidget(li, row)
        except RuntimeError:
            self.analysis_list = None

    def _remove_analysis_atom(self, idx):
        self.analysis_atoms = [a for a in self.analysis_atoms
                               if a["index"] != idx]
        self._analysis_changed()

    def _clear_analysis(self):
        self.analysis_atoms = []
        self._analysis_changed()

    def _on_curve_param(self, *_):
        if self._current_key in (None, "structure"):
            return
        self._load_dos_page(self._current_key)

    def _on_atom_changed(self, *_):
        data = self.cb_atom.currentData()
        st = self._atom_style()
        st["atom"] = data
        if data is not None:
            self._atom_sel = int(data)
        self._plot_atom()
        self._refresh_atom_info()

    def _on_eorbit_changed(self, *_):
        data = self.cb_eorbit.currentData()
        self._eorbit_style()["elem"] = data
        self._plot_eorbit()
        self._refresh_eorbit_info()

    def _update_paths_label(self):
        p = self.paths
        root = self.data_root

        def nm(path):
            return path.name if path else "（未找到）"

        def brk(s):
            # 在路径分隔符后插入零宽空格，使长路径能在任意分隔处换行
            return str(s).replace("\\", "\\\u200b").replace("/", "/\u200b")

        text = (f"根目录：{brk(root) if root else '（未选择）'}\n"
                f"CONTCAR：{nm(p.get('contcar'))}\n"
                f"Data-total：{nm(p.get('total'))}\n"
                f"Data-element：{nm(p.get('element'))}\n"
                f"Data-atom：{nm(p.get('atom'))}")
        self.lbl_paths.setText(text)

    def _pick_data_root(self):
        start = str(self.data_root or BASE_DIR)
        d = QFileDialog.getExistingDirectory(
            self, "选择 data 所在路径（含 Data-atom / Data-element / Data-total）",
            start)
        if not d:
            return
        self._set_data_root(d)
        self._reload_data()
        name = Path(d).name or d
        missing = [k for k in ("total", "element", "atom")
                   if not self.paths.get(k)]
        if missing:
            self.status.setText(
                f"已选择：{name}（未找到 {'、'.join(missing)} 子目录）")
        else:
            elems = "、".join(self._elements) or "无"
            self.status.setText(f"已自动加载：{name}（元素：{elems}）")

    def _reload_data(self):
        self._data = self._load_data()
        self._atom_cache = {}
        self._atom_band_cache = {}
        self._elem_orbit_cache = {}
        self._eorbit_band_cache = {}
        self._gauss_data = None
        self._gauss_smoothed = None
        self._gauss_path = None
        if getattr(self, "gauss_style", None) is not None:
            self.gauss_style["file"] = None
            self.gauss_style["col"] = None
        self._rebuild_sidebar()
        valid = [k for k, _ in self._nav_items()]
        key = self._current_key if self._current_key in valid else "structure"
        self._select_nav(key)
        self.status.setText("数据已重新加载")

    def _toggle_theme(self):
        self._dark = not self._dark
        apply_theme(QApplication.instance(), DARK if self._dark else LIGHT)
        self.theme_btn.setText("☀️" if self._dark else "🌙")
        for sw in self._switches:
            try:
                sw.apply_theme(THEME)
            except RuntimeError:
                pass
        self.dos_view.set_theme(self._dark)
        if self._current_key == "structure":
            self._load_structure_page()
        elif self._current_key is not None:
            if self._current_key == "element":
                self._build_elem_style_panel()
            elif self._current_key == "atom":
                self._build_atom_style_panel()
            elif self._current_key == "eorbit":
                self._build_eorbit_style_panel()
            elif self._current_key == "gauss":
                self._build_gauss_table_panel()
            self._load_dos_page(self._current_key)
        self.update()


def main():
    # QtWebEngine 需要共享 OpenGL 上下文，否则 WebEngine/WebGL 渲染会不稳定、闪烁
    try:
        QApplication.setAttribute(Qt.AA_ShareOpenGLContexts, True)
    except Exception:
        pass
    # Windows 任务栏图标：必须设置显式 AppUserModelID，否则会显示为 python 的图标
    try:
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(
            "pdos.structure.studio.1")
    except Exception:
        pass
    app = QApplication(sys.argv)
    app.setApplicationName("PDOS")
    app.setFont(QFont("Microsoft YaHei UI", 10))
    apply_theme(app, LIGHT)
    win = MainWindow()
    app.setWindowIcon(win.app_icon)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
