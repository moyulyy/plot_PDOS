#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
使用 Origin 的 Python 库（originpro）绘制 VASP DOS 图像
========================================================

参考 ``test/态密度.opju`` 的样式，把 ``Data-total`` / ``Data-element`` /
``Data-atom`` 中的态密度数据导入 Origin，并生成可直接用于论文的 DOS 图：

    * ``TDOS``                  —— 总态密度（上/下自旋）
    * ``System-DOS``            —— TDOS + 各元素 PDOS 合并图（.opju 参考样式）
    * ``Element-Orbital-PDOS``  —— 每个元素的逐轨道 PDOS
    * ``Atom-Orbital-PDOS``     —— 指定原子的逐轨道 PDOS（可选）

前置条件
--------
1. 本机已安装 **Origin 2022 或更新版本**（脚本通过其 Automation Server 通信）；
2. 运行环境已安装 ``originpro``：``pip install originpro``。

本项目建议使用 conda 环境 ``chem_env``（已含 originpro / numpy / Pillow）：

    conda run -n chem_env python origin_dos_plot.py --root test

用法示例
--------
    # 用 test 目录的数据，按参考样式出图并导出 PNG + .opju
    python origin_dos_plot.py --root test --out test/_origin_out

    # 只画总态密度与合并图，并显示 Origin 界面
    python origin_dos_plot.py --root test --pages tdos,system --show

    # 额外画 1、3 号原子的逐轨道 PDOS
    python origin_dos_plot.py --root test --pages atom --atoms 1,3

参考图（test/态密度.opju）要点
------------------------------
    * X 轴：``Energy (eV)``，范围 -6 ~ 6，步长 1
    * Y 轴：``Density of states``，范围 -80 ~ 80，步长 20
    * 曲线：TDOS(灰) / 各元素(配色)，上自旋为正、下自旋为负
    * 图例：每类一条（TDOS、Al、Ni、O …）
    * 标注：金属元素 d 带中心，如 ``ε_d = -1.30 eV``

注意
----
若本机 Origin 未激活（处于演示模式），导出图片会带有 ``demo`` 水印，
这属于 Origin 授权问题，与本脚本无关；请激活 Origin 后重新导出。
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

# ==========================================================================
# 常量
# ==========================================================================
DEFAULT_LINE_WIDTH = 1.0

# 合并图中 TDOS 与各元素的配色（顺序与参考图一致：TDOS 灰 → 绿 → 蓝 → 粉）
TDOS_COLOR = "#808080"
ELEM_COLORS = ("#00D28C", "#1E3CFF", "#FF0080", "#FFB000", "#00C0FF",
               "#B000FF", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22")

# 逐轨道曲线的颜色（13 列，与 pdos_viewer 中 ATOM_DOS_COLS 顺序一致）
ORBIT_COLORS = ("#1F77B4", "#D62728", "#2CA02C", "#9467BD", "#FF7F0E",
                "#17BECF", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22",
                "#4C72B0", "#DD8452", "#937860")

# Data-atom / Data-element 的列顺序（与文件表头一致）
DOS_COLS = ("s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz",
            "dx2-y2", "tot", "sum_s", "sum_p", "sum_d")

# 能量窗口（与参考图一致）
XMIN, XMAX, XSTEP = -6.0, 6.0, 1.0


# ==========================================================================
# 数据解析（与 pdos_viewer 保持一致的最小实现）
# ==========================================================================
def read_dos_file(path):
    """读取 VASPKIT PDOS/TDOS 文件，返回 (energies, {列名: 数值列表})。"""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        return [], {}
    names = lines[0].split()[1:]
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


def parse_fermi(path):
    """从 FERMI_ENERGY 文件里取第一个浮点数。"""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        for tok in line.split():
            try:
                return float(tok)
            except ValueError:
                continue
    return 0.0


def discover_elements(element_dir):
    """从 Data-element 的 PDOS_EIG_UP_*.dat(.bak) 发现元素列表。"""
    elems = []
    if not element_dir:
        return elems
    for f in sorted(Path(element_dir).glob("PDOS_EIG_UP_*.dat*")):
        name = f.name.split("PDOS_EIG_UP_", 1)[-1]
        for suf in (".dat.bak", ".dat"):
            if name.endswith(suf):
                name = name[: -len(suf)]
                break
        if name and name not in elems:
            elems.append(name)
    return elems


def parse_band_center(path):
    """读取 d 带中心文件，返回 (rows, average)；average=(s,p,d) 或 None。"""
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


def _spin_file(directory, tag, key):
    """优先返回 .dat.bak，否则 .dat。"""
    p = Path(directory) / f"PDOS_EIG_{tag}_{key}.dat.bak"
    if not p.exists():
        p = Path(directory) / f"PDOS_EIG_{tag}_{key}.dat"
    return p if p.exists() else None


def _find_dir(root, name):
    root = Path(root)
    if not root.is_dir():
        return None
    for p in root.iterdir():
        if p.is_dir() and p.name.lower() == name.lower():
            return p
    return None


def load_data(root):
    """读取一个 data 根目录下的全部 DOS 数据。"""
    root = Path(root)
    total_dir = _find_dir(root, "Data-total")
    element_dir = _find_dir(root, "Data-element")
    atom_dir = _find_dir(root, "Data-atom")

    data = {
        "root": root,
        "total_dir": total_dir,
        "element_dir": element_dir,
        "atom_dir": atom_dir,
        "fermi": 0.0,
        "tdos": None,
        "elements": [],
        "elem": {},
        "band": {},
        "total_band": None,
        "atom_ids": [],
    }

    if total_dir:
        fe = total_dir / "FERMI_ENERGY"
        if fe.exists():
            data["fermi"] = parse_fermi(fe)
        tdos = total_dir / "TDOS.dat"
        if tdos.exists():
            data["tdos"] = read_dos_file(tdos)
        bc = total_dir / "BAND_CENTER"
        if bc.exists():
            data["total_band"] = parse_band_center(bc)[1]

    data["elements"] = discover_elements(element_dir)
    if element_dir:
        for elem in data["elements"]:
            entry = {}
            up = _spin_file(element_dir, "UP", elem)
            down = _spin_file(element_dir, "DW", elem)
            if up:
                entry["up"] = read_dos_file(up)
            if down:
                entry["down"] = read_dos_file(down)
            data["elem"][elem] = entry
            bc = element_dir / f"BAND_CENTER_p-dbandcenter_{elem}"
            if bc.exists():
                data["band"][elem] = parse_band_center(bc)[1]

    if atom_dir:
        ids = set()
        for f in Path(atom_dir).glob("PDOS_EIG_UP_*.dat*"):
            name = f.stem.replace("PDOS_EIG_UP_", "")
            if name.isdigit():
                ids.add(int(name))
        data["atom_ids"] = sorted(ids)

    return data


def load_atom(atom_dir, n):
    """读取某个原子的 up/down PDOS，返回 {'up': (e, cols), 'down': (e, cols)}。"""
    entry = {}
    if not atom_dir:
        return entry
    up = _spin_file(atom_dir, "UP", n)
    down = _spin_file(atom_dir, "DW", n)
    if up:
        entry["up"] = read_dos_file(up)
    if down:
        entry["down"] = read_dos_file(down)
    return entry


# ==========================================================================
# Origin 辅助函数
# ==========================================================================
def import_origin(show=False):
    """导入 originpro 并连接 Origin；失败时抛出带提示的异常。"""
    try:
        import originpro as op
    except ImportError as exc:  # pragma: no cover - 取决于运行环境
        raise SystemExit(
            "未检测到 originpro。\n"
            "请先安装 Origin（2022+），再执行：pip install originpro\n"
            f"（原始错误：{exc}）")
    try:
        op.set_show(bool(show))
    except Exception:
        pass
    return op


def _nice_step(span, target_ticks=4):
    """在给定跨度下选一个“好看”的刻度间距（与参考图 step=20 一致）。"""
    if not math.isfinite(span) or span <= 0:
        return 1.0
    raw = span / max(1, target_ticks)
    mag = 10.0 ** math.floor(math.log10(raw))
    for m in (1.0, 2.0, 2.5, 5.0, 10.0):
        if raw <= m * mag:
            return m * mag
    return 10.0 * mag


def _nice_symmetric_limit(values, target_ticks=4):
    """把数据幅值取整为对称坐标轴上限，返回 (limit, tick_step)。"""
    vmax = 0.0
    for v in values:
        try:
            vmax = max(vmax, abs(float(v)))
        except (TypeError, ValueError):
            continue
    if not math.isfinite(vmax) or vmax <= 0:
        return 1.0, 1.0
    target = vmax * 1.10
    step = _nice_step(target, target_ticks)
    limit = math.ceil(target / step) * step
    return limit, step


def _window_values(e, y):
    """取能量窗口 [XMIN, XMAX] 内的 y 值。"""
    if not e or not y:
        return []
    e = np.asarray(e, dtype=float)
    y = np.asarray(y, dtype=float)
    m = (e >= XMIN) & (e <= XMAX)
    return y[m].tolist() if m.any() else y.tolist()


def add_curve(layer, wks, coly, color, *, width=DEFAULT_LINE_WIDTH,
              colx=0, type_="l"):
    """往图层加一条曲线并设置颜色/线宽。"""
    p = layer.add_plot(wks, coly=coly, colx=colx, type=type_)
    if p is None:
        return None
    try:
        p.color = color
    except Exception:
        pass
    if width:
        try:
            p.set_float("width", float(width))
        except Exception:
            pass
    return p


def style_axes(layer, xmin=None, xmax=None, xstep=None,
               ymin=None, ymax=None, ystep=None):
    """设置坐标轴标题与范围。"""
    xmin = XMIN if xmin is None else xmin
    xmax = XMAX if xmax is None else xmax
    xstep = XSTEP if xstep is None else xstep
    layer.axis("x").title = "Energy (eV)"
    layer.axis("y").title = "Density of states"
    layer.set_xlim(xmin, xmax, xstep)
    if ymin is not None and ymax is not None:
        if ystep is None:
            ystep = _nice_step(abs(ymax - ymin))
        layer.set_ylim(ymin, ymax, ystep)


def set_legend(layer, entries):
    """把图例文本设置为 entries（[(显示名, 曲线序号), ...]）。

    关闭自动更新，避免导出重绘时 Origin 把未引用的曲线补进图例。
    """
    if not entries:
        return
    text = "\r\n".join(r"\l(%d) %s" % (idx, name) for name, idx in entries)
    leg = layer.label("legend")
    if leg is None:
        try:
            layer.lt_exec("legend;")
        except Exception:
            pass
        leg = layer.label("legend")
    if leg is None:
        return
    try:
        leg.set_int("auto", 0)
    except Exception:
        pass
    try:
        leg.text = text
    except Exception:
        pass


def annotate_band_center(layer, value, ymax):
    """在图内添加 d 带中心标注：ε_d = -1.30 eV。"""
    if value is None:
        return
    x = XMIN + 0.34 * (XMAX - XMIN)      # 与参考图 (-4.04, 68.3) 接近
    y = 0.85 * ymax
    layer.add_label("\u03b5\\-(d) = %.2f eV" % value, x, y)


def _new_graph(op, lname):
    """新建并命名一个图层（设置 long name 便于重开后定位）。"""
    g = op.new_graph(lname)
    try:
        g.lname = lname
    except Exception:
        pass
    return g


def _find_graph(op, name):
    """按 long name / short name 查找图页；找不到返回 None。"""
    if not name:
        return None
    try:
        g = op.find_graph(name)
        if g is not None:
            return g
    except Exception:
        pass
    i = 0
    while True:
        g = op.find_graph(i)
        if g is None:
            return None
        if g.lname == name or g.name == name:
            return g
        i += 1


# ==========================================================================
# 各类图（只负责“构建 + 返回 (图, 待导出路径)”）
# ==========================================================================
def make_tdos_graph(op, data, out, *, ymin=None, ymax=None):
    """总态密度图（TDOS-UP / TDOS-DOWN）。"""
    if not data["tdos"] or not data["tdos"][0]:
        print("  [跳过] 未找到 TDOS 数据")
        return None
    e, cols = data["tdos"]
    wks = op.new_sheet("w", "TDOS")
    wks.from_list(0, list(e), lname="#Energy", units="eV", axis="X")
    idx = 1
    for key in ("TDOS-UP", "TDOS-DOWN"):
        if key in cols:
            wks.from_list(idx, list(cols[key]), lname=key, units="a.u.",
                          axis="Y")
            idx += 1
    g = _new_graph(op, "TDOS")
    layer = g[0]
    idx = 1
    for key in ("TDOS-UP", "TDOS-DOWN"):
        if key in cols:
            add_curve(layer, wks, idx, TDOS_COLOR)
            idx += 1
    ystep = None
    if ymin is None or ymax is None:
        vals = _window_values(e, cols.get("TDOS-UP", []))
        vals += [-v for v in _window_values(e, cols.get("TDOS-DOWN", []))]
        lim, ystep = _nice_symmetric_limit(vals)
        ymin, ymax = -lim, lim
    style_axes(layer, ymin=ymin, ymax=ymax, ystep=ystep)
    set_legend(layer, [("TDOS", 1)])
    return g, Path(out) / "TDOS.png"


def make_system_graph(op, data, out, *, ymin=None, ymax=None):
    """TDOS + 各元素 PDOS 合并图（参考 test/态密度.opju 的样式）。"""
    elements = [el for el in data["elements"] if data["elem"].get(el)]
    if not data["tdos"] and not elements:
        print("  [跳过] 未找到可直接绘制的 TDOS / 元素数据")
        return None

    if data["tdos"] and data["tdos"][0]:
        energy = list(data["tdos"][0])
    else:
        first = data["elem"][elements[0]]
        energy = list((first.get("up") or first.get("down"))[0])

    wks = op.new_sheet("w", "System-DOS")
    wks.from_list(0, energy, lname="#Energy", units="eV", axis="X")
    col = 1
    plotted = []               # (col_index, color, ydata, group_name)
    if data["tdos"] and "TDOS-UP" in data["tdos"][1]:
        wks.from_list(col, list(data["tdos"][1]["TDOS-UP"]),
                      lname="TDOS-UP", units="a.u.", axis="Y")
        plotted.append((col, TDOS_COLOR, data["tdos"][1]["TDOS-UP"], "TDOS"))
        col += 1
    if data["tdos"] and "TDOS-DOWN" in data["tdos"][1]:
        wks.from_list(col, list(data["tdos"][1]["TDOS-DOWN"]),
                      lname="TDOS-DOWN", units="a.u.", axis="Y")
        plotted.append((col, TDOS_COLOR, data["tdos"][1]["TDOS-DOWN"], "TDOS"))
        col += 1
    for i, elem in enumerate(elements):
        color = ELEM_COLORS[i % len(ELEM_COLORS)]
        for spin, tag in (("up", "up"), ("down", "dn")):
            pair = data["elem"][elem].get(spin)
            if not pair or not pair[0] or "tot" not in pair[1]:
                continue
            wks.from_list(col, list(pair[1]["tot"]),
                          lname=f"{elem}-{tag}", units="a.u.", axis="Y")
            plotted.append((col, color, pair[1]["tot"], elem))
            col += 1

    g = _new_graph(op, "System-DOS")
    layer = g[0]
    for cindex, color, _ydata, _name in plotted:
        add_curve(layer, wks, cindex, color)

    # 每类（TDOS / 元素）只在图例中出现一次，引用其首条曲线
    legend_entries = []
    seen = set()
    for i, (_c, _color, _y, name) in enumerate(plotted):
        if name not in seen:
            seen.add(name)
            legend_entries.append((name, i + 1))

    ystep = None
    if ymin is None or ymax is None:
        vals = []
        for _c, _color, ydata, _name in plotted:
            vals.extend(_window_values(energy, ydata))
        lim, ystep = _nice_symmetric_limit(vals)
        ymin, ymax = -lim, lim
    style_axes(layer, ymin=ymin, ymax=ymax, ystep=ystep)
    set_legend(layer, legend_entries)
    annotate_band_center(layer, _system_band_center(data), ymax)
    return g, Path(out) / "System-DOS.png"


def _system_band_center(data):
    """选择用于标注的 d 带中心（|d| 最大的元素，否则用总态密度）。"""
    best = None
    best_abs = 0.0
    for _elem, avg in data["band"].items():
        if avg and avg[2] is not None and abs(avg[2]) > best_abs:
            best_abs = abs(avg[2])
            best = avg[2]
    if best is not None:
        return best
    if data.get("total_band") and data["total_band"][2] is not None:
        return data["total_band"][2]
    return None


def _orbital_graph(op, entry, lname, out, png_name, *, ymin, ymax):
    """逐轨道 PDOS 图（entry 为 {'up': (e, cols), 'down': (e, cols)}）。"""
    if not entry:
        return None
    up = entry.get("up")
    down = entry.get("down")
    energy = list((up or down)[0])
    wks = op.new_sheet("w", lname)
    wks.from_list(0, energy, lname="#Energy", units="eV", axis="X")
    col = 1
    handles = []               # (col, color, ydata, name)
    for ci, name in enumerate(DOS_COLS):
        color = ORBIT_COLORS[ci % len(ORBIT_COLORS)]
        for spin in ("up", "down"):
            pair = entry.get(spin)
            if not pair or name not in pair[1]:
                continue
            wks.from_list(col, list(pair[1][name]),
                          lname=f"{name}-{'up' if spin == 'up' else 'dn'}",
                          units="a.u.", axis="Y")
            handles.append((col, color, pair[1][name], name))
            col += 1
    g = _new_graph(op, lname)
    layer = g[0]
    for cindex, color, _y, _name in handles:
        add_curve(layer, wks, cindex, color)
    # 每个轨道只显示一条图例（引用其 up 曲线）
    legend_entries = []
    seen = set()
    for i, (_c, _color, _y, name) in enumerate(handles):
        if name not in seen:
            seen.add(name)
            legend_entries.append((name, i + 1))
    set_legend(layer, legend_entries)
    ystep = None
    if ymin is None or ymax is None:
        vals = []
        for _c, _color, y, _name in handles:
            vals.extend(_window_values(energy, y))
        lim, ystep = _nice_symmetric_limit(vals)
        ymin, ymax = -lim, lim
    style_axes(layer, ymin=ymin, ymax=ymax, ystep=ystep)
    return g, Path(out) / png_name


def make_eorbit_graph(op, data, elem, out, *, ymin=None, ymax=None):
    """某个元素的逐轨道 PDOS 图。"""
    entry = data["elem"].get(elem, {})
    if not entry:
        return None
    return _orbital_graph(op, entry, f"{elem}-Orbital-PDOS", out,
                          f"Element-Orbital-PDOS_{elem}.png",
                          ymin=ymin, ymax=ymax)


def make_atom_graph(op, data, n, out, *, ymin=None, ymax=None):
    """某个原子的逐轨道 PDOS 图。"""
    entry = load_atom(data["atom_dir"], n)
    if not entry:
        print(f"  [跳过] 原子 {n} 无 PDOS 文件")
        return None
    return _orbital_graph(op, entry, f"Atom{n}-Orbital-PDOS", out,
                          f"Atom-Orbital-PDOS_{n}.png",
                          ymin=ymin, ymax=ymax)


# ==========================================================================
# 导出
# ==========================================================================
def export_graphs(op, tasks, width=1200):
    """导出 (图名, 路径) 列表；缺图时自动重试。"""
    for name, path in tasks:
        g = _find_graph(op, name)
        if g is None:
            print(f"  [警告] 找不到图 {name}，跳过导出")
            continue
        _save_graph(g, path, width)


def _save_graph(g, path, width=1200):
    try:
        path = Path(path).resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        try:
            g.activate()
        except Exception:
            pass
        res = g.save_fig(str(path), type="png", width=int(width))
        if res and path.exists():
            print(f"  已导出图片：{path}")
            return
        # 退回 LabTalk 直接导出（部分 Origin 版本 save_fig 不返回值）
        g.lt_exec(
            f'expgraph type:=png path:="{path.parent}" filename:="{path.stem}" '
            f'overwrite:=replace tr1.Unit:=2 tr1.Width:={int(width)}')
        if path.exists():
            print(f"  已导出图片：{path}")
        else:
            print(f"  [警告] 图片导出失败：{path}")
    except Exception as exc:
        print(f"  [警告] 图片导出失败 {path}: {exc}")


# ==========================================================================
# 主流程
# ==========================================================================
def _parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="使用 Origin(originpro) 绘制 VASP DOS 图像（参考 test/态密度.opju）")
    ap.add_argument("--root", default="test",
                    help="data 根目录（含 Data-total/Data-element/Data-atom），默认 test")
    ap.add_argument("--out", default=None,
                    help="输出目录（PNG 与 .opju），默认 <root>/_origin_out")
    ap.add_argument("--pages", default="tdos,system,eorbit",
                    help="要生成的图，逗号分隔：tdos,system,eorbit,atom（默认前三个）")
    ap.add_argument("--atoms", default="",
                    help="atom 页要绘制的原子编号，逗号分隔（如 1,3,5）")
    ap.add_argument("--xmin", type=float, default=XMIN, help="X 轴下限")
    ap.add_argument("--xmax", type=float, default=XMAX, help="X 轴上限")
    ap.add_argument("--ymin", type=float, default=None, help="Y 轴下限（默认自动）")
    ap.add_argument("--ymax", type=float, default=None, help="Y 轴上限（默认自动）")
    ap.add_argument("--width", type=int, default=1200, help="导出图片宽度（像素）")
    ap.add_argument("--show", action="store_true", help="显示 Origin 界面")
    ap.add_argument("--project", default=None,
                    help="保存的 .opju 路径，默认 <out>/DOS.opju")
    return ap.parse_args(argv)


def plot_origin(op, data, out, pages, atoms, *, xmin, xmax, ymin, ymax,
                width, project):
    """在 Origin 中生成各类图，保存工程后重开再导出图片。

    说明：Origin 在程序化修改图例后，首次导出可能仍使用旧的自动图例；
    先保存并重新打开工程可确保图例等格式正确生效。
    """
    global XMIN, XMAX
    XMIN, XMAX = xmin, xmax
    op.new(False)

    tasks = []          # (图 long name, png 路径)

    def _collect(result):
        if not result:
            return None
        g, png = result
        name = g.lname or g.name
        tasks.append((name, png))
        return g

    if "tdos" in pages:
        print("[1] 总态密度 TDOS")
        _collect(make_tdos_graph(op, data, out, ymin=ymin, ymax=ymax))
    if "system" in pages:
        print("[2] TDOS + 元素 PDOS 合并图（参考样式）")
        _collect(make_system_graph(op, data, out, ymin=ymin, ymax=ymax))
    if "eorbit" in pages:
        for elem in data["elements"]:
            print(f"[3] 元素逐轨道 PDOS：{elem}")
            _collect(make_eorbit_graph(op, data, elem, out,
                                       ymin=ymin, ymax=ymax))
    if "atom" in pages:
        ids = [int(x) for x in re.split(r"[,，\s]+", atoms) if x.strip().isdigit()]
        if not ids:
            ids = data["atom_ids"][:3]        # 默认前 3 个，避免批量出图过慢
            print(f"[4] 原子逐轨道 PDOS（默认取前 3 个：{ids}）")
        else:
            print(f"[4] 原子逐轨道 PDOS：{ids}")
        for n in ids:
            _collect(make_atom_graph(op, data, n, out,
                                     ymin=ymin, ymax=ymax))

    if not tasks:
        print("没有生成任何图。")
        return []

    # 保存工程并重新打开，确保图例等格式在导出时生效
    if not project:
        project = Path(out) / "DOS.opju"
    project = Path(project).resolve()
    project.parent.mkdir(parents=True, exist_ok=True)
    try:
        if op.save(str(project)):
            print(f"已保存 Origin 工程：{project}")
            op.open(str(project))
    except Exception as exc:
        print(f"[警告] 保存/重开工程失败（将直接导出）：{exc}")

    export_graphs(op, tasks, width=width)
    return tasks


def main(argv=None):
    args = _parse_args(argv)
    root = Path(args.root)
    if not root.exists():
        raise SystemExit(f"data 根目录不存在：{root}")
    out = Path(args.out).resolve() if args.out else (root / "_origin_out").resolve()
    project = Path(args.project).resolve() if args.project else out / "DOS.opju"

    print(f"读取数据：{root}")
    data = load_data(root)
    print(f"  元素：{data['elements'] or '（无）'}"
          f"，TDOS：{'有' if data['tdos'] else '无'}"
          f"，原子：{len(data['atom_ids'])} 个"
          f"，费米能级：{data['fermi']:.4f} eV")

    pages = [p.strip().lower() for p in args.pages.split(",") if p.strip()]
    valid = {"tdos", "system", "eorbit", "atom"}
    unknown = [p for p in pages if p not in valid]
    if unknown:
        raise SystemExit(f"不支持的页面：{unknown}（可选 {sorted(valid)}）")

    op = import_origin(args.show)
    print("已连接 Origin，开始绘图 ……")
    tasks = plot_origin(op, data, out, pages, args.atoms,
                        xmin=args.xmin, xmax=args.xmax,
                        ymin=args.ymin, ymax=args.ymax,
                        width=args.width, project=project)
    print(f"完成：共生成 {len(tasks)} 张图，输出目录 {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
