#!/usr/bin/env bash

# ============================================================
# 这个脚本是专门给 VASPKIT 1.5.1 使用的 DOS / d 带中心后处理脚本。
#
# 文件名：cdos_1.5.1.sh
# 目标版本：/opt/vaspkit.1.5.1/bin/vaspkit
#
# 设计目标：
# 1. 不追求短，追求直白。
# 2. 尽量按"从上到下串行执行"的方式写。
# 3. 多写注释，多打印当前进度，方便你检查每一步在做什么。
# 4. 尽量使用简单变量、简单循环、简单判断。
# 5. 这个脚本专门适配 VASPKIT 1.5.1 的 d 带中心输入格式。
#
# 修改记录（2026-08-17）：
# 1. Data-total / Data-element / Data-atom 目录中不保留 cmd_*.log 日志。
#    - 所有 VASPKIT 调用日志改为写到【当前工作目录】的临时文件，
#      成功完成后立即删除；出错时才保留（也在当前目录），方便排错。
#    - 步骤 12 还会做一次"兜底清理"，把三个输出目录里任何残留的
#      cmd_*.log 全部删除，确保最终输出目录干净。
# 2. Data-atom 中的每个原子 PDOS，与 Data-element 一样：
#    生成 .bak 副本文件，并在副本末尾追加 sum_s / sum_p / sum_d 三列。
# 3. 增加稳健性：
#    - set -euo pipefail：管道中 VASPKIT 失败也能被检测到，未定义变量
#      直接报错，任何命令失败立即终止。
#    - export LC_ALL=C：避免某些系统 locale 下 awk 数字格式异常
#      （例如小数点是逗号的地区），保证输出始终用小数点。
#    - 增加 ERR trap：脚本异常退出时打印出错行号，方便定位。
#    - 步骤 1 增加 awk / grep / cp / mv / rm 等基础命令的存在性检查。
#    - 增加 SIGMA、ENERGY_WINDOW 参数合法性检查。
#    - 增加输出目录名合法性检查，防止误删目录。
#    - 把"追加 sum 三列"的重复代码抽成同一个函数 append_sum_columns，
#      元素和原子共用一份逻辑，避免三份代码不一致。
#    - sum 列的 awk 处理增加字段数判断：遇到字段不足的异常行直接原样
#      输出，不计算错误的求和；同时自动去除行尾 \r（兼容 Windows 文件）。
#
# 这个脚本会做的事情：
# 1. 检查 POSCAR、INCAR、VASPKIT 是否存在。
# 2. 读取 POSCAR，显示元素种类、每种元素原子数、总原子数。
# 3. 读取 INCAR，判断 ISPIN 是 1 还是 2。
# 4. 计算每个元素的 PDOS，保存到 Data-element。
# 5. 计算每个元素的 d 带中心，保存到 Data-element。
# 6. 根据变量选择，决定是否计算每个原子的 PDOS。
# 7. 计算总 TDOS，保存到 Data-total。
# 8. 计算总 d 带中心，保存到 Data-total。
#
# 特别说明：
# VASPKIT 1.5.1 的 d 带中心输入格式，与 1.3.x 不同。
#
# 对 1.3.x，常见正确输入更像是：
# echo -e "503\nN\n1\nS\n" | /opt/vaspkit.1.3.5/bin/vaspkit
#
# 但对 1.5.1，这里固定使用：
# echo -e "503\n1\nN\n1\nS\n" | /opt/vaspkit.1.5.1/bin/vaspkit
#
# 也就是说，1.5.1 这里比 1.3.x 多了一组"\n1"。
# 这个脚本已经按 1.5.1 的方式写死，不再混用其它版本的输入格式。
# ============================================================

set -euo pipefail
export LC_ALL=C

# ------------------------------------------------------------
# 异常兜底：脚本异常终止时，打印出错位置所在的行号。
# 注意：被 if 条件包裹的命令失败不会走到这里（那些都有专门处理）。
# ------------------------------------------------------------
trap 'echo "[ERROR] 脚本在第 $LINENO 行附近出错，脚本异常终止。" >&2' ERR

# ============================================================
# 第一部分：你平时最可能修改的变量区
# ============================================================

# DOS 展宽参数。
SIGMA="0.03"

# 能量窗口。
# 含义：起点 终点 点数
ENERGY_WINDOW="-20.0 20.0 2000"

# 输出目录名称。
DATA_ELEMENT_DIR="Data-element"
DATA_ATOM_DIR="Data-atom"
DATA_TOTAL_DIR="Data-total"

# 这个脚本固定使用 VASPKIT 1.5.1。
VASPKIT_BIN="/opt/vaspkit.1.5.1/bin/vaspkit"

# 是否计算每个原子的 PDOS。
# yes = 计算
# no  = 跳过
RUN_ATOM_PDOS="yes"

# ============================================================
# 第二部分：VASPKIT 1.5.1 的 d 带中心输入参数
#
# 对于 VASPKIT 1.5.1：
# 分元素 d 带中心命令应为：
# echo -e "503\n1\nN\n1\nS\n" | /opt/vaspkit.1.5.1/bin/vaspkit
#
# 总 d 带中心命令应为：
# echo -e "503\n1\nN\n1\nall\n" | /opt/vaspkit.1.5.1/bin/vaspkit
#
# 所以这里一共是 5 行输入。
# ============================================================

DBAND_TASK_CODE="503"
DBAND_FIRST_INPUT="1"
DBAND_CHANGE_WINDOW="N"
DBAND_SELECT_MODE="1"

# ============================================================
# 公共函数：append_sum_columns
#
# 作用：把某个 PDOS 文件复制成 .bak 副本，
#       并在副本末尾追加 3 列：
#         sum_s = s
#         sum_p = py + pz + px
#         sum_d = dxy + dyz + dz2 + dxz + dx2
# 参数：$1 = 源 PDOS 文件路径
# 说明：源文件保持不变，只修改 .bak 副本。
#       该函数同时用于 Data-element 和 Data-atom。
# ============================================================
append_sum_columns() {
    local SRC_FILE="$1"
    local BAK_FILE="${SRC_FILE}.bak"
    local TMP_FILE="${BAK_FILE}.tmp"

    if [ ! -f "$SRC_FILE" ]; then
        echo "[ERROR] [后处理] 源文件不存在，无法追加 sum 列：$SRC_FILE"
        return 1
    fi

    echo "[后处理] 源文件：$SRC_FILE"
    echo "[后处理] 第 1 步：创建副本文件：$BAK_FILE"
    cp "$SRC_FILE" "$BAK_FILE"

    echo "[后处理] 第 2 步：在副本文件中追加 sum_s / sum_p / sum_d 三列"
    awk '
    {
        gsub(/\r/, "", $0)
        if (NR == 1) { print $0 " sum_s sum_p sum_d"; next }
        if (NR == 2) { print $0; next }
        if (NF == 0) { print; next }
        if (NF < 10) { print; next }
        SUM_S = $2
        SUM_P = $3 + $4 + $5
        SUM_D = $6 + $7 + $8 + $9 + $10
        printf "%s %12.5f %12.5f %12.5f\n", $0, SUM_S, SUM_P, SUM_D
    }
    ' "$SRC_FILE" > "$TMP_FILE"

    if [ ! -s "$TMP_FILE" ]; then
        echo "[ERROR] [后处理] awk 处理失败，未生成有效内容：$TMP_FILE"
        rm -f "$TMP_FILE"
        return 1
    fi

    mv "$TMP_FILE" "$BAK_FILE"
    echo "[后处理] 已完成：$BAK_FILE"
    echo "[后处理] 该副本现在是 14 列：原 11 列 + sum_s + sum_p + sum_d"
    return 0
}

# ============================================================
# 第三部分：脚本正式开始
# ============================================================

echo "============================================================"
echo "脚本启动：cdos_1.5.1.sh"
echo "说明：这是一个专门给 VASPKIT 1.5.1 使用的简单直白版脚本。"
echo "说明：本脚本会严格使用 1.5.1 对应的 d 带中心输入格式。"
echo "说明：输出目录（Data-element / Data-atom / Data-total）不保留 cmd_*.log。"
echo "============================================================"
echo

# ------------------------------------------------------------
# 步骤 1：检查必要文件和程序
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "步骤 1：检查必要文件和程序"
echo "------------------------------------------------------------"

echo "[步骤 1] 检查 POSCAR 是否存在"
if [ ! -f "POSCAR" ]; then
    echo "[ERROR] 当前目录缺少 POSCAR，脚本终止。"
    exit 1
fi

echo "[步骤 1] 检查 INCAR 是否存在"
if [ ! -f "INCAR" ]; then
    echo "[ERROR] 当前目录缺少 INCAR，脚本终止。"
    exit 1
fi

echo "[步骤 1] 检查 VASPKIT 1.5.1 是否存在：$VASPKIT_BIN"
if [ ! -x "$VASPKIT_BIN" ]; then
    echo "[ERROR] 指定的 VASPKIT 1.5.1 不存在，或者没有执行权限。"
    echo "[ERROR] 路径：$VASPKIT_BIN"
    exit 1
fi

echo "[步骤 1] 检查脚本依赖的基础命令是否存在"
for CMD_NAME in awk grep cp mv rm; do
    if ! command -v "$CMD_NAME" > /dev/null 2>&1; then
        echo "[ERROR] 系统缺少命令：$CMD_NAME"
        echo "[ERROR] 请先安装或确认 PATH 中包含该命令。"
        exit 1
    fi
done
echo "[步骤 1] 基础命令检查通过：awk grep cp mv rm"

echo "[步骤 1] 必要文件和程序检查通过。"
echo

# ------------------------------------------------------------
# 步骤 2：规范化 RUN_ATOM_PDOS 的值
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "步骤 2：检查是否计算每个原子的 PDOS"
echo "------------------------------------------------------------"

if [ "$RUN_ATOM_PDOS" = "yes" ] || [ "$RUN_ATOM_PDOS" = "YES" ] || [ "$RUN_ATOM_PDOS" = "y" ] || [ "$RUN_ATOM_PDOS" = "Y" ]; then
    RUN_ATOM_PDOS="yes"
elif [ "$RUN_ATOM_PDOS" = "no" ] || [ "$RUN_ATOM_PDOS" = "NO" ] || [ "$RUN_ATOM_PDOS" = "n" ] || [ "$RUN_ATOM_PDOS" = "N" ]; then
    RUN_ATOM_PDOS="no"
else
    echo "[ERROR] RUN_ATOM_PDOS 的值无效：$RUN_ATOM_PDOS"
    echo "[ERROR] 只能填写 yes/no 或 y/n。"
    exit 1
fi

echo "[步骤 2] RUN_ATOM_PDOS = $RUN_ATOM_PDOS"
echo

# ------------------------------------------------------------
# 步骤 3：读取 INCAR 中的 ISPIN
# 如果 INCAR 没有显式写 ISPIN，则按默认 1 处理。
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "步骤 3：读取 INCAR 中的 ISPIN"
echo "------------------------------------------------------------"

ISPIN_RAW=$(awk '
{
    line = $0
    gsub(/\r/, "", line)
    sub(/#.*/, "", line)
    upper = toupper(line)

    if (upper ~ /^[[:space:]]*ISPIN[[:space:]]*=/) {
        sub(/^[[:space:]]*ISPIN[[:space:]]*=[[:space:]]*/, "", upper)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", upper)
        split(upper, arr, /[[:space:]]+/)
        print arr[1]
        exit
    }
}
' INCAR)

if [ -z "$ISPIN_RAW" ]; then
    ISPIN="1"
    ISPIN_SOURCE="INCAR 没有显式设置 ISPIN，所以按默认值 1 处理"
else
    if [ "$ISPIN_RAW" = "1" ] || [ "$ISPIN_RAW" = "2" ]; then
        ISPIN="$ISPIN_RAW"
        ISPIN_SOURCE="INCAR 显式设置了 ISPIN=$ISPIN_RAW"
    else
        echo "[ERROR] INCAR 中 ISPIN 的值不合法：$ISPIN_RAW"
        echo "[ERROR] 这里只接受 1 或 2。"
        exit 1
    fi
fi

if [ "$ISPIN" = "2" ]; then
    IS_MAGNETIC="yes"
else
    IS_MAGNETIC="no"
fi

echo "[步骤 3] ISPIN = $ISPIN"
echo "[步骤 3] ISPIN 来源：$ISPIN_SOURCE"
if [ "$IS_MAGNETIC" = "yes" ]; then
    echo "[步骤 3] 当前体系按磁性 / 自旋极化体系处理。"
else
    echo "[步骤 3] 当前体系按非磁性 / 非自旋极化体系处理。"
fi
echo

# ------------------------------------------------------------
# 步骤 4：读取 POSCAR 中的元素信息和原子数信息
# POSCAR 第 6 行：元素名
# POSCAR 第 7 行：每种元素的原子数
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "步骤 4：读取 POSCAR 中的结构信息"
echo "------------------------------------------------------------"

ELEMENT_LINE=$(awk 'NR==6 { gsub(/\r/, ""); print; exit }' POSCAR)
COUNT_LINE=$(awk 'NR==7 { gsub(/\r/, ""); print; exit }' POSCAR)

if [ -z "$ELEMENT_LINE" ]; then
    echo "[ERROR] POSCAR 第 6 行为空，无法读取元素信息。"
    exit 1
fi

if [ -z "$COUNT_LINE" ]; then
    echo "[ERROR] POSCAR 第 7 行为空，无法读取原子数信息。"
    exit 1
fi

ELEMENTS=($ELEMENT_LINE)
COUNTS=($COUNT_LINE)

ELEMENT_TYPE_COUNT=${#ELEMENTS[@]}
COUNT_TYPE_COUNT=${#COUNTS[@]}

if [ "$ELEMENT_TYPE_COUNT" -eq 0 ]; then
    echo "[ERROR] 没有从 POSCAR 第 6 行读取到元素。"
    exit 1
fi

if [ "$COUNT_TYPE_COUNT" -eq 0 ]; then
    echo "[ERROR] 没有从 POSCAR 第 7 行读取到原子数。"
    exit 1
fi

if [ "$ELEMENT_TYPE_COUNT" -ne "$COUNT_TYPE_COUNT" ]; then
    echo "[ERROR] POSCAR 第 6 行和第 7 行的项目数不一致。"
    echo "[ERROR] 元素个数 = $ELEMENT_TYPE_COUNT"
    echo "[ERROR] 原子数项目个数 = $COUNT_TYPE_COUNT"
    exit 1
fi

TOTAL_NAT=0
FORMULA_SUMMARY=""

INDEX=0
while [ "$INDEX" -lt "$ELEMENT_TYPE_COUNT" ]; do
    CURRENT_ELEMENT="${ELEMENTS[$INDEX]}"
    CURRENT_COUNT="${COUNTS[$INDEX]}"

    if ! echo "$CURRENT_COUNT" | grep -Eq '^[0-9]+$'; then
        echo "[ERROR] POSCAR 第 7 行中存在非法原子数：$CURRENT_COUNT"
        exit 1
    fi

    TOTAL_NAT=$((TOTAL_NAT + CURRENT_COUNT))
    FORMULA_SUMMARY="${FORMULA_SUMMARY}${CURRENT_ELEMENT}(${CURRENT_COUNT}) "

    INDEX=$((INDEX + 1))
done

if [ "$TOTAL_NAT" -eq 0 ]; then
    echo "[ERROR] POSCAR 中的原子总数是 0，说明第 6/7 行读取有问题。"
    exit 1
fi

echo "[步骤 4] 结构信息读取完成。"
echo "[步骤 4] 元素种类：${ELEMENTS[*]}"
echo "[步骤 4] 元素组成：$FORMULA_SUMMARY"

INDEX=0
while [ "$INDEX" -lt "$ELEMENT_TYPE_COUNT" ]; do
    CURRENT_ELEMENT="${ELEMENTS[$INDEX]}"
    CURRENT_COUNT="${COUNTS[$INDEX]}"
    echo "[步骤 4] 元素 $CURRENT_ELEMENT 的原子数 = $CURRENT_COUNT"
    INDEX=$((INDEX + 1))
done

echo "[步骤 4] 原子总数 = $TOTAL_NAT"
echo

# ------------------------------------------------------------
# 步骤 5：打印本次运行总览，并校验 SIGMA / ENERGY_WINDOW
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "步骤 5：打印本次运行配置"
echo "------------------------------------------------------------"

# 校验 SIGMA 是否为合法数字
if ! echo "$SIGMA" | grep -Eq '^[+-]?[0-9]+([.][0-9]+)?$'; then
    echo "[ERROR] SIGMA 不是合法数字：$SIGMA"
    echo "[ERROR] 请检查脚本顶部 SIGMA 变量的值。"
    exit 1
fi

# 校验 ENERGY_WINDOW 是否正好是 3 个数字（起点 终点 点数）
EW_COUNT=0
for EW_ITEM in $ENERGY_WINDOW; do
    EW_COUNT=$((EW_COUNT + 1))
    if ! echo "$EW_ITEM" | grep -Eq '^[+-]?[0-9]+([.][0-9]+)?$'; then
        echo "[ERROR] ENERGY_WINDOW 中含有非数字项：$EW_ITEM"
        echo "[ERROR] 请检查脚本顶部 ENERGY_WINDOW 变量的值。"
        exit 1
    fi
done

if [ "$EW_COUNT" -ne 3 ]; then
    echo "[ERROR] ENERGY_WINDOW 必须是 3 个数字（起点 终点 点数）。"
    echo "[ERROR] 当前读到的项数 = $EW_COUNT，请检查脚本顶部配置。"
    exit 1
fi

echo "[步骤 5] SIGMA = $SIGMA eV"
echo "[步骤 5] ENERGY_WINDOW = $ENERGY_WINDOW"
echo "[步骤 5] VASPKIT 路径 = $VASPKIT_BIN"
echo "[步骤 5] ISPIN = $ISPIN"
echo "[步骤 5] RUN_ATOM_PDOS = $RUN_ATOM_PDOS"
echo "[步骤 5] VASPKIT 1.5.1 的 d 带中心输入格式 = 503 -> 1 -> N -> 1 -> 目标"
echo

# ------------------------------------------------------------
# 步骤 6：准备输出目录
# 注意：这一步会删除旧的 Data-element / Data-atom / Data-total
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "步骤 6：准备输出目录"
echo "------------------------------------------------------------"

# 防止输出目录名被配错成空或根目录，避免误删
for DIR_NAME in "$DATA_ELEMENT_DIR" "$DATA_ATOM_DIR" "$DATA_TOTAL_DIR"; do
    if [ -z "$DIR_NAME" ] || [ "$DIR_NAME" = "/" ]; then
        echo "[ERROR] 输出目录名配置不合法：$DIR_NAME"
        echo "[ERROR] 请检查脚本顶部 DATA_ELEMENT_DIR / DATA_ATOM_DIR / DATA_TOTAL_DIR。"
        exit 1
    fi
done

echo "[步骤 6] 删除旧输出目录：$DATA_ELEMENT_DIR $DATA_ATOM_DIR $DATA_TOTAL_DIR"
rm -rf "$DATA_ELEMENT_DIR" "$DATA_ATOM_DIR" "$DATA_TOTAL_DIR"

echo "[步骤 6] 重新创建输出目录"
mkdir -p "$DATA_ELEMENT_DIR" "$DATA_ATOM_DIR" "$DATA_TOTAL_DIR"

for DIR_NAME in "$DATA_ELEMENT_DIR" "$DATA_ATOM_DIR" "$DATA_TOTAL_DIR"; do
    if [ ! -d "$DIR_NAME" ]; then
        echo "[ERROR] 输出目录创建失败：$DIR_NAME"
        exit 1
    fi
done

echo "[步骤 6] 输出目录准备完成。"
echo

# ------------------------------------------------------------
# 步骤 7：计算每个元素的 PDOS
# 输出放到 Data-element
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "步骤 7：开始计算每个元素的 PDOS"
echo "------------------------------------------------------------"

INDEX=0
while [ "$INDEX" -lt "$ELEMENT_TYPE_COUNT" ]; do
    CURRENT_ELEMENT="${ELEMENTS[$INDEX]}"
    # 注意：日志写到【当前目录】的临时文件，成功完成后删除，
    # 因此 Data-element 目录中不会保留 cmd_*.log。
    ELEMENT_PDOS_LOG="cmd_PDOS_${CURRENT_ELEMENT}.log"

    echo
    echo "[步骤 7] 当前正在处理元素：$CURRENT_ELEMENT"
    echo "[步骤 7] 先删除可能残留的 PDOS 中间文件"
    rm -f SELECTED_ATOMS_LIST PDOS_EIG.dat PDOS_EIG_UP.dat PDOS_EIG_DW.dat

    echo "[步骤 7] 调用 VASPKIT 计算元素 $CURRENT_ELEMENT 的 PDOS"
    echo "[步骤 7] 对应临时日志文件（当前目录）：$ELEMENT_PDOS_LOG"
    if ! echo -e "118\n1\n${ENERGY_WINDOW} ${SIGMA}\n1\n${CURRENT_ELEMENT}\n" | "$VASPKIT_BIN" > "$ELEMENT_PDOS_LOG" 2>&1; then
        echo "[ERROR] 元素 $CURRENT_ELEMENT 的 PDOS 计算失败。"
        echo "[ERROR] 请查看日志：$ELEMENT_PDOS_LOG"
        exit 1
    fi

    echo "[步骤 7] 元素 $CURRENT_ELEMENT 的 PDOS 命令执行完成，开始整理输出文件。"

    if [ -f "SELECTED_ATOMS_LIST" ]; then
        mv SELECTED_ATOMS_LIST "${DATA_ELEMENT_DIR}/SELECTED_ATOMS_LIST_${CURRENT_ELEMENT}"
        echo "[步骤 7] 已保存 -> ${DATA_ELEMENT_DIR}/SELECTED_ATOMS_LIST_${CURRENT_ELEMENT}"
    else
        echo "[WARN] 没有找到 SELECTED_ATOMS_LIST，跳过保存。"
    fi

    if [ -f "PDOS_EIG.dat" ]; then
        mv PDOS_EIG.dat "${DATA_ELEMENT_DIR}/PDOS_EIG_${CURRENT_ELEMENT}.dat"
        echo "[步骤 7] 已保存非自旋 PDOS -> ${DATA_ELEMENT_DIR}/PDOS_EIG_${CURRENT_ELEMENT}.dat"

        # 后处理：对 Data-element/PDOS_EIG_元素.dat 生成 .bak 副本，
        # 并在副本末尾追加 sum_s / sum_p / sum_d 三列（与 Data-atom 共用函数）。
        append_sum_columns "${DATA_ELEMENT_DIR}/PDOS_EIG_${CURRENT_ELEMENT}.dat"
    else
        if [ -f "PDOS_EIG_UP.dat" ] && [ -f "PDOS_EIG_DW.dat" ]; then
            mv PDOS_EIG_UP.dat "${DATA_ELEMENT_DIR}/PDOS_EIG_UP_${CURRENT_ELEMENT}.dat"
            mv PDOS_EIG_DW.dat "${DATA_ELEMENT_DIR}/PDOS_EIG_DW_${CURRENT_ELEMENT}.dat"
            echo "[步骤 7] 已保存自旋 PDOS -> ${DATA_ELEMENT_DIR}/PDOS_EIG_UP_${CURRENT_ELEMENT}.dat"
            echo "[步骤 7] 已保存自旋 PDOS -> ${DATA_ELEMENT_DIR}/PDOS_EIG_DW_${CURRENT_ELEMENT}.dat"

            append_sum_columns "${DATA_ELEMENT_DIR}/PDOS_EIG_UP_${CURRENT_ELEMENT}.dat"
            append_sum_columns "${DATA_ELEMENT_DIR}/PDOS_EIG_DW_${CURRENT_ELEMENT}.dat"
        else
            echo "[ERROR] 没有找到元素 $CURRENT_ELEMENT 的 PDOS 输出文件。"
            echo "[ERROR] 请查看日志：$ELEMENT_PDOS_LOG"
            exit 1
        fi
    fi

    # 处理成功：删除临时日志（按需求不在输出目录保留 cmd_*.log）
    rm -f "$ELEMENT_PDOS_LOG"
    echo "[步骤 7] 临时日志已删除：$ELEMENT_PDOS_LOG"

    echo "[步骤 7] 元素 $CURRENT_ELEMENT 的 PDOS 已完成。"
    INDEX=$((INDEX + 1))
done

echo
echo "[步骤 7] 所有元素的 PDOS 都已完成。"
echo

# ------------------------------------------------------------
# 步骤 8：计算每个元素的 d 带中心（VASPKIT 1.5.1 专用格式）
#
# 这一段特意拆成 5 个阶段：
# 1. 构造输入文本
# 2. 打印输入文本
# 3. 执行命令
# 4. 检查输出
# 5. 移动文件
#
# 对于元素 S，这里真正执行的命令应等价于：
# echo -e "503\n1\nN\n1\nS\n" | /opt/vaspkit.1.5.1/bin/vaspkit >> cmd.log 2>&1
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "步骤 8：开始计算每个元素的 d 带中心（适配 VASPKIT 1.5.1）"
echo "------------------------------------------------------------"

INDEX=0
while [ "$INDEX" -lt "$ELEMENT_TYPE_COUNT" ]; do
    CURRENT_ELEMENT="${ELEMENTS[$INDEX]}"

    echo
    echo "============================================================"
    echo "[步骤 8] 当前正在处理元素：$CURRENT_ELEMENT"
    echo "============================================================"

    echo "[步骤 8] 先删除可能残留的 d 带中心中间文件"
    rm -f SELECTED_ATOMS_LIST BAND_CENTER OCCUPATION_NUMBER FERMI_ENERGY cmd.log

    # --------------------------------------------------------
    # 阶段 1/5：构造输入文本
    # --------------------------------------------------------
    echo "[步骤 8 - 阶段 1/5] 开始构造元素 $CURRENT_ELEMENT 的 d 带中心输入文本"
    ELEMENT_DBAND_LINE_1="$DBAND_TASK_CODE"
    ELEMENT_DBAND_LINE_2="$DBAND_FIRST_INPUT"
    ELEMENT_DBAND_LINE_3="$DBAND_CHANGE_WINDOW"
    ELEMENT_DBAND_LINE_4="$DBAND_SELECT_MODE"
    ELEMENT_DBAND_LINE_5="$CURRENT_ELEMENT"
    echo "[步骤 8 - 阶段 1/5] 输入文本构造完成。"

    # --------------------------------------------------------
    # 阶段 2/5：打印输入文本
    # --------------------------------------------------------
    echo "[步骤 8 - 阶段 2/5] 打印即将送入 VASPKIT 的输入文本"
    echo "[步骤 8 - 阶段 2/5] 这一段会严格按下面 5 行输入送给 VASPKIT 1.5.1："
    echo "[步骤 8 - 阶段 2/5]   第 1 行：$ELEMENT_DBAND_LINE_1"
    echo "[步骤 8 - 阶段 2/5]   第 2 行：$ELEMENT_DBAND_LINE_2"
    echo "[步骤 8 - 阶段 2/5]   第 3 行：$ELEMENT_DBAND_LINE_3"
    echo "[步骤 8 - 阶段 2/5]   第 4 行：$ELEMENT_DBAND_LINE_4"
    echo "[步骤 8 - 阶段 2/5]   第 5 行：$ELEMENT_DBAND_LINE_5"
    echo "[步骤 8 - 阶段 2/5] 对应等价命令为："
    echo "[步骤 8 - 阶段 2/5] echo -e \"${ELEMENT_DBAND_LINE_1}\\n${ELEMENT_DBAND_LINE_2}\\n${ELEMENT_DBAND_LINE_3}\\n${ELEMENT_DBAND_LINE_4}\\n${ELEMENT_DBAND_LINE_5}\\n\" | $VASPKIT_BIN >> cmd.log"

    # --------------------------------------------------------
    # 阶段 3/5：执行命令
    # --------------------------------------------------------
    echo "[步骤 8 - 阶段 3/5] 开始执行元素 $CURRENT_ELEMENT 的 d 带中心命令"
    # cmd.log 只是当前目录里的临时日志，成功完成后会删除。
    if ! echo -e "${ELEMENT_DBAND_LINE_1}\n${ELEMENT_DBAND_LINE_2}\n${ELEMENT_DBAND_LINE_3}\n${ELEMENT_DBAND_LINE_4}\n${ELEMENT_DBAND_LINE_5}\n" | "$VASPKIT_BIN" >> cmd.log 2>&1; then
        echo "[ERROR] 元素 $CURRENT_ELEMENT 的 d 带中心计算失败。"
        echo "[ERROR] 请查看当前目录中的日志文件：cmd.log"
        exit 1
    fi
    echo "[步骤 8 - 阶段 3/5] 命令执行完成。"

    # --------------------------------------------------------
    # 阶段 4/5：检查输出
    # --------------------------------------------------------
    echo "[步骤 8 - 阶段 4/5] 开始检查元素 $CURRENT_ELEMENT 的输出文件"

    if [ ! -f "BAND_CENTER" ]; then
        echo "[ERROR] 元素 $CURRENT_ELEMENT 的 d 带中心没有生成 BAND_CENTER 文件。"
        echo "[ERROR] 请查看当前目录中的日志文件：cmd.log"
        exit 1
    fi
    echo "[步骤 8 - 阶段 4/5] 已检测到 BAND_CENTER。"

    if [ -f "SELECTED_ATOMS_LIST" ]; then
        echo "[步骤 8 - 阶段 4/5] 已检测到 SELECTED_ATOMS_LIST。"
    else
        echo "[WARN] [步骤 8 - 阶段 4/5] 没有检测到 SELECTED_ATOMS_LIST。"
    fi

    if [ -f "OCCUPATION_NUMBER" ]; then
        echo "[步骤 8 - 阶段 4/5] 已检测到 OCCUPATION_NUMBER。"
    else
        echo "[WARN] [步骤 8 - 阶段 4/5] 没有检测到 OCCUPATION_NUMBER。"
    fi

    if [ -f "FERMI_ENERGY" ]; then
        echo "[步骤 8 - 阶段 4/5] 已检测到 FERMI_ENERGY。"
    else
        echo "[WARN] [步骤 8 - 阶段 4/5] 没有检测到 FERMI_ENERGY。"
    fi

    echo "[步骤 8 - 阶段 4/5] 将 BAND_CENTER 内容追加到 cmd.log 末尾，便于出错时查看。"
    cat BAND_CENTER >> cmd.log
    echo "[步骤 8 - 阶段 4/5] 输出检查完成。"

    # --------------------------------------------------------
    # 阶段 5/5：移动文件
    # --------------------------------------------------------
    echo "[步骤 8 - 阶段 5/5] 开始移动元素 $CURRENT_ELEMENT 的 d 带中心结果文件"

    if [ -f "SELECTED_ATOMS_LIST" ]; then
        mv SELECTED_ATOMS_LIST "${DATA_ELEMENT_DIR}/SELECTED_ATOMS_LIST_p-dbandcenter_${CURRENT_ELEMENT}"
        echo "[步骤 8 - 阶段 5/5] 已保存 -> ${DATA_ELEMENT_DIR}/SELECTED_ATOMS_LIST_p-dbandcenter_${CURRENT_ELEMENT}"
    else
        echo "[WARN] [步骤 8 - 阶段 5/5] 没有找到 SELECTED_ATOMS_LIST，跳过移动。"
    fi

    mv BAND_CENTER "${DATA_ELEMENT_DIR}/BAND_CENTER_p-dbandcenter_${CURRENT_ELEMENT}"
    echo "[步骤 8 - 阶段 5/5] 已保存 -> ${DATA_ELEMENT_DIR}/BAND_CENTER_p-dbandcenter_${CURRENT_ELEMENT}"

    if [ -f "OCCUPATION_NUMBER" ]; then
        mv OCCUPATION_NUMBER "${DATA_ELEMENT_DIR}/OCCUPATION_NUMBER_p-dbandcenter_${CURRENT_ELEMENT}"
        echo "[步骤 8 - 阶段 5/5] 已保存 -> ${DATA_ELEMENT_DIR}/OCCUPATION_NUMBER_p-dbandcenter_${CURRENT_ELEMENT}"
    else
        echo "[WARN] [步骤 8 - 阶段 5/5] 没有找到 OCCUPATION_NUMBER，跳过移动。"
    fi

    if [ -f "FERMI_ENERGY" ]; then
        mv FERMI_ENERGY "${DATA_ELEMENT_DIR}/FERMI_ENERGY_p-dbandcenter_${CURRENT_ELEMENT}"
        echo "[步骤 8 - 阶段 5/5] 已保存 -> ${DATA_ELEMENT_DIR}/FERMI_ENERGY_p-dbandcenter_${CURRENT_ELEMENT}"
    else
        echo "[WARN] [步骤 8 - 阶段 5/5] 没有找到 FERMI_ENERGY，跳过移动。"
    fi

    # 按需求：cmd.log 不在输出目录保留，成功完成后直接删除。
    if [ -f "cmd.log" ]; then
        rm -f cmd.log
        echo "[步骤 8 - 阶段 5/5] 临时日志 cmd.log 已删除（不保留在输出目录）。"
    fi

    echo "[步骤 8] 元素 $CURRENT_ELEMENT 的 d 带中心已完成。"
    INDEX=$((INDEX + 1))
done

echo
echo "[步骤 8] 所有元素的 d 带中心都已完成。"
echo

# ------------------------------------------------------------
# 步骤 9：根据 RUN_ATOM_PDOS 的值，决定是否计算每个原子的 PDOS
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "步骤 9：决定是否计算每个原子的 PDOS"
echo "------------------------------------------------------------"

if [ "$RUN_ATOM_PDOS" = "yes" ]; then
    echo "[步骤 9] 你选择了 yes，因此现在开始计算每个原子的 PDOS。"
    echo

    ATOM_INDEX=1
    while [ "$ATOM_INDEX" -le "$TOTAL_NAT" ]; do
        # 注意：日志写到【当前目录】的临时文件，成功完成后删除，
        # 因此 Data-atom 目录中不会保留 cmd_*.log。
        ATOM_PDOS_LOG="cmd_PDOS_${ATOM_INDEX}.log"

        echo "[步骤 9] 当前正在处理原子编号：$ATOM_INDEX"
        echo "[步骤 9] 先删除可能残留的 PDOS 中间文件"
        rm -f SELECTED_ATOMS_LIST PDOS_EIG.dat PDOS_EIG_UP.dat PDOS_EIG_DW.dat

        echo "[步骤 9] 调用 VASPKIT 计算原子 $ATOM_INDEX 的 PDOS"
        echo "[步骤 9] 对应临时日志文件（当前目录）：$ATOM_PDOS_LOG"
        if ! echo -e "118\n1\n${ENERGY_WINDOW} ${SIGMA}\n1\n${ATOM_INDEX}\n" | "$VASPKIT_BIN" > "$ATOM_PDOS_LOG" 2>&1; then
            echo "[ERROR] 原子 $ATOM_INDEX 的 PDOS 计算失败。"
            echo "[ERROR] 请查看日志：$ATOM_PDOS_LOG"
            exit 1
        fi

        if [ -f "SELECTED_ATOMS_LIST" ]; then
            mv SELECTED_ATOMS_LIST "${DATA_ATOM_DIR}/SELECTED_ATOMS_LIST_${ATOM_INDEX}"
            echo "[步骤 9] 已保存 -> ${DATA_ATOM_DIR}/SELECTED_ATOMS_LIST_${ATOM_INDEX}"
        else
            echo "[WARN] 没有找到 SELECTED_ATOMS_LIST，跳过保存。"
        fi

        if [ -f "PDOS_EIG.dat" ]; then
            mv PDOS_EIG.dat "${DATA_ATOM_DIR}/PDOS_EIG_${ATOM_INDEX}.dat"
            echo "[步骤 9] 已保存非自旋 PDOS -> ${DATA_ATOM_DIR}/PDOS_EIG_${ATOM_INDEX}.dat"

            # 后处理：与 Data-element 一样，生成 .bak 副本，
            # 并在副本末尾追加 sum_s / sum_p / sum_d 三列。
            append_sum_columns "${DATA_ATOM_DIR}/PDOS_EIG_${ATOM_INDEX}.dat"
        else
            if [ -f "PDOS_EIG_UP.dat" ] && [ -f "PDOS_EIG_DW.dat" ]; then
                mv PDOS_EIG_UP.dat "${DATA_ATOM_DIR}/PDOS_EIG_UP_${ATOM_INDEX}.dat"
                mv PDOS_EIG_DW.dat "${DATA_ATOM_DIR}/PDOS_EIG_DW_${ATOM_INDEX}.dat"
                echo "[步骤 9] 已保存自旋 PDOS -> ${DATA_ATOM_DIR}/PDOS_EIG_UP_${ATOM_INDEX}.dat"
                echo "[步骤 9] 已保存自旋 PDOS -> ${DATA_ATOM_DIR}/PDOS_EIG_DW_${ATOM_INDEX}.dat"

                append_sum_columns "${DATA_ATOM_DIR}/PDOS_EIG_UP_${ATOM_INDEX}.dat"
                append_sum_columns "${DATA_ATOM_DIR}/PDOS_EIG_DW_${ATOM_INDEX}.dat"
            else
                echo "[ERROR] 没有找到原子 $ATOM_INDEX 的 PDOS 输出文件。"
                echo "[ERROR] 请查看日志：$ATOM_PDOS_LOG"
                exit 1
            fi
        fi

        # 处理成功：删除临时日志（按需求不在输出目录保留 cmd_*.log）
        rm -f "$ATOM_PDOS_LOG"
        echo "[步骤 9] 临时日志已删除：$ATOM_PDOS_LOG"

        echo "[步骤 9] 原子 $ATOM_INDEX 的 PDOS 已完成。"
        echo

        ATOM_INDEX=$((ATOM_INDEX + 1))
    done

    echo "[步骤 9] 所有原子的 PDOS 都已完成。"
else
    echo "[步骤 9] 你选择了 no，因此跳过每个原子的 PDOS 计算。"
fi

echo

# ------------------------------------------------------------
# 步骤 10：计算总 TDOS
# 输出放到 Data-total
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "步骤 10：开始计算总 TDOS"
echo "------------------------------------------------------------"

echo "[步骤 10] 先删除可能残留的 TDOS 中间文件"
rm -f TDOS.dat TDOS_UP.dat TDOS_DW.dat FERMI_ENERGY

# 注意：日志写到【当前目录】的临时文件，成功完成后删除，
# 因此 Data-total 目录中不会保留 cmd_*.log。
TDOS_LOG="cmd_TDOS.log"

echo "[步骤 10] 调用 VASPKIT 计算总 TDOS"
echo "[步骤 10] 对应临时日志文件（当前目录）：$TDOS_LOG"
if ! echo -e "117\n1\n${ENERGY_WINDOW} ${SIGMA}\n" | "$VASPKIT_BIN" > "$TDOS_LOG" 2>&1; then
    echo "[ERROR] TDOS 计算失败。"
    echo "[ERROR] 请查看日志：$TDOS_LOG"
    exit 1
fi

TDOS_FILE_COUNT=0

if [ -f "TDOS.dat" ]; then
    mv TDOS.dat "${DATA_TOTAL_DIR}/TDOS.dat"
    echo "[步骤 10] 已保存 -> ${DATA_TOTAL_DIR}/TDOS.dat"
    TDOS_FILE_COUNT=$((TDOS_FILE_COUNT + 1))
fi

if [ -f "TDOS_UP.dat" ]; then
    mv TDOS_UP.dat "${DATA_TOTAL_DIR}/TDOS_UP.dat"
    echo "[步骤 10] 已保存 -> ${DATA_TOTAL_DIR}/TDOS_UP.dat"
    TDOS_FILE_COUNT=$((TDOS_FILE_COUNT + 1))
fi

if [ -f "TDOS_DW.dat" ]; then
    mv TDOS_DW.dat "${DATA_TOTAL_DIR}/TDOS_DW.dat"
    echo "[步骤 10] 已保存 -> ${DATA_TOTAL_DIR}/TDOS_DW.dat"
    TDOS_FILE_COUNT=$((TDOS_FILE_COUNT + 1))
fi

if [ "$TDOS_FILE_COUNT" -eq 0 ]; then
    echo "[ERROR] 没有找到任何 TDOS 输出文件。"
    echo "[ERROR] 请查看日志：$TDOS_LOG"
    exit 1
fi

if [ -f "FERMI_ENERGY" ]; then
    mv FERMI_ENERGY "${DATA_TOTAL_DIR}/FERMI_ENERGY"
    echo "[步骤 10] 已保存 -> ${DATA_TOTAL_DIR}/FERMI_ENERGY"
else
    echo "[WARN] 没有找到 FERMI_ENERGY，跳过保存。"
fi

# 处理成功：删除临时日志（按需求不在输出目录保留 cmd_*.log）
rm -f "$TDOS_LOG"
echo "[步骤 10] 临时日志已删除：$TDOS_LOG"

echo "[步骤 10] 总 TDOS 已完成。"
echo

# ------------------------------------------------------------
# 步骤 11：计算总 d 带中心（VASPKIT 1.5.1 专用格式）
#
# 这里也拆成 5 个阶段：
# 1. 构造输入文本
# 2. 打印输入文本
# 3. 执行命令
# 4. 检查输出
# 5. 移动文件
#
# 这里真正执行的命令应等价于：
# echo -e "503\n1\nN\n1\nall\n" | /opt/vaspkit.1.5.1/bin/vaspkit >> cmd.log 2>&1
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "步骤 11：开始计算总 d 带中心（适配 VASPKIT 1.5.1）"
echo "------------------------------------------------------------"

echo "[步骤 11] 先删除可能残留的 d 带中心中间文件"
rm -f SELECTED_ATOMS_LIST BAND_CENTER OCCUPATION_NUMBER FERMI_ENERGY cmd.log

# --------------------------------------------------------
# 阶段 1/5：构造输入文本
# --------------------------------------------------------
echo "[步骤 11 - 阶段 1/5] 开始构造总 d 带中心输入文本"
TOTAL_DBAND_LINE_1="$DBAND_TASK_CODE"
TOTAL_DBAND_LINE_2="$DBAND_FIRST_INPUT"
TOTAL_DBAND_LINE_3="$DBAND_CHANGE_WINDOW"
TOTAL_DBAND_LINE_4="$DBAND_SELECT_MODE"
TOTAL_DBAND_LINE_5="all"
echo "[步骤 11 - 阶段 1/5] 输入文本构造完成。"

# --------------------------------------------------------
# 阶段 2/5：打印输入文本
# --------------------------------------------------------
echo "[步骤 11 - 阶段 2/5] 打印即将送入 VASPKIT 的输入文本"
echo "[步骤 11 - 阶段 2/5] 这一段会严格按下面 5 行输入送给 VASPKIT 1.5.1："
echo "[步骤 11 - 阶段 2/5]   第 1 行：$TOTAL_DBAND_LINE_1"
echo "[步骤 11 - 阶段 2/5]   第 2 行：$TOTAL_DBAND_LINE_2"
echo "[步骤 11 - 阶段 2/5]   第 3 行：$TOTAL_DBAND_LINE_3"
echo "[步骤 11 - 阶段 2/5]   第 4 行：$TOTAL_DBAND_LINE_4"
echo "[步骤 11 - 阶段 2/5]   第 5 行：$TOTAL_DBAND_LINE_5"
echo "[步骤 11 - 阶段 2/5] 对应等价命令为："
echo "[步骤 11 - 阶段 2/5] echo -e \"${TOTAL_DBAND_LINE_1}\\n${TOTAL_DBAND_LINE_2}\\n${TOTAL_DBAND_LINE_3}\\n${TOTAL_DBAND_LINE_4}\\n${TOTAL_DBAND_LINE_5}\\n\" | $VASPKIT_BIN >> cmd.log"

# --------------------------------------------------------
# 阶段 3/5：执行命令
# --------------------------------------------------------
echo "[步骤 11 - 阶段 3/5] 开始执行总 d 带中心命令"
# cmd.log 只是当前目录里的临时日志，成功完成后会删除。
if ! echo -e "${TOTAL_DBAND_LINE_1}\n${TOTAL_DBAND_LINE_2}\n${TOTAL_DBAND_LINE_3}\n${TOTAL_DBAND_LINE_4}\n${TOTAL_DBAND_LINE_5}\n" | "$VASPKIT_BIN" >> cmd.log 2>&1; then
    echo "[ERROR] 总 d 带中心计算失败。"
    echo "[ERROR] 请查看当前目录中的日志文件：cmd.log"
    exit 1
fi
echo "[步骤 11 - 阶段 3/5] 命令执行完成。"

# --------------------------------------------------------
# 阶段 4/5：检查输出
# --------------------------------------------------------
echo "[步骤 11 - 阶段 4/5] 开始检查总 d 带中心输出文件"
if [ ! -f "BAND_CENTER" ]; then
    echo "[ERROR] 总 d 带中心没有生成 BAND_CENTER 文件。"
    echo "[ERROR] 请查看当前目录中的日志文件：cmd.log"
    exit 1
fi
echo "[步骤 11 - 阶段 4/5] 已检测到 BAND_CENTER。"

echo "[步骤 11 - 阶段 4/5] 将 BAND_CENTER 内容追加到 cmd.log 末尾，便于出错时查看。"
cat BAND_CENTER >> cmd.log
echo "[步骤 11 - 阶段 4/5] 输出检查完成。"

# --------------------------------------------------------
# 阶段 5/5：移动文件
# --------------------------------------------------------
echo "[步骤 11 - 阶段 5/5] 开始移动总 d 带中心结果文件"
mv BAND_CENTER "${DATA_TOTAL_DIR}/BAND_CENTER"
echo "[步骤 11 - 阶段 5/5] 已保存 -> ${DATA_TOTAL_DIR}/BAND_CENTER"

# 按需求：cmd.log 不在输出目录保留，成功完成后直接删除。
if [ -f "cmd.log" ]; then
    rm -f cmd.log
    echo "[步骤 11 - 阶段 5/5] 临时日志 cmd.log 已删除（不保留在输出目录）。"
fi

echo "[步骤 11] 总 d 带中心已完成。"
echo

# ------------------------------------------------------------
# 步骤 12：兜底清理输出目录中的 cmd_*.log，并打印最终总结
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "步骤 12：兜底清理输出目录中的 cmd_*.log"
echo "------------------------------------------------------------"

CLEANED_COUNT=0
for CLEAN_DIR in "$DATA_ELEMENT_DIR" "$DATA_ATOM_DIR" "$DATA_TOTAL_DIR"; do
    if [ ! -d "$CLEAN_DIR" ]; then
        continue
    fi
    for CLEAN_LOG in "$CLEAN_DIR"/cmd_*.log; do
        if [ -f "$CLEAN_LOG" ]; then
            echo "[步骤 12] 删除残留日志：$CLEAN_LOG"
            rm -f "$CLEAN_LOG"
            CLEANED_COUNT=$((CLEANED_COUNT + 1))
        fi
    done
done

if [ "$CLEANED_COUNT" -eq 0 ]; then
    echo "[步骤 12] 三个输出目录中没有发现 cmd_*.log，无需清理。"
else
    echo "[步骤 12] 共删除 $CLEANED_COUNT 个残留日志。"
fi

echo
echo "============================================================"
echo "脚本执行完成：cdos_1.5.1.sh"
echo "============================================================"
echo "[总结] 使用的 VASPKIT = $VASPKIT_BIN"
echo "[总结] 分元素 PDOS 目录 = $DATA_ELEMENT_DIR"
if [ "$RUN_ATOM_PDOS" = "yes" ]; then
    echo "[总结] 每原子 PDOS 目录 = $DATA_ATOM_DIR"
else
    echo "[总结] 每原子 PDOS = 本次已跳过"
fi
echo "[总结] 总 TDOS / 总 d 带中心目录 = $DATA_TOTAL_DIR"
echo "[总结] 已完成：分元素 PDOS、分元素 d 带中心、总 TDOS、总 d 带中心"
echo "[总结] Data-element 和 Data-atom 的每个 PDOS 都已生成 .bak 副本，并追加 sum_s / sum_p / sum_d 三列"
echo "[总结] 三个输出目录中已清理全部 cmd_*.log 日志"
echo "[总结] 本脚本的 d 带中心输入格式已固定为 VASPKIT 1.5.1 版本格式：503 -> 1 -> N -> 1 -> 目标"
echo "============================================================"
