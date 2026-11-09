"""amgen.decorate — 化学装饰：轨道级元素替换 + 价态/配位过滤 + 晶格标定"""
from collections import deque
import numpy as np

MAG_SET_3D = ["V", "Cr", "Mn", "Fe", "Co", "Ni"]
MAG_SET_4D = ["Nb", "Mo", "Ru", "Rh", "Pd"]
MAG_SET_5D = ["Ta", "W", "Re", "Os", "Ir"]
MAG_SET = MAG_SET_3D + MAG_SET_4D + MAG_SET_5D

GRP_NUM = {"V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10,
           "Nb": 5, "Mo": 6, "Ru": 8, "Rh": 9, "Pd": 10,
           "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9}
OS_MAG = {"V": [2, 3, 4, 5], "Cr": [2, 3, 4, 6], "Mn": [2, 3, 4],
          "Fe": [2, 3], "Co": [2, 3], "Ni": [2],
          "Nb": [3, 5], "Mo": [4], "Ru": [3, 4], "Rh": [3], "Pd": [2],
          "Ta": [5], "W": [4], "Re": [4], "Os": [4], "Ir": [3, 4]}
MOMENT = {"V": 2.0, "Cr": 3.0, "Mn": 3.5, "Fe": 3.0, "Co": 2.0, "Ni": 1.5,
          "Nb": 3.0, "Mo": 3.0, "Ru": 2.5, "Rh": 2.0, "Pd": 1.0,
          "Ta": 3.0, "W": 3.0, "Re": 2.5, "Os": 2.0, "Ir": 1.5}
COV = {"H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
       "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 0.99, "Ge": 1.20, "As": 1.19,
       "Se": 1.20, "Br": 1.14, "Te": 1.38, "I": 1.33,
       "Li": 1.28, "Na": 1.66, "K": 2.03, "Rb": 2.16, "Cs": 2.32, "Mg": 1.41, "Ca": 1.76,
       "Al": 1.21, "Ga": 1.22, "Sr": 1.95, "Ba": 1.98,
       "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.50, "Fe": 1.42,
       "Co": 1.38, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
       "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42,
       "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
       "La": 2.07, "Ce": 2.04, "Tb": 1.94, "Hf": 1.87, "Ta": 1.70, "W": 1.62,
       "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Th": 2.06}
OS_COMMON = {
    "Li": 1, "Na": 1, "K": 1, "Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2,
    "B": 3, "Al": 3, "Ga": 3, "C": -4, "Si": -4, "Ge": 2, "Sn": 2,
    "N": -3, "P": -3, "As": -3, "Sb": -3, "Bi": 3,
    "O": -2, "S": -2, "Se": -2, "Te": -2, "F": -1, "Cl": -1, "Br": -1, "I": -1,
    "Cu": 2, "Ag": 1, "Au": 3, "Zn": 2, "Cd": 2,
    "Sc": 3, "Y": 3, "La": 3, "Ce": 3, "Tb": 3, "Zr": 4, "Hf": 4, "Th": 4,
    "Ti": 4, "Nb": 5, "Ta": 5, "Mo": 4, "W": 4,
    "V": 3, "Cr": 3, "Mn": 2, "Fe": 2, "Co": 2, "Ni": 2,
    "Ru": 3, "Rh": 3, "Pd": 2, "Re": 4, "Os": 4, "Ir": 3}
ANIONS = {"O", "S", "Se", "Te", "F", "Cl", "Br", "I", "N", "P", "As", "Sb"}


def substitute_orbits(struct, orbit_map):
    """按 {wyckoff_letter: 新元素} 逐轨道替换（整轨道所有位点同元素）"""
    from .core import find_orbits
    s = struct.copy()
    orbits = find_orbits(struct)
    for wk, el in orbit_map.items():
        for i in orbits[wk]:
            s.replace(i, el)
    return s


def cn_of(struct, idx, exclude_el=None, rcut=3.2):
    """位点 idx 的配位数（可排除同元素近邻）"""
    ns = struct.get_neighbors(struct[idx], rcut)
    if exclude_el:
        ns = [n for n in ns if str(n.specie.symbol) != exclude_el]
    return len(ns)


def valence_charge_need(struct, mag_idx, nM):
    """非着色位点固定价态之和 -> 每个磁性位需要的氧化态"""
    fixed = 0.0
    idxset = set(mag_idx)
    for i in range(len(struct)):
        if i not in idxset:
            fixed += OS_COMMON.get(str(struct[i].specie.symbol), 0)
    return -fixed / nM


def os_ok(M, need):
    """元素 M 是否存在满足电荷中性且 d 电子数 2..8 的氧化态"""
    for o in OS_MAG.get(M, []):
        if abs(o - need) < 0.51 and 2 <= GRP_NUM[M] - o <= 8:
            return o
    return None


def scale_by_radius(struct, old_el, new_el):
    """按共价半径比均匀缩放晶格体积"""
    s = struct.copy()
    f = COV.get(new_el, 1.4) / COV.get(old_el, 1.4)
    s.scale_lattice(s.volume * f ** 3)
    return s, f


def afm_g_pattern(struct, midx, rcut=3.2):
    """最近邻二部图 AFM 染色（G 型参考构型），返回 ±1 数组"""
    pos = struct.cart_coords
    nM = len(midx)
    adj = [[] for _ in range(nM)]
    inv = np.linalg.inv(struct.lattice.matrix)
    for a in range(nM):
        for b in range(a + 1, nM):
            dd = pos[midx[a]] - pos[midx[b]]
            dd -= np.round(dd @ inv) @ struct.lattice.matrix
            if np.linalg.norm(dd) < rcut:
                adj[a].append(b)
                adj[b].append(a)
    color = np.full(nM, 99)
    for r0 in range(nM):
        if color[r0] != 99:
            continue
        color[r0] = 1
        dq = deque([r0])
        while dq:
            u = dq.popleft()
            for v in adj[u]:
                if color[v] == 99:
                    color[v] = -color[u]
                    dq.append(v)
                elif color[v] == color[u]:
                    color[v] = -color[u]
    color[color == 99] = 1
    return color


def compress_magmom(mag):
    """[0,0,5,5,-5,-5,...] -> '6*0 2*5 2*-5 ...' 紧凑写法"""
    parts, i = [], 0
    while i < len(mag):
        j = i
        while j + 1 < len(mag) and abs(mag[j + 1] - mag[i]) < 1e-8:
            j += 1
        n = j - i + 1
        parts.append(f"{n}*{mag[i]:g}" if n > 1 else f"{mag[i]:g}")
        i = j + 1
    return " ".join(parts)


def full_magmom(n_sites, mag_idx, signs, moment=2.0):
    """生成全原子长度的 MAGMOM 列表（非磁位 0），moment 为标量 μB"""
    arr = np.zeros(n_sites)
    for i, sg in zip(mag_idx, signs):
        arr[i] = moment * sg
    return arr
