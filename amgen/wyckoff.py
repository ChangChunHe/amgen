"""amgen.wyckoff — 群分解: SG -> Wyckoff 表查询 + 按配方生成结构（pyxtal 后端）

配方格式: {元素: [字母, ...]}, 同一字母可被不同元素重复使用
（如 Cs3M2X9 中 Cs@f 与 M@f 是两个独立轨道实例）。
"""
import itertools
import numpy as np
from pyxtal import pyxtal
from pyxtal.symmetry import Group


def wyckoff_table(sg):
    """SG 号 -> [(letter, multiplicity, n_free_params)]，按重数升序"""
    g = Group(sg)
    table = [(w.letter, w.multiplicity, w.get_dof()) for w in g.Wyckoff_positions]
    return sorted(table, key=lambda t: t[1])


def _mult_map(sg):
    return {L: m for L, m, _ in wyckoff_table(sg)}


def generate(sg, recipe, seed=0):
    """按配方生成结构。

    recipe: {'Mn': ['f'], 'Ni': ['a', 'c'], ...} —— 字母为 Wyckoff 字母。
    返回 pymatgen.Structure（自由参数为 pyxtal 随机值，晶格为随机六角/三方格子）。
    """
    mult = _mult_map(sg)
    species, sites, numions = [], [], []
    for el, letters in recipe.items():
        species.append(el)
        sites.append([f"{mult[L]}{L}" for L in letters])
        numions.append(sum(mult[L] for L in letters))
    x = pyxtal()
    x.from_random(3, sg, species, numions, sites=sites, random_state=seed)
    return x.to_pymatgen()


def min_distance(struct):
    """全结构最小原子间距（Å）"""
    dm = struct.distance_matrix
    dm[dm == 0] = np.inf
    return float(dm.min())


def prescale_volume(struct, packing=0.65):
    """按共价半径估算体积并等比缩放晶格（CHGNet 前的必要步骤:
    pyxtal 随机晶格的体积可能与合理值差 50%, 粗弛豫救不回来）"""
    from .decorate import COV
    v = sum(4.0 / 3.0 * np.pi * COV.get(str(x.specie.symbol), 1.4) ** 3
            for x in struct) / packing
    struct.scale_lattice(v)
    return struct


def generate_ok(sg, recipe, seeds=range(20), min_d=1.6):
    """带距离检查的种子重试 + 体积预标定: 返回 (结构 or None)"""
    last = None
    for seed in seeds:
        try:
            s = generate(sg, recipe, seed=seed)
            prescale_volume(s)
            last = s
            if min_distance(s) >= min_d:
                return s
        except Exception:
            continue
    return last


def orbit_combos(sg, max_atoms=30, need_even_orbit=True, max_orbits=6):
    """枚举 SG 下的 Wyckoff 轨道组合（字母可重复，总重数 <= max_atoms）。

    返回 [((letter, mult), ...), ...]，要求至少含一个偶数重轨道（磁染色用）。
    组合按总原子数升序产出，供装饰/染色漏斗使用。
    """
    table = [(L, m) for L, m, _ in wyckoff_table(sg)]
    combos = []
    for k in range(1, max_orbits + 1):
        for combo in itertools.combinations_with_replacement(table, k):
            n = sum(m for _, m in combo)
            if n > max_atoms:
                continue
            if need_even_orbit and not any(m % 2 == 0 and m >= 4 for _, m in combo):
                continue
            combos.append(combo)
    return sorted(combos, key=lambda c: sum(m for _, m in c))
