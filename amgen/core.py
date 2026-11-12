"""amgen.core — Wyckoff 轨道分析 + 补偿磁染色枚举（群置换 + 全局翻号去重）"""
import itertools
import numpy as np
import spglib
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def find_orbits(struct, symprec=1e-3):
    """返回 {wyckoff_letter: [site_indices]}（对称化结构口径）。

    这是按字母合并的旧接口；同一字母被多个独立轨道复用时会被合并。
    需要轨道实例级枚举时用 find_orbit_instances。
    """
    ss = SpacegroupAnalyzer(struct, symprec=symprec).get_symmetrized_structure()
    groups = {}
    for i, wk in enumerate(ss.wyckoff_letters):
        groups.setdefault(wk, []).append(i)
    return groups


def find_orbit_instances(struct, symprec=1e-3):
    """返回独立 Wyckoff 轨道实例，而不是按字母合并的轨道。

    同一 Wyckoff 字母可以被不同元素/不同自由参数复用；此时它们是两个
    独立实例。pymatgen 的 ``equivalent_indices`` 正是按实际等价关系分组，
    因此这里以它为口径。返回顺序按最小位点索引稳定排序，便于复现。
    """
    ss = SpacegroupAnalyzer(struct, symprec=symprec).get_symmetrized_structure()
    letters = ss.wyckoff_letters
    instances = []
    for raw in ss.equivalent_indices:
        idxs = tuple(sorted(int(i) for i in raw))
        instances.append({
            "letter": str(letters[idxs[0]]),
            "indices": idxs,
            "element": str(struct[idxs[0]].specie.symbol),
            "multiplicity": len(idxs),
        })
    instances.sort(key=lambda x: (min(x["indices"]), x["letter"]))
    return instances


def site_perms(coords, lattice, symprec=1e-4, species_numbers=None):
    """空间群操作在给定位点集合上的置换表示（用于染色去重）。

    species_numbers 为 None 时按同种原子处理，得到抽象几何框架的对称操作；
    传入原子序数列表时，只保留同时保持化学种类的操作。角色枚举需要后者，
    否则不同元素的轨道会被错误地视为可互换。
    """
    n = len(coords)
    if species_numbers is None:
        species_numbers = [25] * n
    cell = (lattice.matrix, [list(c) for c in coords], list(species_numbers))
    ds = spglib.get_symmetry_dataset(cell, symprec=symprec)
    if ds is None:
        return [np.arange(n)]
    perms = set()
    for rot, tra in zip(ds.rotations, ds.translations):
        p, ok = [], True
        for i in range(n):
            x = rot @ coords[i] + tra
            hit = None
            for k in range(n):
                d = x - coords[k]
                d -= np.round(d)
                if np.linalg.norm(d @ lattice.matrix) < 0.05:
                    hit = k
                    break
            if hit is None:
                ok = False
                break
            p.append(hit)
        if ok:
            perms.add(tuple(p))
    return [np.array(p) for p in perms]


def canonical(signs, perms):
    """染色的 canonical 代表：群置换 × 全局翻号下取最小元"""
    s = np.array(signs)
    best = None
    for p in perms:
        for flip in (1, -1):
            t = tuple(flip * s[p])
            if best is None or t < best:
                best = t
    return best


def orbit_colorings(coords, lattice):
    """一个轨道的全部不等价补偿染色（+/- 各半），返回 pattern 字符串列表"""
    n = len(coords)
    assert n % 2 == 0, "奇数重轨道无法单独补偿"
    perms = site_perms(coords, lattice)
    seen, pats = set(), []
    for idxs in itertools.combinations(range(n), n // 2):
        signs = tuple(1 if i in idxs else -1 for i in range(n))
        c = canonical(signs, perms)
        if c in seen:
            continue
        seen.add(c)
        pats.append("".join("+" if x > 0 else "-" for x in signs))
    return pats


def joint_colorings(struct, mag_elements, max_mult=12, max_combos=300):
    """多轨道联合染色：磁性元素占据的偶数重轨道各枚举染色，跨轨道笛卡尔积。

    按独立轨道实例枚举，同一 Wyckoff 字母的复用实例分开处理。
    产出 (desc, odd_note, mag_idx, signs) 迭代器；组合爆炸或无可染色轨道时
    产出 ("SKIP", reason, None, None)。
    """
    orbits = find_orbit_instances(struct)
    even, odd = [], []
    for o in orbits:
        wk, idxs, el = o["letter"], o["indices"], o["element"]
        if el not in mag_elements:
            continue
        if len(idxs) % 2 == 0 and 2 <= len(idxs) <= max_mult:
            even.append((wk, el, idxs))
        else:
            odd.append(f"{el}({wk})^{len(idxs)}")
    odd_note = ";".join(odd)
    if not even:
        yield ("SKIP", "无可染色偶数重轨道|" + odd_note, None, None)
        return
    pat_lists = []
    for wk, el, idxs in even:
        coords = [struct[i].frac_coords for i in idxs]
        pat_lists.append((wk, el, idxs, orbit_colorings(coords, struct.lattice)))
    n = int(np.prod([len(p[3]) for p in pat_lists]))
    if n > max_combos:
        yield ("SKIP", f"组合爆炸 {n}|" + odd_note, None, None)
        return
    for combo in itertools.product(*[p[3] for p in pat_lists]):
        midx, signs, desc = [], [], []
        for (wk, el, idxs, _), pat in zip(pat_lists, combo):
            midx += idxs
            signs += [1 if c == "+" else -1 for c in pat]
            desc.append(f"{el}({wk})^{len(idxs)}:{pat}")
        yield ("+".join(desc), odd_note, midx, signs)
