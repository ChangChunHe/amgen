"""amgen.core — Wyckoff 轨道分析 + 补偿磁染色枚举（群置换 + 全局翻号去重）"""
import itertools
import numpy as np
import spglib
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def find_orbits(struct, symprec=1e-3):
    """返回 {wyckoff_letter: [site_indices]}（对称化结构口径）"""
    ss = SpacegroupAnalyzer(struct, symprec=symprec).get_symmetrized_structure()
    groups = {}
    for i, wk in enumerate(ss.wyckoff_letters):
        groups.setdefault(wk, []).append(i)
    return groups


def site_perms(coords, lattice, symprec=1e-4):
    """空间群操作在给定位点集合上的置换表示（用于染色去重）"""
    n = len(coords)
    cell = (lattice.matrix, [list(c) for c in coords], [25] * n)
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

    产出 (desc, mag_idx, signs, odd_note) 迭代器；组合爆炸或无可染色轨道时
    产出 ("SKIP", reason, None, None)。
    """
    orbits = find_orbits(struct)
    even, odd = [], []
    for wk, idxs in sorted(orbits.items()):
        el = str(struct[idxs[0]].specie.symbol)
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
