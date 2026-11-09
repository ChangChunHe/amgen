#!/usr/bin/env python3
"""纯枚举版复现验证：群分解配方 -> 生成 -> (前染色 fsg 门) -> CHGNet 预弛豫
-> 对称化回母群 -> (复判 fsg 门) -> 与已确认 DFT 结构对照

四个家族配方见 amgen/recipes.py（Mg6Mn6Ni6 / BaFe2O4 / Cs3Cr2Br9 / Mg2Cr2B2O7）。

用法: python reproduce_families.py [--write-vasp OUTDIR]
"""
import sys
from collections import Counter, defaultdict
from pathlib import Path
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from amgen import core, classify, decorate, export, wyckoff, relax
from amgen.recipes import RECIPES

SM = StructureMatcher(ltol=0.2, stol=0.5, angle_tol=10, primitive_cell=False)


def framework_match(a, b, symprec=1e-2, lat_tol=0.25):
    """同原型判定：同 SG + 各元素轨道重数谱相同 + 晶格常数接近（原点选择不敏感）"""
    sa = SpacegroupAnalyzer(a, symprec)
    sb = SpacegroupAnalyzer(b, symprec)
    if sa.get_space_group_number() != sb.get_space_group_number():
        return False

    def sig(s):
        ss = SpacegroupAnalyzer(s, symprec).get_symmetrized_structure()
        d = defaultdict(Counter)
        for i, w in enumerate(ss.wyckoff_letters):
            d[str(s[i].specie.symbol)][w] += 1
        return {el: sorted(c.values()) for el, c in d.items()}

    if sig(a) != sig(b):
        return False
    for x, y in ((a.lattice.a, b.lattice.a), (a.lattice.c, b.lattice.c)):
        if abs(x / y - 1) > lat_tol:
            return False
    return True


def wave_gate(struct, mag_el, want_wave, tag):
    """在磁元素轨道上枚举染色并 fsg 判定；返回 (通过的pattern, idxs, signs) 或 None"""
    orbits = core.find_orbits(struct)
    hits = []
    for wk, idxs in sorted(orbits.items()):
        if str(struct[idxs[0]].specie.symbol) != mag_el:
            continue
        if len(idxs) % 2:
            continue
        coords = [struct[i].frac_coords for i in idxs]
        for pat in core.orbit_colorings(coords, struct.lattice):
            signs = [1 if c == "+" else -1 for c in pat]
            r = classify.fsg_classify(struct, idxs, signs)
            if r["texture"] == want_wave:
                hits.append((f"{mag_el}({wk})^{len(idxs)}:{pat}", idxs, signs))
    if hits:
        for desc, _, _ in hits:
            print(f"  [{tag}] {want_wave} 命中: {desc}")
        return hits[0]
    print(f"  [{tag}] 未找到 {want_wave} 染色!")
    return None


def run_family(name, rc, write=None):
    print("=" * 64)
    print(f"{name}: SG{rc['sg']} 配方 {rc['recipe']}  目标 {rc['target_wave']}")

    # 1. 群分解生成
    raw = wyckoff.generate_ok(rc["sg"], rc["recipe"], min_d=1.5)
    if raw is None:
        print("  生成失败 (pyxtal 全部种子碰撞)")
        return False
    print(f"  生成: {raw.composition.reduced_formula}, n={len(raw)}, "
          f"min_d={wyckoff.min_distance(raw):.2f} Å")

    # 2. 前染色门（未弛豫框架上）
    g1 = wave_gate(raw, rc["mag_el"], rc["target_wave"], "前染色")
    if g1 is None:
        return False

    # 3. CHGNet 预弛豫 + 对称化回母群
    s2, sg2 = relax.relax_and_symmetrize(raw, rc["sg"])
    if s2 is None:
        print(f"  弛豫后无法回到 SG{rc['sg']} (终态 SG={sg2})")
        return False
    print(f"  CHGNet 弛豫完成, 对称化回 SG{sg2}, min_d={wyckoff.min_distance(s2):.2f} Å")

    # 4. 复判门（弛豫后结构上重新染色）
    g2 = wave_gate(s2, rc["mag_el"], rc["target_wave"], "复判")
    if g2 is None:
        return False

    # 5. 与已确认结构对照
    ref = Structure.from_file(rc["ref"])
    if SM.fit(s2, ref):
        print("  [PASS] 与参考结构 StructureMatcher 一致")
        ok = True
    elif framework_match(s2, ref):
        print("  [PASS] 与参考结构同原型框架")
        ok = True
    else:
        print("  [FAIL] 与参考结构不匹配")
        ok = False

    # 6. 可选: 落盘 DFT 输入（用复判命中的染色）
    if write and g2 is not None:
        desc, idxs, signs = g2
        mag = decorate.full_magmom(len(s2), idxs, signs,
                                   decorate.MOMENT.get(rc["mag_el"], 3.0))
        export.write_vasp(Path(write) / name / rc["target_wave"].split("-")[0],
                          s2, mag)
        print(f"  DFT 输入已写入 {write}/{name}/ ({desc})")
    return ok


if __name__ == "__main__":
    write = None
    if "--write-vasp" in sys.argv:
        write = sys.argv[sys.argv.index("--write-vasp") + 1]
    only = None
    if "--only" in sys.argv:
        only = sys.argv[sys.argv.index("--only") + 1]
    results = {}
    for name, rc in RECIPES.items():
        if only and name != only:
            continue
        try:
            results[name] = run_family(name, rc, write)
        except Exception as e:
            print(f"  [ERROR] {type(e).__name__}: {e}")
            results[name] = False
    print("=" * 64)
    for k, v in results.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    sys.exit(0 if all(results.values()) else 1)
