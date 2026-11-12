#!/usr/bin/env python3
"""自动枚举 M/非磁轨道角色，并筛选目标波形的自旋涂色。

示例：
    python examples/enumerate_roles.py POSCAR --target-wave i-wave
    python examples/enumerate_roles.py --sg 180 --orbits f f j --target-wave i-wave

给定 --sg/--orbits 时，先用同一种占位元素生成抽象框架；给定结构文件时，
直接在该结构的轨道实例上扫描。扫描不改变元素，只把自旋画到候选 M 轨道上。
"""
import argparse
import sys
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from amgen import roles


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("structure", nargs="?",
                   help="POSCAR/CIF 等结构文件；与 --sg/--orbits 二选一")
    p.add_argument("--sg", type=int, help="空间群号，用于生成抽象框架")
    p.add_argument("--orbits", nargs="+",
                   help="Wyckoff 字母序列，如 f f j；重复表示复用轨道实例")
    p.add_argument("--target-wave", default="i-wave",
                   help="目标波形，如 i-wave/g-wave；all 表示不过滤")
    p.add_argument("--moment", type=float, default=3.0, help="fsg 用的磁矩")
    p.add_argument("--symprec", type=float, default=1e-3,
                   help="空间群识别容差；弛豫后的结构可放宽到 0.01–0.1")
    p.add_argument("--resymmetrize-sg", type=int,
                   help="跑 fsg 前先对称化回该空间群（消除弛豫噪声）")
    p.add_argument("--no-require-altermagnet", action="store_true",
                   help="不过滤 is_altermagnet（默认要求是交错磁）")
    p.add_argument("--candidate-elements", nargs="+",
                   help="只允许这些元素所在的轨道作为 M，如 Fe Cr Mn")
    p.add_argument("--max-orbit-multiplicity", type=int, default=12,
                   help="单个候选磁轨道的最大重数")
    p.add_argument("--max-role-assignments", type=int, default=2000)
    p.add_argument("--max-colorings", type=int, default=5000)
    p.add_argument("--max-hits", type=int, default=50)
    return p.parse_args()


def load_structure(args):
    if args.structure:
        return Structure.from_file(args.structure)
    if args.sg is None or not args.orbits:
        raise SystemExit("需要结构文件，或同时给出 --sg 与 --orbits")
    return roles.generate_framework(args.sg, args.orbits)


def main():
    args = parse_args()
    struct = load_structure(args)
    sg = SpacegroupAnalyzer(struct, symprec=args.symprec).get_space_group_number()
    print(f"framework: {len(struct)} sites, SG{sg}, "
          f"{struct.composition.reduced_formula}")

    instances = roles.find_orbit_instances(struct, symprec=args.symprec)
    print("orbit instances:")
    for o in instances:
        print(f"  {o.label:>16s}  n={o.multiplicity:<2d}  "
              f"indices={o.indices}")

    target = None if args.target_wave == "all" else args.target_wave
    n_hits = 0
    for hit in roles.scan_roles(
            struct, target_wave=target, moment=args.moment,
            symprec=args.symprec,
            candidate_elements=args.candidate_elements,
            max_orbit_multiplicity=args.max_orbit_multiplicity,
            max_role_assignments=args.max_role_assignments,
            max_colorings=args.max_colorings, max_hits=args.max_hits,
            resymmetrize_sg=args.resymmetrize_sg,
            require_altermagnet=not args.no_require_altermagnet):
        n_hits += 1
        print(f"\n[hit {n_hits}] M = {', '.join(hit.magnetic_orbits)}")
        print(f"  pattern : {hit.pattern}")
        print(f"  wave    : {hit.wave}  ({hit.mag_order})")
        print(f"  fsg     : {hit.msg}")

    print(f"\nhits: {n_hits} (target={args.target_wave})")


if __name__ == "__main__":
    main()
