#!/usr/bin/env python3
"""端到端候选枚举：SG → 轨道组合 → 角色 → 化学装饰。

示例：
    python examples/enumerate_candidates.py --sg 180 --max-atoms 18 \\
        --magnetic-elements Mn Fe Cr --nonmagnetic-elements Mg Al Si P S \\
        --max-frameworks 200 --max-candidates 20
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from amgen import enumerate as am_enum


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sg", type=int, help="单空间群；不填则用 --sgs-start/--sgs-end")
    p.add_argument("--sgs-start", type=int, default=168)
    p.add_argument("--sgs-end", type=int, default=194)
    p.add_argument("--target-wave", default="i-wave",
                   help="i-wave / g-wave / all")
    p.add_argument("--max-atoms", type=int, default=24)
    p.add_argument("--max-orbits", type=int, default=6)
    p.add_argument("--magnetic-elements", nargs="+",
                   default=["V", "Cr", "Mn", "Fe", "Co", "Ni"])
    p.add_argument("--nonmagnetic-elements", nargs="+",
                   help="非磁元素池；不填则保留占位元素")
    p.add_argument("--candidate-elements", nargs="+",
                   help="角色枚举时允许作为 M 的元素；默认全部")
    p.add_argument("--min-d", type=float, default=1.6)
    p.add_argument("--min-cn", type=int, help="磁位最小配位数")
    p.add_argument("--charge-neutral", action="store_true",
                   help="开启形式价态中性过滤")
    p.add_argument("--max-magnetic-fraction", type=float,
                   help="磁原子数/非磁原子数上限")
    p.add_argument("--max-nonmag-combos", type=int, default=100)
    p.add_argument("--max-frameworks", type=int, default=2000)
    p.add_argument("--max-candidates", type=int, default=500)
    return p.parse_args()


def main():
    args = parse_args()
    if args.sg is not None:
        sgs = [args.sg]
    else:
        sgs = range(args.sgs_start, args.sgs_end + 1)
    target = None if args.target_wave == "all" else args.target_wave

    n = 0
    for cand in am_enum.enumerate_candidates(
            sgs=sgs, target_wave=target, max_atoms=args.max_atoms,
            max_orbits=args.max_orbits, magnetic_elements=args.magnetic_elements,
            nonmagnetic_elements=args.nonmagnetic_elements,
            candidate_elements=args.candidate_elements, min_d=args.min_d,
            min_cn=args.min_cn, charge_neutral=args.charge_neutral,
            max_magnetic_fraction=args.max_magnetic_fraction,
            max_nonmag_combos=args.max_nonmag_combos,
            max_frameworks=args.max_frameworks,
            max_candidates=args.max_candidates):
        n += 1
        print(f"[{n}] SG{cand.sg} {cand.composition} "
              f"{cand.magnetic_element} {cand.wave} "
              f"orbits={cand.magnetic_orbits} pattern={cand.pattern}")
    print(f"\ncandidates: {n}")


if __name__ == "__main__":
    main()
