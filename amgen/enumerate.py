"""amgen.enumerate — 端到端 i/g-wave 候选枚举。

把分散的步骤连起来：
SG → Wyckoff 轨道组合 → 生成抽象框架 → 自动角色枚举 → 化学装饰 → 过滤。
"""
from dataclasses import dataclass

from . import decorate, roles, wyckoff


@dataclass(frozen=True)
class Candidate:
    sg: int
    letters: tuple
    structure: object
    composition: str
    magnetic_element: str
    nonmagnetic_map: tuple
    magnetic_orbits: tuple
    pattern: str
    wave: str


def enumerate_candidates(sgs=None, *, target_wave="i-wave",
                         max_atoms=30, max_orbits=6,
                         magnetic_elements=None,
                         nonmagnetic_elements=None,
                         candidate_elements=None,
                         min_d=1.2, min_cn=None, charge_neutral=False,
                         max_magnetic_fraction=None,
                         max_nonmag_combos=100,
                         max_frameworks=2000, max_candidates=500,
                         role_scan_kwargs=None):
    """枚举全部 SG/轨道组合下的目标波形候选。

    sgs 默认取六方/三方 hP 设置；magnetic_elements 默认 3d 磁芯
    {V,Cr,Mn,Fe,Co,Ni}。nonmagnetic_elements 为 None 时非磁轨道保留占位
    元素（只做对称候选）；给定元素池时对非磁轨道做笛卡尔积并跑化学过滤。
    max_* 是硬上限，max_candidates 为命中上限（达到后停止）。
    """
    if sgs is None:
        sgs = range(168, 195)  # 六方晶族；三方 hP 可自行传 range(143, 168)
    if magnetic_elements is None:
        magnetic_elements = decorate.MAG_SET_3D
    role_scan_kwargs = role_scan_kwargs or {}

    n_frameworks = 0
    n_candidates = 0
    for sg in sgs:
        combos = wyckoff.orbit_combos(sg, max_atoms=max_atoms,
                                      max_orbits=max_orbits)
        for combo in combos:
            letters = tuple(letter for letter, _ in combo)
            try:
                framework = roles.generate_framework(sg, letters, min_d=min_d)
            except Exception:
                continue
            if framework is None:
                continue
            if wyckoff.min_distance(framework) < min_d:
                continue
            n_frameworks += 1
            if n_frameworks > max_frameworks:
                raise RuntimeError(
                    f"框架数超过 max_frameworks={max_frameworks}；请收紧 SG/"
                    "max_atoms/max_orbits。")

            for hit in roles.scan_roles(
                    framework, target_wave=target_wave,
                    candidate_elements=candidate_elements,
                    resymmetrize_sg=sg, **role_scan_kwargs):
                if n_candidates >= max_candidates:
                    return
                nonmag_slots = [
                    i for i, is_mag in enumerate(hit.assignment) if not is_mag]
                try:
                    if nonmagnetic_elements is None:
                        fills = [None]
                    else:
                        fills = list(decorate.iter_nonmagnetic_assignments(
                            len(nonmag_slots), nonmagnetic_elements,
                            max_combos=max_nonmag_combos))
                except RuntimeError:
                    continue
                for fill in fills:
                    nonmag_map = (None if fill is None else
                                  dict(zip(nonmag_slots, fill)))
                    for mag_el in magnetic_elements:
                        try:
                            s = decorate.decorate_assignment(
                                framework, hit.assignment, mag_el, nonmag_map)
                        except Exception:
                            continue
                        if not decorate.decorated_ok(
                                s, hit.mag_indices, mag_el, min_cn=min_cn,
                                charge_neutral=charge_neutral,
                                max_magnetic_fraction=max_magnetic_fraction):
                            continue
                        n_candidates += 1
                        yield Candidate(
                            sg=sg,
                            letters=letters,
                            structure=s,
                            composition=s.composition.reduced_formula,
                            magnetic_element=mag_el,
                            nonmagnetic_map=tuple(
                                sorted(nonmag_map.items())
                                if nonmag_map is not None else ()),
                            magnetic_orbits=hit.magnetic_orbits,
                            pattern=hit.pattern,
                            wave=hit.wave,
                        )
                        if n_candidates >= max_candidates:
                            return
