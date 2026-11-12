"""amgen.roles — 轨道实例级 M/非磁角色枚举。

现有 ``core.joint_colorings`` 的前提是“哪些轨道有磁性”已经由化学装饰
决定。本模块补上这一层：对每个独立 Wyckoff 轨道实例枚举它是磁性 M 还是
非磁，再在磁性子集上枚举补偿自旋涂色，最后交给 fsg 判定波形。

注意：这里得到的是群论候选，不代表该轨道在真实材料中一定磁化，也不代表
目标波形一定是磁基态；仍需化学装饰、CHGNet 预弛豫和 DFT 竞争能量。
"""
from dataclasses import dataclass
from itertools import combinations, product
import numpy as np

from . import classify, core, relax, wyckoff


@dataclass(frozen=True)
class OrbitInstance:
    orbit_id: int
    letter: str
    indices: tuple
    element: str
    multiplicity: int

    @property
    def label(self):
        return (
            f"{self.element}@{self.letter}^{self.multiplicity}"
            f"#{self.orbit_id}")


@dataclass(frozen=True)
class RoleHit:
    assignment: tuple
    magnetic_orbits: tuple
    mag_indices: tuple
    signs: tuple
    pattern: str
    wave: str
    texture: str
    mag_order: str
    msg: str
    is_altermagnet: bool


def find_orbit_instances(struct, symprec=1e-3):
    """返回轨道实例列表；同一 Wyckoff 字母的复用实例不会被合并。"""
    return [
        OrbitInstance(
            orbit_id=i,
            letter=o["letter"],
            indices=tuple(o["indices"]),
            element=o["element"],
            multiplicity=o["multiplicity"],
        )
        for i, o in enumerate(
            core.find_orbit_instances(struct, symprec=symprec))
    ]


def role_assignments(instances, *, min_magnetic_orbits=1,
                     max_magnetic_orbits=None, min_magnetic_sites=2,
                     max_magnetic_sites=None, min_nonmagnetic_orbits=0,
                     candidate_elements=None, max_orbit_multiplicity=12,
                     max_role_assignments=2000):
    """枚举 M/非磁角色分配。

    只把偶数重、至少 2 个位点的轨道列为候选 M 轨道，因为当前补偿染色
    在每个 M 轨道内等量分配 +/-。candidate_elements 可把候选 M 轨道限制在
    指定元素上，例如只允许 Fe/Cr/Mn 被测试为磁性；max_orbit_multiplicity
    防止 24 重等大轨道直接触发组合爆炸。返回长度为轨道实例数的 bool tuple；
    True 表示该轨道在候选结构中被指定为磁性 M。
    """
    n = len(instances)
    if isinstance(candidate_elements, str):
        candidate_elements = [candidate_elements]
    allowed = None if candidate_elements is None else set(candidate_elements)
    eligible = [
        i for i, o in enumerate(instances)
        if 2 <= o.multiplicity <= max_orbit_multiplicity
        and o.multiplicity % 2 == 0
        and (allowed is None or o.element in allowed)
    ]
    max_k = len(eligible) if max_magnetic_orbits is None else min(
        max_magnetic_orbits, len(eligible))
    count = 0
    for k in range(max(0, min_magnetic_orbits), max_k + 1):
        available_nonmag = n - k
        if available_nonmag < min_nonmagnetic_orbits:
            continue
        for chosen in combinations(eligible, k):
            chosen_set = set(chosen)
            n_sites = sum(instances[i].multiplicity for i in chosen)
            if n_sites < min_magnetic_sites:
                continue
            if max_magnetic_sites is not None and n_sites > max_magnetic_sites:
                continue
            count += 1
            if count > max_role_assignments:
                raise RuntimeError(
                    "角色分配数超过 max_role_assignments="
                    f"{max_role_assignments}；请收紧轨道数/磁轨道数，"
                    "或显式提高上限。")
            yield tuple(i in chosen_set for i in range(n))


def joint_role_colorings(struct, instances, assignment, *,
                         max_colorings=5000):
    """给定一个角色分配，联合枚举所有磁性轨道上的补偿自旋涂色。

    每个 M 轨道内部保持等量 +/-；跨轨道状态由完整结构的空间群操作和
    全局翻号去重。返回 ``(mag_idx, signs, pattern)``。
    """
    if len(assignment) != len(instances):
        raise ValueError("assignment 长度必须等于轨道实例数")
    m_orbits = [o for o, is_mag in zip(instances, assignment) if is_mag]
    if not m_orbits:
        return

    pat_lists = []
    for o in m_orbits:
        if o.multiplicity % 2:
            raise ValueError(
                f"{o.label} 为奇数重轨道；当前补偿染色不支持。")
        coords = [struct[i].frac_coords for i in o.indices]
        pat_lists.append((o, core.orbit_colorings(coords, struct.lattice)))

    total = 1
    for _, pats in pat_lists:
        total *= len(pats)
    if total > max_colorings:
        raise RuntimeError(
            f"该角色分配产生 {total} 个联合染色，超过 max_colorings="
            f"{max_colorings}；请提高上限或先收紧角色分配。")

    m_set = {i for o in m_orbits for i in o.indices}
    n_set = set(range(len(struct))) - m_set
    perms = core.site_perms(
        [site.frac_coords for site in struct],
        struct.lattice,
        species_numbers=[site.specie.Z for site in struct],
    )
    role_perms = [
        p for p in perms
        if set(int(i) for i in p[list(m_set)]) == m_set
        and set(int(i) for i in p[list(n_set)]) == n_set
    ]
    if not role_perms:
        role_perms = [np.arange(len(struct))]

    seen = set()
    for combo in product(*[pats for _, pats in pat_lists]):
        full = np.zeros(len(struct), dtype=int)
        mag_idx, signs, desc = [], [], []
        for (o, _), pat in zip(pat_lists, combo):
            mag_idx.extend(o.indices)
            orbit_signs = [1 if c == "+" else -1 for c in pat]
            signs.extend(orbit_signs)
            full[list(o.indices)] = orbit_signs
            desc.append(f"{o.label}:{pat}")
        key = core.canonical(full, role_perms)
        if key in seen:
            continue
        seen.add(key)
        yield tuple(mag_idx), tuple(signs), "+".join(desc)


def apply_role_species(struct, instances, assignment, magnetic_species,
                       nonmagnetic_species=None):
    """把角色分配落实为元素装饰；非磁轨道不指定时保留原元素。"""
    if len(assignment) != len(instances):
        raise ValueError("assignment 长度必须等于轨道实例数")
    s = struct.copy()
    for o, is_mag in zip(instances, assignment):
        target = magnetic_species if is_mag else nonmagnetic_species
        if target is None:
            continue
        for i in o.indices:
            s.replace(i, target)
    return s


def generate_framework(sg, letters, placeholder="Si", min_d=1.5,
                       return_best=True):
    """按 Wyckoff 字母序列生成抽象几何框架。

    letters 中的重复字母表示同一 Wyckoff 字母被复用的独立轨道实例。
    所有位点先用同一种占位元素生成，便于后续进行不依赖化学的角色枚举。
    return_best=True 时，若没有种子满足 min_d，返回最小间距最大的种子；
    False 时严格失败并抛出。
    """
    letters = list(letters)
    if not letters:
        raise ValueError("letters 不能为空")
    s = wyckoff.generate_ok(sg, {placeholder: letters}, min_d=min_d,
                            return_best=return_best)
    if s is None:
        raise RuntimeError(f"SG{sg} 的框架生成失败")
    return s


def _wave_matches(result, target_wave):
    if target_wave is None:
        return True
    target = target_wave.removesuffix("-wave")
    for key in ("wave", "texture"):
        value = (result.get(key) or "").removesuffix("-wave")
        if value == target:
            return True
    return False


def scan_roles(struct, target_wave="i-wave", *, moment=3.0, symprec=1e-3,
               min_magnetic_orbits=1, max_magnetic_orbits=None,
               min_magnetic_sites=2, max_magnetic_sites=None,
               min_nonmagnetic_orbits=0, candidate_elements=None,
               max_orbit_multiplicity=12, max_role_assignments=2000,
               max_colorings=5000, max_hits=100, resymmetrize_sg=None,
               require_altermagnet=True):
    """在一个框架上自动扫描 M/非磁角色和自旋涂色。

    返回 ``RoleHit`` 迭代器。默认只保留 i-wave；target_wave=None 时保留所有
    fsg 成功分类的涂色。require_altermagnet=True 时只保留 fsg 判定为交错磁
    的命中。candidate_elements 可把候选 M 轨道限制到指定元素；
    max_orbit_multiplicity 限制单轨道大小；symprec 控制轨道识别容差。
    resymmetrize_sg 非 None 时先用 relax.resymmetrize 把结构清理回该母群，
    再跑 fsg，避免弛豫噪声把波形漏判。max_* 是防止组合爆炸的硬上限，超过时
    抛出 RuntimeError；max_hits 为命中上限，达到后停止（截断）。
    """
    if resymmetrize_sg is not None:
        struct, actual = relax.resymmetrize(struct, resymmetrize_sg)
        if struct is None:
            return

    instances = find_orbit_instances(struct, symprec=symprec)
    hits = 0
    for assignment in role_assignments(
            instances,
            min_magnetic_orbits=min_magnetic_orbits,
            max_magnetic_orbits=max_magnetic_orbits,
            min_magnetic_sites=min_magnetic_sites,
            max_magnetic_sites=max_magnetic_sites,
            min_nonmagnetic_orbits=min_nonmagnetic_orbits,
            candidate_elements=candidate_elements,
            max_orbit_multiplicity=max_orbit_multiplicity,
            max_role_assignments=max_role_assignments):
        for mag_idx, signs, pattern in joint_role_colorings(
                struct, instances, assignment, max_colorings=max_colorings):
            result = classify.fsg_classify(
                struct, mag_idx, signs, moment=moment)
            if require_altermagnet and not result["is_altermagnet"]:
                continue
            if not _wave_matches(result, target_wave):
                continue
            if hits >= max_hits:
                return
            hits += 1
            yield RoleHit(
                assignment=assignment,
                magnetic_orbits=tuple(
                    o.label for o, is_mag in zip(instances, assignment) if is_mag),
                mag_indices=mag_idx,
                signs=signs,
                pattern=pattern,
                wave=result["wave"],
                texture=result["texture"],
                mag_order=result["mag_order"],
                msg=result["msg"],
                is_altermagnet=result["is_altermagnet"],
            )
