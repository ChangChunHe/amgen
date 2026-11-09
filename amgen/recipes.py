"""四个已确认家族的显式配方（回归测试基准 + 枚举空间存在的直接证明）

配方 = SG 号 + {元素: [Wyckoff 字母]}，全部由群分解得到，不依赖任何数据库结构文件。
字母以 spglib/pymatgen 常用原点选择为准；复判与对照对原点选择不敏感
（framework_match 只比较各元素的轨道重数谱）。
"""

BASE = "/share/home/scuthecc/surface_AM"

RECIPES = {
    "Mg6Mn6Ni6": dict(
        sg=180,
        recipe={"Mn": ["f"], "Mg": ["j"], "Ni": ["a", "c"]},
        mag_el="Mn", target_wave="i-wave",
        ref=f"{BASE}/hex-protype/dft2/A2B_hP18_180_fj_ac-001__f6_+++---/Mn__MgMnNi/i/POSCAR",
    ),
    "BaFe2O4": dict(
        sg=182,
        recipe={"Ba": ["b"], "Fe": ["f"], "O": ["g", "c"]},
        mag_el="Fe", target_wave="i-wave",
        ref=f"{BASE}/hex-protype/mp_iwave_scan/dft/BaFe2O4/i/POSCAR",
    ),
    "Cs3Cr2Br9": dict(
        sg=194,
        recipe={"Cs": ["d", "f"], "Cr": ["e"], "Br": ["h", "k"]},
        mag_el="Cr", target_wave="g-wave",
        ref=f"{BASE}/hex-protype/a3m2x9/Cs3Cr2Br9/FM/POSCAR",
    ),
    "Mg2Cr2B2O7": dict(
        sg=188,
        recipe={"Mg": ["e", "f"], "Cr": ["g"], "B": ["h"], "O": ["a", "l"]},
        mag_el="Cr", target_wave="i-wave",
        ref=f"{BASE}/hex-protype/mechanism/bonding_group2/Mg2Cr2B2O7_i/POSCAR",
    ),
}
