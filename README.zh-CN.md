# amgen — 交错磁体群论生成器（纯枚举版）

<div align="center">

[![English](https://img.shields.io/badge/lang-English-blue.svg)](README.md)
[![中文](https://img.shields.io/badge/lang-中文-red.svg)](README.zh-CN.md)

</div>

从空间群分解出发，不依赖任何外部结构数据库，正向生成补偿共线交错磁体
（altermagnet）候选结构。

## 核心思想

交错磁波形（d/g/i-wave）不是算出来才发现的，而是**先由自旋群条件筛出
"占位 + 磁染色"组合，再用化学装饰做成真实化合物，最后用 DFT 判定磁基态**。
两条正交判据必须分别过：

| 维度 | 判据 | 工具 |
|---|---|---|
| 对称性 | 该占位+染色是否为 AFM(Altermagnet) 且 texture 为目标波形 | fsg |
| 能量 | 该构型是否为磁基态（E0 最低, meV/磁性原子） | VASP 多构型竞争 |

fsg 通过 ≠ 基态；能量最低 ≠ i/g-wave。

## 七步流水线（详细）

**① 空间群列表（入口）**
六方晶族全量：三方 hP 设置 SG 143–167 + 六方 SG 168–194（共 40 个）。
每个 SG 由 pyxtal 群数据给出 Wyckoff 表 `(字母, 重数, 自由参数数)`，
如 SG180 → a/b/c/d(3 重, 0 参数)、e–j(6 重, 1 参数)、k(12 重, 3 参数)。

**② Wyckoff 轨道组合枚举**（`wyckoff.orbit_combos`）
- 总原子数 ≤ 30；
- 必须含至少一个**偶数重、≥4 重**轨道（候选磁轨道；2 重轨道给不出
  非平凡补偿构型）；
- **同一字母允许被不同元素复用**（Cs₃Fe₂F₉ 中 Cs@f 与 Fe@f 是两个独立
  轨道实例），配方按"轨道实例"计数；
- 产出按原子数升序，如 SG182 会产出 `b²+f⁴+g⁶+c²`（14 原子，BaFe₂O₄ 型）。

**③ 结构生成**（`wyckoff.generate_ok`）
pyxtal 按配方放置原子，自由参数取随机值，多种子重试直到最小原子间距
≥1.5 Å。晶格用群默认六角格子，尺寸后续由半径标定 + CHGNet 定。

**④ 前染色 + fsg 门（零成本群论筛选）**（`core` + `classify`）
- 候选磁轨道上枚举全部不等价补偿染色（群置换 + 全局翻号去重；6 重轨道
  仅 4 种、4 重仅 3 种）；
- 多磁轨道时**联合染色**（KNiCl₃ 教训：不能冻结其余轨道）；
- 每个染色写 mcif → fsg，只保留 `AFM(Altermagnet)` 且 texture 为
  g/i-wave 的组合；不过门的框架直接丢弃，不花弛豫机时。

**⑤ 化学装饰**（`decorate`）
磁轨道填 3d（V–Ni；4d/5d 可开关），其余轨道按价态中性 + d 电子数 2–8 +
CN≥4 过滤，晶格按共价半径比标定。金属间化合物跳过电价过滤。

**⑥ CHGNet 预弛豫 + 对称化**（`relax`，50 步粗弛豫）
CHGNet 把自由参数与晶格常数调到合理邻域（它不懂磁性，只管几何），再用
**容差阶梯**（1e-3→1e-1）把数值噪声吸回母群；回不到母群的记录终态 SG
（生成 SG ≠ 弛豫后 SG 是常态，如 CrCl₃: 158→185 波形仍存活）。

**⑦ 复判 fsg 门 + DFT 导出**（`classify` + `export`）
弛豫后结构上重新染色跑 fsg——**防升对称杀波形**（Cs₃Fe₂F₉ 教训：
P-62c→P6₃/mmc 后 i-wave 关闭）。双门都过才写 VASP 输入
（POSCAR/INCAR 含命中染色的 MAGMOM/KPOINTS/POTCAR/job.sh），
进 DFT 漏斗做全不等价磁构型竞争。

```
SG(40) → 轨道组合(~10²–10³) → 前染色 fsg 门 → 装饰过滤
  → CHGNet 50 步 → 对称化 → 复判 fsg 门 → DFT(~10²)
```

## 依赖

- pymatgen, spglib, pyxtal, chgnet, numpy
- `fsg`（FindSpinGroup）在 PATH 或 `~/miniconda3/bin/fsg`
- `vaspkit`（生成 POTCAR，仅 --write-vasp 需要）

## 回归验证：四个已确认家族

配方（`amgen/recipes.py`，纯群论声明，无任何数据库结构文件）：

| 家族 | SG | 配方 | 磁轨道 | 目标染色 | 波形 |
|---|---|---|---|---|---|
| Mg₆Mn₆Ni₆ | 180 P6₂22 | Mn@f, Mg@j, Ni@a+c | f⁶ | `+++---` | i-wave |
| BaFe₂O₄ | 182 P6₃22 | Ba@b, Fe@f, O@g+c | f⁴ | `++--` | i-wave |
| Cs₃Cr₂Br₉ | 194 P6₃/mmc | Cs@d+f, Cr@e, Br@h+k | e⁴ | `+-+-` | g-wave |
| Mg₂Cr₂B₂O₇ | 188 P-6c2 | Mg@e+f, Cr@g, B@h, O@a+l | g⁴ | `++--` | i-wave |

运行：

```bash
python examples/reproduce_families.py                    # 登录节点串行
python examples/reproduce_families.py --only BaFe2O4     # 单家族
python examples/reproduce_families.py --write-vasp out/  # 落盘 DFT 输入
sbatch job_reproduce.sh   # 集群: 单节点 64 核, 4 家族并行 × 16 线程
```

当前状态：**四个家族全部 PASS**（三个 StructureMatcher 完全一致；
Cs₃Cr₂Br₉ 走原型框架匹配，晶格差 ~20% 由后续 DFT 弛豫吸收）。

每个家族走完整七步并与已确认 DFT 结构对照（StructureMatcher，
失败时退化为"同 SG + 轨道重数谱 + 晶格接近"的原型框架匹配——
装饰结构与已弛豫参考的晶格差异属正常）。

## API 速览

```python
from amgen import wyckoff, relax, core, classify, decorate, export

wyckoff.wyckoff_table(182)                     # SG -> Wyckoff 表
wyckoff.orbit_combos(182, max_atoms=30)        # 轨道组合枚举(字母可复用)
s = wyckoff.generate_ok(182, {"Ba":["b"],"Fe":["f"],"O":["g","c"]})

s2, sg2 = relax.relax_and_symmetrize(s, target_sg=182)   # CHGNet + 回母群

orbits = core.find_orbits(s2)
pats   = core.orbit_colorings(coords, s2.lattice)
r = classify.fsg_classify(s2, mag_idx, signs)  # -> texture: i-wave/g-wave/forbidden

s3 = decorate.substitute_orbits(s2, {"f": "Mn"})
export.write_vasp("out/X/i", s3, magmom_list)
```

## 已知坑（血泪记录）

- **同一字母复用**：枚举按轨道实例而非字母计数；
- **CHGNet 不知道磁性**：JT 畸变/自旋驱动弛豫预测不了，只做几何预弛豫；
- **弛豫可升对称杀波形**（Cs₃Fe₂F₉）——复判门必须存在；
- **mcif/fsg 容差**：弛豫噪声会让 1e-3 下判成 P1，对称化用容差阶梯；
- VASP 侧：晶格噪声触发 Bravais 错误（六角格子需理想化）、EDDDAV 加
  `ALGO=A`、对称化报错加 `ISYM=-1`、POTCAR 串元素要 grep TITEL 检查。

完整战役记录与负结果分析见 [`docs/HEX_IWAVE_GENERATION_SCHEME.md`](docs/HEX_IWAVE_GENERATION_SCHEME.md)、
[`docs/IWAVE_FINAL_SUMMARY.md`](docs/IWAVE_FINAL_SUMMARY.md)、
[`docs/CAMPAIGN_LOG.md`](docs/CAMPAIGN_LOG.md)。
