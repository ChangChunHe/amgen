# 基于 hex-protype 的 i 波交错磁体生成方案（详细版 v1.0）

> 定位：本文件是**可执行的方法方案**（从原型池到 DFT 判决的完整流水线），不是文章正文。
> 论文体例的表述见 `SUPPLEMENTARY_METHODS_iwave.md`；战役流水账见 `CAMPAIGN_LOG.md`；结果花名册见 `IWAVE_FINAL_SUMMARY.md`。
> 数据口径：VASP-PBE / PAW，ENCUT = 520 eV，EDIFF = 1e-6；能量一律 E0(σ→0)；波形一律 fsg（FindSpinGroup）对 mcif 的 leading spin texture。
> 整理日期：2026-10-24。

---

## 0. 一页摘要

**目标**：在六角/三方框架里，系统生成"无 SOC 自旋纹理为 6 阶（i 波）"的补偿共线磁体，并用 DFT 判定它确实是磁基态。

**核心思想**：i 波不是"算出来才发现的"，而是**先由自旋群条件筛出可行的占位 + 磁染色组合，再用化学替换把它做成真实化合物，最后用 DFT 验证磁基态**。群论给条件命题，DFT 验证前件。

**三层漏斗**：

```
对称性层（零 DFT）：AFLOW 六角原型 168–194（248 个）
      → Wyckoff 轨道 + 补偿磁染色 → fsg 判 order-6        [i 波规则]
化学层（零 DFT）：3d 磁芯替换 + 价态/配位/半径/电荷过滤     [可合成性]
能量层（DFT）：全不等价磁构型竞争（E0，meV/磁性原子）      [磁基态]
      → 基态 fsg 复判 → 稳定性/表征（声子/凸包/MAE/能带）
```

**本战役实际漏斗**（详见 §10）：
248 原型 → 63 个"3d + 偶数重轨道"原型 → 单轨道 1223 染色（18 个 i 波）
→ 多轨道联合 8649 构型（238 个 i 波候选）→ dft3 238 个候选弛豫（52 收敛、12 过门）
→ 全构型竞争 → 纯 3d 幸存 0；靠"替换 + 5d 框架解禁 + 装饰伪影复活"另开四线，
最终确认 **14 个 i 波磁基态 + 1 个暂定**，另收获 **2 个 g 波基态**。

**一句话结论**：i 波在六角原型里**存在但极稀有**；能不能活到基态，几乎完全由"化学"（磁芯 d 填充、A 位尺寸、框架电负性）决定，而不是由几何单独决定。

---

## 1. 方案总纲：两个正交维度

生成 i 波必须同时满足两个独立条件，缺一不可：

| 维度 | 问题 | 判据 | 工具 |
|---|---|---|---|
| **对称性** | 这个"占位 + 磁染色"在群论上是不是 i 波？ | fsg：`AFM(Altermagnet)` + `order 6` | fsg、spglib |
| **能量** | 这个磁构型是不是该化合物的磁基态？ | E0 最低（meV/磁性原子） | VASP 多构型竞争 |

因此流水线分成两条并行支线：

- **生成支线**（§3–§6）：原型 → 染色 → 装饰 → DFT 输入；
- **判决支线**（§7–§9）：i 弛豫 → fsg 门控 → 竞争构型 → 能量判决 → 稳定性。

**最容易犯的错误**是把两条支线混为一谈：fsg 通过 ≠ 基态；能量最低 ≠ i 波。两个判据必须分别过。

---

## 2. 判据：什么算 i 波

### 2.1 角谐波定义

共线 z 向磁矩、无 SOC 时，自旋劈裂

$$
\Delta(\mathbf{k}) = \varepsilon_\uparrow(\mathbf{k}) - \varepsilon_\downarrow(\mathbf{k})
$$

在面内极角 φ 下展开为 $\Delta = \sum_l A_l \cos(l\varphi + \varphi_0)$。leading 阶数 l 决定波型：

| 波型 | l | 角依赖 | 节点 |
|---|---|---|---|
| s | 0 | 常数 | 无（铁磁） |
| d | 2 | cos2φ | 2 条节点线 |
| g | 4 | cos4φ | 4 条（CrSb/MnTe 型） |
| **i** | **6** | **cos6φ** | **6 条 / 12 花瓣** |

本方案只追 l=6；d/g 是竞争构型的副产物（同一框架内常见）。

### 2.2 自旋群四条硬条件（i 波的必要条件）

对 cos6φ σz 型序参量，磁（自旋空间）群必须：

1. **保留 unprimed 的 C3 或 C6**：C3 保证 l 是 3 的倍数，C6 进一步排除 l=3，使 leading 落到 l=6；
2. **不允许 primed 转动或 primed 镜面**：primed C6 = (−1|6) 会把 cos6φ 变号，精确禁戒 i 波；
3. **不允许纯反平移** (−1|t)：否则所有 k 点自旋简并（退化为普通 AFM）；
4. **自旋翻转只能挂在带平移的算符上**：3D 中是 6₃ 螺旋 + c-glide；2D 中是层群 g-glide。

已确认的体相 i 波全部是 **type-I MSG**（所有空间群操作 unprimed），OSSG 形如

```
G_s = P 1 | 6_3   −1 | c   −1 | m   ∞_{001} m | 1
```

### 2.3 为什么母体池必须是 C3/C6 框架

C4 不变性要求 f(k)=f(C4⁻¹k)，90° 扇区内 cos6φ 无法闭合 → **正方/四方框架无 i 波**。单镜面单斜同理。故母体池限定在六角 168–194（三方 143–167 可作 g 波，见 tri-protype 支线）。

### 2.4 fsg 判据（入口与终判统一）

```bash
fsg <structure>.mcif
```

必须同时看到：

```
Magnetic order: Collinear; AFM(Altermagnet)
Leading spin texture: without SOC i-wave (order 6, collinear); ...
Net moment: 0 μB
```

工程要点：

- fsg 的 allowed/forbidden 是**对称性陈述，不是劈裂量级**；量级要靠能带。
- `Leading spin texture` 行可能输出到 **stderr**，必须 stdout+stderr 合并解析。
- 正则捕获整句 `without SOC (.+?)(?:\s*\(order|;|$)`；用 `[igdp]` 字符类会漏 f-wave/forbidden。
- mcif 文件名不能以 `.` 开头（fsg 不识别隐藏文件）。
- 磁矩全塌（<0.02 μB/原子）时 fsg 报 `degenerates to nonmagnetic` → 记"无法判定"，不是"非交错磁"。

---

## 3. 阶段 A：母体池（AFLOW 六角原型）

### 3.1 抓取与校验

1. 从 AFLOW Prototype Encyclopedia 的 hexagonal 页面抓空间群 **168–194** 的全部原型（27 个群），去重 **248 个**；
2. 每个原型下载 `aflow.cif`（redesignate 的旧 POSCAR 链接已失效，用 CIF）；
3. pymatgen 读入 → POSCAR，按 `hexagonal_prototypes/poscar/sg<SG>/` 落盘；
4. spglib（symprec=1e-3）校验：文件内 SG = 标签 SG，Pearson 原子数 = Wyckoff 多重度之和；
5. 部分占位（mixed occupancy）原型按主成分有序化后重校验。

**校验结果**：248/248 转换成功；239 个直接通过；9 个部分占位有序化后吻合。

### 3.2 可枚举子集

| 条件 | 数量 |
|---|---|
| 全部六角原型 168–194 | 248 |
| 含 3d 磁性元素（V–Ni） | 66 |
| 同时具备"3d 元素 + 偶数重（2–12）单位轨道" | **63**（主枚举对象） |

### 3.3 原型池的定位

- 它是"**已被实验验证可存在的高对称占位模板**"，不是最终可计算结构；
- 自由参数（a、c/a、内坐标 x）要给定初值，最终几何交给 DFT（ISIF=3）；
- 原型不编码磁性，磁序必须由染色枚举决定；
- 原型库是**代表子集**，不穷尽某空间群的全部 Wyckoff 组合；要穷尽需 pyXtal/ISODISTORT 级枚举。

---

## 4. 阶段 B：对称性染色枚举

### 4.1 提取轨道与位点置换群

对每个原型：

1. `SpacegroupAnalyzer(symprec=1e-3).get_symmetrized_structure()` 取 Wyckoff letter 分组；
2. 只保留"元素 ∈ 磁性可承载集、多重度为偶数、2 ≤ N ≤ 12"的轨道；
3. 用 spglib（symprec=1e-4）取空间群操作在磁位点上诱导的**置换群** P：对每个 {R|t}，检查是否把磁位点集合映到自身，记录置换。

磁性可承载集（早期全量）：`{V,Cr,Mn,Fe,Co,Ni,Cu,Nb,Mo,Ru,Rh,Pd,Ag,Ta,W,Re,Os,Ir,Pt}`。
**后期硬规则**：磁芯只允许 3d `{V,Cr,Mn,Fe,Co,Ni}`；4d/5d 只能做非磁框架；双磁轨道装饰禁用；Tb 不计。

### 4.2 canonical 平衡染色

N 重轨道上，补偿染色 = N/2 个 ↑、N/2 个 ↓。在置换群 P × 全局翻转 Z2 下取规范代表：

$$
[s] = \min_{\pi\in P,\ \sigma=\pm1} \sigma\, s_{\pi}
$$

只枚举互不等价的代表。N=6 时恰好 10 个（FM 1 + 5:1 型 1 + 4:2 型 3 + 3:3 型 5），其中 i 波只占 1 个。

### 4.3 单轨道扫描（第一层）

`color_scan.py`：逐轨道独立染色，其余轨道冻结；其余非磁原子保留真实元素、零磁矩；fsg 在**完整原型结构**上判 texture。

| 项 | 数量 |
|---|---|
| 枚举组合（63 原型 × 偶数重轨道） | **1223** |
| i-wave 命中 | **18** |
| 其中 d / g / s / forbidden | 500 / 100 / 383 / 222 |

`color_scan_all.py`（放宽倍率上限/轨道范围）得到 **4248 组合、77 个 i 波**，作为交叉核对。

> 这一步的意义是**给出 i 波的"存在性证书"**：i 波在哪些 SG、哪些 Wyckoff 轨道、哪些 ± 图案上可实现。它不涉及化学。

### 4.4 多轨道联合染色（最终方法）

**动机（KNiCl₃ 教训）**：磁性元素若占据多个 Wyckoff 轨道，只给一个轨道染色、其余置零，会构造出"人工单轨道子空间"——该子空间对称性自洽、fsg 报 i 波，但真实基态在未计算的多轨道磁空间里。KNiCl₃ 的单轨道 i 构型看着是基态，双轨道 DFT 证明真实基态是 **g 波**（a=+−, b=++−−，低 7.4 meV）。

`color_scan2_new.py` 的做法：

1. 对每个**装饰后的候选**，找出全部磁性可承载轨道；
2. 每个偶数重轨道枚举全部 canonical 平衡染色；
3. 跨轨道做笛卡尔积 = 联合磁构型；
4. 组合数上限 300，超过则记 `n_combos` 跳过；
5. 奇数重轨道置零并在输出中标记（亚铁磁风险）。

| 项 | 数量 |
|---|---|
| 参与联合枚举的装饰候选 | **251** |
| 联合磁构型 | **8649** |
| i-wave 构型 / 含 i 波候选 | **238 / 238** |

i 波候选按空间群分布（联合枚举）：

| SG | 177 | 180 | 186 | 188 | 189 | 190 | 194 |
|---|---|---|---|---|---|---|---|
| 候选数 | 8 | 16 | 23 | 16 | 33 | 126 | 16 |

### 4.5 已确立的 i 波规则（体相）

| 类别 | 空间群 | 轨道 / 图案 | 典型 MSG |
|---|---|---|---|
| 手性 type-1 | 178 P6₁22 | 12c `++++++------` | 178.155 t1 |
| 手性 type-1 | 179 P6₅22 | 12c `++++++------` | 179.161 t1 |
| 手性 type-1 | 180 P6₂22 | 12k `++++++------`；6f `+++---` | 180.167 t1 |
| 手性 type-1 | 181 P6₄22 | 12k `++++++------`；6f `+++---` | 181.173 t1 |
| 多轨道 | 185/186/188/189/190/191/193/194 | 多轨道联合设置（2/4/6/12 重任意组合） | P6₃cm 型 type-I |

经验边界：

- 单轨道 ≤12 重框架下，**三角 143–177 与六角 185–194 的单轨道 i 波为零**（三角支线 2363 染色实测 0 个 i，见 tri-protype）；185–194 的 i 波全部落在多轨道/退化设置；
- 178–181 靠单轨道 6f/12k 就能出 i，是"干净的实现"；
- 最低门槛不是"重数越大越好"，而是"轨道上恰好存在 unprimed 6₃ + 反平移"。

---

## 5. 阶段 C：化学装饰

### 5.1 磁芯与初始磁矩

| M | V | Cr | Mn | Fe | Co | Ni |
|---|---|---|---|---|---|---|
| μ (μB) | 2.0 | 3.0 | 3.5 | 3.0 | 2.0 | 1.5 |

（4d/5d 作为非磁框架时不写磁矩；后期 4d¹⁰ 框架如 Cd²⁺ 也可用。）

### 5.2 五道过滤

对每个 i 波命中轨道，按顺序施加：

1. **配位数 CN**：以候选磁原子为心，数 3.2 Å 内非同名原子邻居，要求 **CN ≥ 4**（现行 2D 用 ≥4，体相主池历史用 ≥4）。
2. **电荷中性**：非着色原子固定电荷和 Q_fixed，磁轨道原子数 n_M，则
   $$
   \mathrm{os}_M \approx -Q_\mathrm{fixed}/n_M
   $$
   在 ±0.51 容差内选允许氧化态 os_M ∈ OS(M)，并要求 d 电子数 2 ≤ n_d = Z_M − os_M ≤ 8。
   金属间化合物（无阴离子）跳过电价过滤；Zintl/金属间区间可放行 `need ∈ [−4, +1.5]`（见 §5.3）。
3. **晶格标定**：按共价半径比缩放体积
   $$
   V' = V\left(\frac{r_\mathrm{cov}(M)}{r_\mathrm{cov}(M_\mathrm{native})}\right)^3
   $$
   让 M–配体键长初值合理；最终由 ISIF=3 弛豫。
4. **化学式去重**：按 (SG, 约化化学式) 去重。
5. **fsg 入口复判**：在**完整结构（磁 + 配体）**上跑 fsg，只放行 i-wave；配体导致的对称性破缺在此剔除。

### 5.3 两个装饰伪影与修复（bucket C 复活经验）

| 伪影 | 现象 | 修复 |
|---|---|---|
| CN 判定顺序错 | Ca/Na 长键框架先按 3.2 Å 数 CN 被误杀 | **先按共价半径比标定晶格，再数 CN** |
| 电荷簿记过严 | Cu₃P 系按 Cu:+2 算 `need` 超界被误杀 | 放行 Zintl/金属间区间 `need ∈ [−4, +1.5]`；Cu 按 +1 试 |

修正后复活 18 个候选，直接产出 SG185 b4 三个 i 基态（§11）。

### 5.4 化学禁令（现行）

- **双磁性装饰禁用**：两个磁轨道同时给磁矩 → 451 个候选实测无 viable（收敛 81 个，基态全为非磁畸变或 NM）；
- **4d/5d 不做磁芯**：Re/Ru 系 i 波纹理成立但 FM 低 0.85/2.26 eV；XO₂ 六角 n12 六种全部 NM；
- **允许 3d 磁芯 + 4d/5d 非磁框架**：这是头号富矿（SG189 IrO₄/OsO₄ 家族）；
- Tb 不计（f 电子自旋群处理复杂、POTCAR/磁矩解析易错）。

### 5.5 候选输出

```
candidates_new/<prototype_label>__<wyckoff><mult>_<pattern>/<M>__<formula>/
├── POSCAR
├── MAGMOM.FM        # 全同号
├── MAGMOM.AM        # 目标 i 波染色（符号表来自 §4）
├── MAGMOM.AFM-G     # 按 M–M 接触图 BFS 二色（3.2 Å 截断）
├── meta.json        # 轨道、图案、元素、价态、CN、缩放比、fsg 复判
└── ...
```

| 池 | 候选数 |
|---|---|
| 早期单轨道装饰 `candidates/` | 119 |
| 多轨道联合装饰 `candidates_new/` | 267（其中 251 有联合构型、238 含 i 波） |
| 伪影复活 `candidates_revive/` | 18 |

---

## 6. 阶段 D：DFT 输入与计算协议

### 6.1 标准 INCAR（体相弛豫）

```
EDIFF   = 1e-06
EDIFFG  = -0.01
ENCUT   = 520
IBRION  = 2
ICHARG  = 2
ISIF    = 3
ISMEAR  = 0
ISPIN   = 2
ISTART  = 0
LCHARG  = F
LREAL   = Auto
LWAVE   = F
NCORE   = 4
NSW     = 100
LORBIT  = 11
MAGMOM  = <逐原子初值>
```

要点：

- **LORBIT=11 必开**，否则 OUTCAR 没有磁矩块，后面 mcif 生成直接失败；
- **NCORE=4**（部分节点 `NPAR=4` 会 segfault）；
- 磁构型计算 **ISYM=−1**（关闭对称性，避免把磁序劈裂对称化掉）；
- 小胞若 SCF 极慢，`LREAL=.FALSE.`（VASP 对 8–30 原子小胞的 Auto 投影效率差）；
- k 点密度 **3000/atom**（等价 KSPACING），或直接用 vaspkit 0.04；
- 所有能量比较必须在**同一二进制、同一 ENCUT、同一 ISYM/LREAL** 下进行；换版本 = 旧数据作废。

### 6.2 生成方式

```bash
pyvasp prep_single_vasp POSCAR -a job_name=<name>,ENCUT=520,NSW=100,ISIF=3,ISPIN=2
# 若 POTCAR 映射失败，追加 sym_potcar_map=K_sv,Ca_sv,...
```

随后统一覆盖：`EDIFF=1e-6`、`NPAR→NCORE=4`、追加 `LORBIT=11` 与 `MAGMOM`。

### 6.3 提交纪律

```bash
sbatch --chdir=<绝对路径> <绝对路径>/job_vasp.sh
```

- `--chdir` 必须绝对路径（相对路径会导致 WorkDir 翻倍、作业秒败）；
- 单次 ≤ 50 作业（大量作业用打包脚本串跑）；
- 批量脚本必须幂等（已有 OSZICAR 的目录跳过）；
- 新设置先 1–2 个冒烟，再放量。

### 6.4 故障处方表

| 错误签名 | 处方 |
|---|---|
| `ZBRENT`（电子步） | 加阻尼 `AMIX=0.2, BMIX=1e-4, NELMDL=5`；或 MAGMOM 减半 |
| `ZBRENT: can not reach required accuracy`（离子步） | `cp CONTCAR POSCAR` + `POTIM=0.3` 续算 |
| `ZHEGV / EDDDAV` | `ISYM=-1` + `ALGO=A` |
| `ZPOTRF`（原子过近） | 同上，降 POTIM |
| `Try calculating without symmetrizing` | `ISYM=-1`（或 ISYM=0）重算 |
| `NELM` 不够 | `NELM=200` |
| signal 9 / 节点闪退 | 换节点重交，不改物理设置 |
| 0 字节 CHGCAR/CONTCAR | 清理后重交；ICHARG=11 前必须确认 CHGCAR 完整 |

---

## 7. 阶段 E：磁基态判决（本方案的核心）

### 7.1 构型空间要算满

对 N 个等价磁位点，全部 ± 构型在"空间群置换群 P × 全局翻转"下约化为少数代表：

| N | 不等价构型数 | 结构 |
|---|---|---|
| 2 | 2 | FM、AFM（`+-`；在 2 位点上 i 波与 AFM 同型） |
| 4 | 4 | FM、FiM(3:1)、AFM(d2)、AFM(g2/`+-+−`) |
| 6 | **10** | FM 1 + 5:1 FiM 1 + 4:2 FiM 3 + 3:3 补偿 5（i、gd、d2、d3、G） |
| 12 | 更多 | 建议用程序枚举 + 能量预筛 |

**硬教训**：只算 i/gd/FM 三件套会漏掉 d2/d3/F42 这些 3:3 竞争者。多个曾判"i 基态"的候选在补算全构型后被 d3/d2/FiM 推翻（wave12 六候选全灭；subst 的 Ti₆Mn₆Ni₆ 被 d3 低 18.1 meV、Mg₆Mn₆Cu₆ 被 d3 低 0.8 meV）。

**等价性检验**：`tmp/magcfg_equiv.py` 从弛豫结构提取置换群，对符号向量做 canonical 化，相同即对称等价，只算一个。

### 7.2 起始几何协议

1. 先从原型/母相几何**只弛豫 i 构型**；
2. i 收敛后，用 **i 的 CONTCAR 作为所有竞争构型的初始 POSCAR** 再弛豫（各构型起点一致，对比干净）；
3. 扩展构型（F51/F42a–c 等）同样从已优化结构分支。

> 反面教材：所有构型各自从原始原型冷启动，会各自落入不同亚稳盆地，能量差被几何噪声污染（3D 旧池 6 个确认结构里 4 个曾因此判错）。

### 7.3 能量口径

- 取 OUTCAR 最后一个 `energy  without entropy`（E0，σ→0）；
- 归一化：ΔE = (E_comp − E_i)/N_M，单位 **meV/磁性原子**；
- **严禁 F 对 E0**：两者差 −σS 项（0.1–0.18 eV），曾导致批量假杀；
- 近简并：|ΔE| < 20 meV/M 用更紧的静态单点（NSW=0、EDIFF 收紧）复核；< 10 meV/M 记"未定"，不作结论；
- 金属性体系注意**电子多稳态**：同几何同磁矩可有多个自洽解（母相差 1.18 eV），跨协议（弛豫 vs 静态）只能定性复核。

### 7.4 收敛口径与资源收割

- 收敛判据：OSZICAR **最后离子步 |ΔE| ≤ 1e-3** 即视为收敛（用户口径；力收敛常卡住，不强求 `reached required accuracy`）；
- 弛豫能量单调下降 → 竞争构型连续两个离子步 E0 低于 i（阈值 0.5 meV）即可 scancel（`tmp/kill_loser.py`）；
- 单作业 > 1 h 仍缓慢优化 → 取消，记 FAIL，不阻塞批次；
- 状态写 `STATUS.json` / `GATE.json` / `VERDICT.json`，支持断点续跑。

### 7.5 磁矩核验（必做）

判决前读 OUTCAR 最后一个 `magnetization (x)` 块：

- 磁矩是否保持初设量级？FiM 初值可能弛豫成 AM，AM 也可能塌成 NM；
- 大 margin 判决尤其要查：Ba₃Mn₂(OsO₄)₃ 的 FM 初值滑到 0.24 μB 低自旋解，导致"i 赢"建立在坍缩对手上；
- 含 f 元素 POTCAR 的磁矩块多一组 f 列，正则要兼容。

**判决可信度分级**：

| 级别 | 条件 |
|---|---|
| 铁案 | 全构型算满 + 双方高自旋 + margin > 20 meV/M |
| 成立但需复核 | 全构型算满 + margin 10–20 meV/M |
| 未定 | margin < 10 meV/M 或对手磁矩坍缩 |
| 暂定 | 关键对手未钉扎（如高自旋 FM 未证） |

---

## 8. 阶段 F：fsg 终判

1. 生成 mcif：`vasp2mcif.py`（CONTCAR + Outcar.magnetization 的 `tot`）或 `outcar2mcif.py`（正则解析末块，兼容 f 列）；
2. `fsg --details <mcif>`，解析 `Magnetic order` / `Leading spin texture` / `MSG with SOC` / `Net moment`；
3. 门控：`texture == i-wave` 且 `|Σm| < 0.5 μB`；
4. 基态构型再判一次（弛豫后纹理可能变：tri-protype 支线 48 个收敛的 AM-g 中，36 个保持 g 波且磁矩 >0.3 μB、4 个变 forbidden、8 个 NM/伪磁塌缩——这是正常损耗，不是 bug）。

---

## 9. 阶段 G：稳定性与物理表征

### 9.1 动力学（声子）

- phonopy 有限位移法；体相超胞 **2×2×2**（z 不是真空，必须扩胞——早期误用 3×3×1 是 2D 习惯）；
- KPOINTS 全 1×1×1（Γ 点）；MAGMOM 用 `10*3 10*-2 24*0.0` 简写；
- 判据：无虚频；有虚频时分解位移方向，判断真实不稳定 vs 磁态假象。

### 9.2 热力学（凸包）

- 参考相：Materials Project API 覆盖元素/二元/三元竞争相；
- 含磁元素参考相必须 ISPIN=2（bcc Fe 等）；
- ΔE = E_cand − Σ xᵢEᵢ；ΔE < 0 稳定，< 50 meV/atom 可称亚稳；
- 大胞高能参考相（K₃As₇、K₃As₁₁）可舍弃但需说明。

### 9.3 MAE / SOC

- `vasp_ncl`，ENCUT=640，ICHARG=11（读 scf CHGCAR），**LMAXMIX=4**（d/f 必须）；
- MAGMOM 必须三分量 `x y z`；
- 固定 |M| 扫方向，拟合 E = E₀ + K sin²φ；
- **六角 y 轴陷阱**：常规晶轴取向下 y 未必是独立高对称方向，判易轴前核对取向；配 rotCHGCAR/fresh 双对照。

### 9.4 能带与劈裂

- 沿磁对称感知路径（AlterSeeK-Path）算无 SOC 能带；
- 劈裂取**同指数（孪生带）配对**，EF ± 0.5 eV 窗口；禁止最近距离配对（交叉处会伪造大劈裂）；
- 教科书判据：对称线（Γ–M）劈裂 ≈ 0（≤1 meV），一般 k 点有限（数十–数百 meV）。

### 9.5 花瓣/费米面

- 在费米面环上采样 ΔE(φ)，角向 FFT 得 c_l，i-purity = |c6|/Σ|c_l|；
- 要求 l=6 主导，极坐标 12 花瓣；
- 自旋分辨费米面（kz=0）应呈六重花瓣。

---

## 10. 数量漏斗（本战役真实数据）

| 阶段 | 数量 | 说明 |
|---|---|---|
| AFLOW 六角原型 168–194 | 248 | 239 直接通过 spglib，9 个有序化后通过 |
| 含 3d + 偶数重轨道 | 63 | 主枚举池 |
| 单轨道染色组合 | 1223 | i-wave **18** |
| 全轨道单染色（放宽） | 4248 | i-wave **77** |
| 装饰候选（有联合构型） | 251 | 8649 个联合构型 |
| 联合枚举 i 波候选 | **238** | SG190 最多（126） |
| dft3 实际弛豫 | 238 | i 收敛 52，门控 PASS 12，取消 186 |
| 纯 3d i 基态（dft3） | **0** | 12 个门过者能量全输（gd/FM） |
| 替换战役 subst（9 个） | 9 | 初判 7 个 i 基态 → 全构型复核后 **6 个** |
| wave12 组合/翻身 | 23 | 全构型复核后 **0 个**（全被 d2/d3/FM/F42 推翻） |
| 桶 A 门控残部 | 10 | LiMSb 7/7、UM₂F₆ 5/5 全关 |
| 桶 C 伪影复活 | 18 | 新增 i×3、g×2 |
| wave15 A 位扫描 | 8 | 新增 i×1（Sr₃Fe₂）、暂定 i×1（Ba₃Mn₂(OsO₄)₃） |
| **最终** | **14 铁案 + 1 暂定** | 另 g 波基态 2 个 |

**家族构成**：

- SG180 P6₂22 f6 `+++---` ×8（MgMnNi 母系：母相、Al₆Mn₆Ni₆、Ca₆Mn₆Ni₆、Sc₆Mn₆Ni₆、Mg₆Mn₆Al₆、Mg₆Mn₆Ga₆、Cd₆Mn₆Ni₆、V₆Mn₆Ni₆）；
- SG189 d2 `+-` ×3（Ca₃Mn₂(IrO₄)₃、Sr₃Mn₂(IrO₄)₃、Sr₃Fe₂(IrO₄)₃）+ 暂定 Ba₃Mn₂(OsO₄)₃；
- SG185 b4 孤立磁四聚体 ×3（Na₇Fe₂As₃、Mn-Na₇Mn₂As₃、Cr₂Cu₇P₃）；
- g 波基态 ×2（Mn₂Cu₇P₃、Fe₂Cu₇P₃）。

---

## 11. 机制认知（决定"下一步往哪打"）

1. **i 波的头号杀手是 d3/d2**：10 个不等价构型里 5 个是 3:3 补偿，i 只是其一；判定必须算满。
2. **化学是波形唯一裁判**：同一几何骨架，换 A 位/框架位就能把 i 推赢或推输（MgMnNi 家族 9 个替换 → 7 i + 2 d）。
3. **位点化学 > 元素本身**：Cd 换 Mg 位 → i 基态；Cd 换 Cu 位 → d3 基态。
4. **A 位杠杆方向依赖**：Ca→Sr/Ba 在 Mn 系把 i 推到 +4 eV，在 Fe/Co/Ni 系却让磁性全灭（NM 坍缩）。
5. **孤立磁四聚体（M–M > 3.8 Å）是另一套物理**：交换全靠超交换，FM / i(++−−) / g(+-+-) 三态同场竞争；传统 3.2 Å 二分图给出的"AFM-G"会退化成 FM（假竞争态），竞争集应改为 {FM, i, g2}。
6. **电子多稳态**：金属性 i 波体系同几何可有多自洽解，跨协议能量不可混比。
7. **4d/5d 的正确用法**：不做磁芯，做框架（IrO₄/OsO₄）；4d¹⁰（Cd²⁺）也可做 A 位。

---

## 12. 编排、目录与状态文件

### 12.1 推荐编排（i 先行 + 门控 + 自动分支）

`orchestrate_v3.py`：

1. 阶段 1：提交全部 i 构型；
2. 阶段 2：i 收敛 → fsg 复检 → PASS（i-wave 且 |净磁矩|<0.5）则**从 i/CONTCAR 分支**提交 gd + FM（以及全构型）；
3. 失败自动重交 1 次，再失败记 FAIL 跳过，不阻塞其他候选；
4. 汇总 `dft3_summary.csv`。

配套：`dft3_watchdog.py`（错误嗅探 + INCAR 修复 + 重置重试）、`tmp/kill_loser.py`（竞争构型两次低于 i 即杀）、`tmp/*_watch.py`（各战役专用）。

### 12.2 目录

```
hex-protype/
├── hexagonal_prototypes/poscar/sg<SG>/    # 248 原型
├── color_scan*.py                         # 染色枚举
├── decorate_new.py                        # 化学装饰
├── prep_dft3_i.py                         # i 构型 DFT 输入
├── orchestrate_v3.py / dft3_watchdog.py   # 编排/看门狗
├── candidates_new/                        # 装饰候选
├── dft/, dft2/, dft3/                     # 主力 DFT 批次
├── dft_revive/, wave12/.../wave15/        # 复活/替换/组合/A 位战役
├── subst_iwave/                           # MgMnNi 家族替换 + 能带/花瓣/费米面
├── out/                                   # CSV 与汇总
└── tmp/                                   # 一次性脚本与门控
```

单候选状态：

```
<structure>/
├── POSCAR/INCAR/KPOINTS/POTCAR/job_vasp.sh
├── i/  gd/  FM/  d2/  d3/  G/  F51/  F42a-c/
├── GATE.json      # {texture, net, keep}
├── VERDICT.json   # {verdict, Ei, Ecomp, margin_meV}
└── STATUS.json    # 流水线状态
```

---

## 13. 硬规则与陷阱清单（按代价排序）

1. **能量口径**：一律 E0 对 E0；禁止 F 对 E0；同批必须同二进制/ENCUT/ISYM/LREAL。
2. **构型要算满**：i/gd/FM 三件套不可靠；C6 位点必须覆盖 10 个代表。
3. **竞争态从 i-CONTCAR 出发**：冷启动会污染能量差。
4. **禁止单轨道冻结**：多轨道磁元素必须联合枚举。
5. **fsg 解析**：stdout+stderr 合并；整句捕获；不吃隐藏文件名；LORBIT=11 必开。
6. **收敛判据**：只认 OSZICAR 末离子步 |ΔE|≤1e-3（或 `reached required accuracy`）；`General timing` 不算。
7. **提交**：`sbatch --chdir` 绝对路径。
8. **OUTCAR 能量正则**：存在 `energy  without entropy=`（等号前无空格）变体，正则要 `\s*=\s*` 且吃 E 记号。
9. **磁矩核验**：大 margin 判决前必查 OUTCAR 末块；对手坍缩的低自旋解不算真对手。
10. **孤立磁四聚体**：3.2 Å 二分图的 AFM-G 是假构型，竞争集用 {FM, i, g2}。
11. **0 字节中间文件**：崩溃残留，重跑前清理。
12. **>1h 慢弛豫**：取消并记录，不占用机时。
13. **特殊分数坐标**：网格扫描落在 1/4、1/3 等特殊点会引入偶然对称操作，把 i 波改判 forbidden；实例选择要跳过特殊分数。

---

## 14. 复现命令清单

```bash
cd /share/home/scuthecc/surface_AM/hex-protype
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1

# 1) 原型 → POSCAR（248 个，spglib 校验）
#    历史脚本已并入 out/ 与 hexagonal_prototypes/，如需重抓：AFLOW aflow.cif → pymatgen

# 2) 单轨道染色
python color_scan.py            # -> out/color_scan.csv        (1223, i=18)
python color_scan_all.py        # -> out/color_scan_all.csv    (4248, i=77)

# 3) 化学装饰（多轨道池）
python decorate_new.py --only-metal-sites   # -> candidates_new/

# 4) 多轨道联合染色 + fsg
python color_scan2_new.py       # -> out/color_scan2_new.csv   (8649, i=238)

# 5) 生成 i 构型 DFT 输入
python prep_dft3_i.py           # -> dft3/<hit>/<cand>/i/

# 6) 编排（i 先行 → 门控 → 竞争构型）
nohup python orchestrate_v3.py > out/orchestrate_v3.log 2>&1 &
nohup python dft3_watchdog.py  > out/watchdog.log      2>&1 &

# 7) mcif + fsg 终判
python ~/bin/vasp2mcif.py <calc_dir>
fsg --details <calc_dir>/<name>.mcif

# 8) 全目录复扫
python tmp/scan_iwave.py        # -> tmp/scan_iwave.csv
```

### 14.1 新增一批候选的标准操作顺序（SOP）

| 步 | 动作 | 产物 | 通过条件 |
|---|---|---|---|
| 1 | 选原型子集（按 SG / 磁性元素 / 轨道重数） | 原型 POSCAR 列表 | 248 池内 |
| 2 | 单轨道染色 + fsg | `out/color_scan*.csv` | i-wave 命中 |
| 3 | 化学装饰（3d 磁芯） | `candidates_new/` | 电荷/CN/去重/入口 fsg 过 |
| 4 | 多轨道联合染色 + fsg | `out/color_scan2_new.csv` | i-wave 且净磁矩 0 |
| 5 | 生成 i 构型输入 | `dft3/*/i/` | INCAR 有 LORBIT=11、ISYM=−1 |
| 6 | 提交 i 弛豫（≤50/批） | 队列 | 指数正常 |
| 7 | 看门狗：收敛 → fsg 门控 | `GATE.json` | texture=i-wave 且 \|Σm\|<0.5 |
| 8 | 从 i/CONTCAR 分支**全不等价**竞争构型 | `gd/ FM/ d2/ d3/ ...` | 按 N 枚举（2/4/10…） |
| 9 | 能量判决 + 磁矩核验 | `VERDICT.json` | i 最低且对手高自旋 |
| 10 | 稳定性/表征 | 声子/凸包/MAE/能带 | 无虚频、ΔE_hull、K、劈裂 |
| 11 | 归档 | SUMMARY/花名册 | — |

**成本量级**（64 核、ENCUT 520）：单候选 i 弛豫 5–40 min（8–30 原子）；全构型竞争（6 位点 10 个）约 1–6 h；一个"i 波家族"（母相 + 6–10 个替换）≈ 1–2 天队列。

### 14.2 高产出配方（按已验证的性价比排序）

| 配方 | 骨架/轨道 | 产出 | 备注 |
|---|---|---|---|
| **A. SG180 f6 `+++---`** | Mg₂Ni 型（A2B_hP18_180_fj_ac） | 8 个 i 基态（最多） | Mn 磁芯 + Ni/Al/Ga/Cu 框架位替换；最成熟 |
| **B. SG189 d2 `+-`** | Ca₃Mn₂(IrO₄)₃ 型 | 3 个 i + 1 暂定 | 3d 磁芯 + 5d 非磁框架；A 位 Ca→Sr 是强杠杆 |
| **C. SG185 b4 孤立四聚体** | Na₇Cu₂As₃ / Cu₃P 衍生 | 3 个 i + 2 个 g | 先标度后数 CN；电荷放行 Zintl 区间；竞争集 {FM, i, g2} |
| **D. SG178–181 12c/6f** | 手性 P6ₓ22 | 单轨道即可出 i（干净） | 12c `++++++------`、6f `+++---`；化学窗口待扫 |
| **E. SG190 b2/f4 多轨道** | 126 个 i 候选（最大池） | 尚未系统判决 | 多轨道联合，注意不要单轨道冻结 |

---

## 15. 局限与下一步

**局限**

- 原型库是代表子集，非 Wyckoff 组合穷举；未覆盖的占位可能藏 i 波。
- 单轨道/多轨道枚举都受倍率上限（12）与组合上限（300）约束。
- 化学过滤只能保证"价态合理"，不能保证可合成；最终仍需凸包/声子。
- i 波在六角框架里确实稀有：238 个群论候选只活下 0 个（纯 3d），说明"几何可行"离"能量可行"很远。

**下一步优先级**

1. **承重墙**：新 i 基态的声子 + 凸包 + MAE（Sr₃Fe₂(IrO₄)₃、Na₇Fe₂As₃、Mn-Na₇Mn₂As₃、Cr₂Cu₇P₃）。
2. **Ba₃Mn₂(OsO₄)₃ 复核**：钉扎高自旋 FM（加大 MAGMOM、从 i CHGCAR 翻转启动），否则只能算暂定。
3. **g 波支线**：tri-protype 的三角 g 波池（611 候选、65 已算、37 个弛豫后仍 g 波）需要补 FM/AFM 竞争才能判决。
4. **borate Ca 替换线**：Sr₂Co₂B₂O₇_Ca 只差 4.4 meV，是唯一未打完的二轮方向。
5. **SOC 纹理家族**：本族 with-SOC 多为 p 波，SOC 主导的 i/p 波交错磁完全未探索。
6. **2D i12 层群池**：24,233 装饰剪枝池中的 i-wave 子集未铺 DFT（但 2D i 波几何上更稀缺）。
