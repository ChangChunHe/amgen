# amgen — Symmetry-First Generator for Altermagnets (Pure Enumeration)

<div align="center">

[![English](https://img.shields.io/badge/lang-English-blue.svg)](README.md)
[![中文](https://img.shields.io/badge/lang-中文-red.svg)](README.zh-CN.md)

</div>

Generates compensated collinear altermagnet candidates **from space-group
decomposition alone** — no external structure database required.

## Core idea

Magnetic wave patterns (d/g/i-wave) are not discovered after the fact — they are
**selected upfront by spin-group constraints on "occupancy + spin coloring",
then turned into real compounds by chemical decoration, and finally validated
as magnetic ground states by DFT**. Two orthogonal criteria must both pass:

| Dimension | Criterion | Tool |
|---|---|---|
| Symmetry | Is this occupancy + coloring an `AFM(Altermagnet)` with the target texture? | fsg |
| Energetics | Is this configuration the magnetic ground state (lowest E0, meV/magnetic atom)? | VASP multi-config competition |

Passing fsg ≠ ground state; lowest energy ≠ i/g-wave.

## The seven-stage pipeline

**① Space-group list (entry point)**
Full hexagonal crystal family: trigonal hP settings SG 143–167 plus hexagonal
SG 168–194 (40 groups). For each SG, pyxtal provides the Wyckoff table
`(letter, multiplicity, # free parameters)` — e.g. SG180 → a/b/c/d (mult 3,
0 params), e–j (mult 6, 1 param), k (mult 12, 3 params).

**② Wyckoff orbit-combination enumeration** (`wyckoff.orbit_combos`)
- Total atoms ≤ 30;
- At least one **even-multiplicity orbit with multiplicity ≥ 4** (the
  candidate magnetic orbit; 2-fold orbits yield no nontrivial compensated
  configurations);
- **The same letter may be reused by different elements** (in Cs₃Fe₂F₉,
  Cs@f and Fe@f are two independent orbit instances) — recipes count
  *orbit instances*, not letters;
- Output sorted by atom count, e.g. SG182 yields `b²+f⁴+g⁶+c²`
  (14 atoms, the BaFe₂O₄ type).

**③ Structure generation** (`wyckoff.generate_ok`)
pyxtal places atoms by recipe; free parameters start random, seeds are retried
until minimum interatomic distance ≥ 1.5 Å; the lattice volume is pre-scaled
from covalent radii (`prescale_volume`) — random pyxtal lattices can be ~50%
off in volume, beyond what a coarse relaxation can repair.

**④ Pre-coloring + fsg gate (zero-cost group-theory filter)**
(`core` + `classify`)
- Enumerate all inequivalent compensated colorings on the candidate magnetic
  orbit (deduplicated under group permutations × global spin flip; only 4
  colorings for a 6-fold orbit, 3 for a 4-fold);
- **Joint coloring** across multiple magnetic orbits (KNiCl₃ lesson: never
  freeze the remaining orbits);
- Each coloring → mcif → fsg; keep only `AFM(Altermagnet)` with g/i-wave
  texture. Frameworks failing the gate are discarded at zero relaxation cost.

**⑤ Chemical decoration** (`decorate`)
Magnetic orbits filled with 3d metals (V–Ni; 4d/5d optional), remaining orbits
filtered by charge neutrality + d-electron count 2–8 + CN ≥ 4; lattice scaled
by covalent-radius ratio. Intermetallics skip the valence filter.

**⑥ CHGNet pre-relaxation + re-symmetrization** (`relax`, 50-step coarse relax)
CHGNet brings free parameters and lattice constants into a reasonable basin
(it is magnetism-blind — geometry only), then a **symprec tolerance ladder**
(1e-3→1e-1) pulls numerical noise back into the parent group. Cases that cannot
return to the parent SG are logged with their final SG (generated SG ≠ relaxed
SG is normal — e.g. CrCl₃: 158→185 with the wave pattern intact).

**⑦ Re-judgment fsg gate + DFT export** (`classify` + `export`)
Coloring and fsg are re-run on the relaxed structure — **guards against
symmetry-upgrade killing the wave** (Cs₃Fe₂F₉ lesson: P-62c→P6₃/mmc closed the
i-wave channel). Only structures passing both gates get VASP inputs
(POSCAR / INCAR with the hit-coloring MAGMOM / KPOINTS / POTCAR / job.sh)
and enter the DFT funnel for full magnetic-configuration competition.

```
SG(40) → orbit combos(~10²–10³) → pre-coloring fsg gate → decoration filters
  → CHGNet 50 steps → re-symmetrize → re-judgment fsg gate → DFT(~10²)
```

## Automatic M / nonmagnetic role enumeration

The recipes in `recipes.py` explicitly declare which orbit is magnetic. When
that is unknown, `roles` removes the assumption:

- `roles.find_orbit_instances()` keeps repeated Wyckoff letters as separate
  orbit instances (unlike letter-only grouping);
- `roles.role_assignments()` enumerates even-multiplicity orbit subsets as
  candidate magnetic `M` roles;
- `roles.joint_role_colorings()` enumerates balanced `+/-` colorings on every
  selected orbit and removes duplicates under the structure's symmetry;
- `roles.scan_roles()` calls fsg and returns only the requested wave type.

On a known framework:

```bash
python examples/enumerate_roles.py POSCAR --target-wave i-wave
```

Or generate an abstract framework first (repeated letters means separate
orbit instances):

```bash
python examples/enumerate_roles.py --sg 180 --orbits f f j --target-wave i-wave
```

On a real decorated structure, restrict candidate magnetic orbits to the
chemically plausible elements and relax the symmetry tolerance when needed:

```bash
python examples/enumerate_roles.py POSCAR --target-wave i-wave \
  --candidate-elements Fe Cr Mn --symprec 0.01
```

The scan answers "for which orbit assignment does this geometry admit the
target spin texture?", not "which atoms are actually magnetic". Without
`--candidate-elements` every even orbit up to multiplicity 12 is tested, so
O/B/Mg orbits can produce symmetry-valid but chemically implausible hits. Element decoration,
post-relaxation fsg re-judgment, and DFT magnetic-ground-state competition are
still required. `max_role_assignments` and `max_colorings` are hard guards;
`max_hits` stops after the requested number (truncation). `scan_roles` by
default keeps only `is_altermagnet=True` hits and can symmetrize the framework
back to a target SG (`resymmetrize_sg`) before running fsg.

## End-to-end candidate enumeration

`enumerate_candidates` wires the whole loop together:
SG → Wyckoff orbit combinations → abstract framework → automatic role scan →
chemical decoration → filters.

```bash
python examples/enumerate_candidates.py --sg 180 --max-atoms 18 \
  --magnetic-elements Mn Fe Cr --nonmagnetic-elements Mg Al Si P S \
  --max-frameworks 200 --max-candidates 20
```

```python
from amgen import enumerate as am_enum
for cand in am_enum.enumerate_candidates(
        sgs=[180, 182], target_wave="i-wave",
        magnetic_elements=["Mn", "Fe", "Cr"],
        nonmagnetic_elements=["Mg", "Al", "Si", "P", "S"],
        min_cn=4, charge_neutral=True):
    print(cand.sg, cand.composition, cand.magnetic_element, cand.pattern)
```

`nonmagnetic_elements=None` keeps the placeholder species (symmetry-only
candidates). Element decoration is orbit-instance aware: repeated Wyckoff
letters are replaced independently.

## Dependencies

- pymatgen, spglib, pyxtal, chgnet, numpy
- `fsg` (FindSpinGroup) on PATH or at `~/miniconda3/bin/fsg`
- `vaspkit` (POTCAR generation, only for `--write-vasp`)

## Regression test: four confirmed families

Recipes live in `amgen/recipes.py` — pure group-theory declarations, no
database structure files:

| Family | SG | Recipe | Magnetic orbit | Target coloring | Wave |
|---|---|---|---|---|---|
| Mg₆Mn₆Ni₆ | 180 P6₂22 | Mn@f, Mg@j, Ni@a+c | f⁶ | `+++---` | i-wave |
| BaFe₂O₄ | 182 P6₃22 | Ba@b, Fe@f, O@g+c | f⁴ | `++--` | i-wave |
| Cs₃Cr₂Br₉ | 194 P6₃/mmc | Cs@d+f, Cr@e, Br@h+k | e⁴ | `+-+-` | g-wave |
| Mg₂Cr₂B₂O₇ | 188 P-6c2 | Mg@e+f, Cr@g, B@h, O@a+l | g⁴ | `++--` | i-wave |

Run:

```bash
python examples/reproduce_families.py                    # serial, login node
python examples/reproduce_families.py --only BaFe2O4     # single family
python examples/reproduce_families.py --write-vasp out/  # also write DFT inputs
sbatch job_reproduce.sh   # cluster: one 64-core node, 4 families × 16 threads
```

Each family runs the full seven-stage pipeline and is compared against the
confirmed DFT structure (StructureMatcher; falls back to prototype-framework
matching — same SG + orbit-multiplicity spectrum + similar lattice — since a
decorated structure legitimately differs in lattice from a relaxed reference).

Current status: **all four families PASS** (three exact StructureMatcher
matches; Cs₃Cr₂Br₉ framework-matched with ~20% lattice offset that DFT
relaxation absorbs).

## API quick tour

```python
from amgen import wyckoff, relax, core, classify, decorate, export, roles

wyckoff.wyckoff_table(182)                     # SG -> Wyckoff table
wyckoff.orbit_combos(182, max_atoms=30)        # orbit-combination enumeration
s = wyckoff.generate_ok(182, {"Ba": ["b"], "Fe": ["f"], "O": ["g", "c"]})

# unknown magnetic orbit: enumerate M/nonmagnetic roles, then colorings
for hit in roles.scan_roles(s, "i-wave", candidate_elements=["Fe"]):
    print(hit.magnetic_orbits, hit.pattern, hit.wave)

s2, sg2 = relax.relax_and_symmetrize(s, target_sg=182)   # CHGNet + back to parent SG

orbits = core.find_orbits(s2)
pats   = core.orbit_colorings(coords, s2.lattice)
r = classify.fsg_classify(s2, mag_idx, signs)  # -> texture: i-wave / g-wave / forbidden

s3 = decorate.substitute_orbits(s2, {"f": "Mn"})
export.write_vasp("out/X/i", s3, magmom_list)
```

## Known pitfalls (hard-won)

- **Letter reuse**: enumerate by orbit instance, not by letter;
- **CHGNet is magnetism-blind**: Jahn–Teller and spin-driven distortions are
  invisible to it — geometry pre-relaxation only;
- **Relaxation can upgrade symmetry and kill the wave** (Cs₃Fe₂F₉) — the
  re-judgment gate is mandatory;
- **mcif/fsg tolerance**: relaxation noise misjudges as P1 at symprec 1e-3 —
  use the tolerance ladder for symmetrization;
- VASP side: lattice noise triggers the Bravais-lattice error (idealize
  hexagonal cells); EDDDAV → retry with `ALGO=A`; symmetrization errors →
  `ISYM=-1`; POTCAR element mix-ups are not caught by VASP — grep TITEL.

## Repository layout

```
amgen/
  amgen/{wyckoff,relax,core,classify,decorate,export,recipes,roles,enumerate}.py
  examples/reproduce_families.py     # four-family regression test
  examples/enumerate_roles.py        # automatic M/nonmagnetic role scan
  examples/enumerate_candidates.py   # end-to-end candidate enumeration
  job_reproduce.sh                   # cluster submission (parallel)
```

Full campaign logs and negative-result analyses (in Chinese):
[`docs/HEX_IWAVE_GENERATION_SCHEME.md`](docs/HEX_IWAVE_GENERATION_SCHEME.md),
[`docs/IWAVE_FINAL_SUMMARY.md`](docs/IWAVE_FINAL_SUMMARY.md),
[`docs/CAMPAIGN_LOG.md`](docs/CAMPAIGN_LOG.md).
