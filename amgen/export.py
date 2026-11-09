"""amgen.export — DFT 输入文件生成（POSCAR/INCAR/KPOINTS/POTCAR/job.sh）"""
import subprocess
from pathlib import Path
from .decorate import compress_magmom

INCAR_RELAX = """EDIFF = 1e-06
EDIFFG = -0.01
ENCUT = {encut}
IBRION = 2
ICHARG = 2
ISIF = 3
ISMEAR = 0
ISPIN = 2
ISTART = 0
LCHARG = F
LREAL = Auto
LWAVE = F
NCORE = 4
NSW = {nsw}
LORBIT = 11
MAGMOM = {magmom}
"""

JOB_TMPL = """#!/bin/bash -l
#SBATCH -J {name}
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --partition=cpuXeon6458
#SBATCH --exclude=c22n03,c22n04,c01n05,c01n04
#SBATCH --chdir={workdir}

module load vasp/6.4.0
mpirun -n ${{SLURM_NPROCS}} /share/home/scuthecc/bin/vasp_std > vasp.out 2>&1
"""


def gamma_kpoints(a, b, c, density=0.25):
    """按实空间晶格常数给 Γ 中心网格（约 kppa~2000-3000 口径）"""
    import math
    return [max(1, int(round(density * 2 * math.pi / x / 0.628))) for x in (a, b, c)]


def write_vasp(dirpath, struct, magmom, encut=520, nsw=100, kmesh=None,
               potcar=True, job_name=None, submit_script=True):
    """写一个弛豫目录。magmom 为全原子长度列表。返回 Path。"""
    d = Path(dirpath)
    d.mkdir(parents=True, exist_ok=True)
    struct.to(filename=str(d / "POSCAR"), fmt="poscar")
    (d / "INCAR").write_text(INCAR_RELAX.format(
        encut=encut, nsw=nsw, magmom=compress_magmom(list(magmom))))
    if kmesh is None:
        lat = struct.lattice
        kmesh = gamma_kpoints(lat.a, lat.b, lat.c)
    (d / "KPOINTS").write_text(
        "Gamma %dx%dx%d\n0\nGamma\n%d %d %d\n0 0 0\n"
        % (kmesh[0], kmesh[1], kmesh[2], kmesh[0], kmesh[1], kmesh[2]))
    if potcar:
        r = subprocess.run("vaspkit -task 103", shell=True, cwd=d,
                           capture_output=True, text=True, timeout=120)
        if not (d / "POTCAR").exists():
            raise RuntimeError(f"vaspkit POTCAR 失败: {d}\n{r.stdout}{r.stderr}")
    if submit_script:
        name = job_name or d.name
        (d / "job.sh").write_text(JOB_TMPL.format(name=name[:12], workdir=d))
    return d
