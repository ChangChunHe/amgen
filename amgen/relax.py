"""amgen.relax — CHGNet 预弛豫 + 对称化容差阶梯

CHGNet 不懂磁性: 只用于把自由参数/晶格调到合理邻域并淘汰键长荒谬者,
波形/稳定性判决一律交给后面的 DFT。弛豫引入的数值噪声用容差阶梯吸回母群。
"""
import warnings
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

warnings.filterwarnings("ignore")

_RELAXER = None


def _get_relaxer():
    global _RELAXER
    if _RELAXER is None:
        from chgnet.model import StructOptimizer
        _RELAXER = StructOptimizer()
    return _RELAXER


def chgnet_relax(struct, fmax=0.2, steps=50, relax_cell=True):
    """CHGNet 粗弛豫, 返回弛豫后 Structure; 失败返回 None"""
    try:
        relaxer = _get_relaxer()
        res = relaxer.relax(struct, fmax=fmax, steps=steps,
                            relax_cell=relax_cell, verbose=False)
        return res["final_structure"]
    except Exception:
        return None


def resymmetrize(struct, target_sg, ladder=(1e-3, 3e-3, 1e-2, 3e-2, 5e-2, 1e-1)):
    """容差阶梯对称化: 取仍能给出 target_sg 的最严容差, 返回精修结构。
    若所有容差都回不到 target_sg, 返回 (None, 实际SG号)。
    """
    for tol in ladder:
        try:
            sga = SpacegroupAnalyzer(struct, symprec=tol)
            if sga.get_space_group_number() == target_sg:
                return sga.get_refined_structure(), target_sg
        except Exception:
            continue
    try:
        actual = SpacegroupAnalyzer(struct, symprec=0.1).get_space_group_number()
    except Exception:
        actual = -1
    return None, actual


def relax_and_symmetrize(struct, target_sg, fmax=0.2, steps=50):
    """CHGNet 弛豫 + 回母群对称化。返回 (structure or None, 终态SG号)"""
    s = chgnet_relax(struct, fmax=fmax, steps=steps)
    if s is None:
        return None, -1
    return resymmetrize(s, target_sg)
