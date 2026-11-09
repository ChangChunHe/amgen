"""amgen.classify — mcif 生成 + fsg 波形判定"""
import os, re, shutil, subprocess
from pymatgen.io.cif import CifWriter

FSG = shutil.which("fsg") or os.path.expanduser("~/miniconda3/bin/fsg")

# 角动量阶次 -> 波型标签（fsg texture 字符串已含标签，这里只做兜底映射）
ORDER2WAVE = {1: "p-wave", 2: "d-wave", 3: "f-wave", 4: "g-wave",
              5: "h-wave", 6: "i-wave"}


def fsg_classify(struct, mag_idx, signs, moment=2.0, symprec_note=""):
    """给结构加共线 z 向磁矩后跑 fsg。

    返回 dict: mag_order, texture ('i-wave' 或 'forbidden' 等), order, msg, wave
    signs 为 ±1 列表（与 mag_idx 对应）；moment 可以是标量或 {元素: 磁矩} 字典。
    """
    s = struct.copy()
    moms = [0.0] * len(s)
    for i, sg in zip(mag_idx, signs):
        m = moment.get(str(s[i].specie.symbol), 2.0) if isinstance(moment, dict) else moment
        moms[i] = m * sg
    s.add_site_property("magmom", moms)
    tmp = "/tmp/amgen_%d.mcif" % os.getpid()
    CifWriter(s, write_magmoms=True).write_file(tmp)
    try:
        out = subprocess.run([FSG, tmp], capture_output=True, text=True,
                             timeout=120).stdout
    except Exception:
        return {"mag_order": "ERROR", "texture": "?", "order": 0, "msg": "", "wave": "?"}
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    mo = re.search(r"Magnetic order: (.+)", out)
    tx = re.search(r"Leading spin texture: without SOC (.+?) \(order (\d+)", out)
    mg = re.search(r"MSG with SOC: (\S+) (\S+) \(type (\d)\)", out)
    order = int(tx.group(2)) if tx else 0
    texture = tx.group(1).strip() if tx else "forbidden"
    return {
        "mag_order": mo.group(1).strip() if mo else "?",
        "texture": texture,
        "order": order,
        "msg": ("%s %s t%s" % mg.groups()) if mg else "",
        "wave": texture if texture.endswith("wave") else ORDER2WAVE.get(order, "forbidden"),
        "is_altermagnet": "Altermagnet" in (mo.group(1) if mo else ""),
    }
