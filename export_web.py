"""
Convert data/materials.csv into the compact, column-oriented JSON the web page loads.

    python scripts/export_web.py            # reads data/materials.csv, writes data/materials.json
"""
import json
import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "data" / "materials.csv"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "data" / "materials.json"

# Colour groups: at most 8 hues + a neutral "Other", so families stay distinguishable.
# The detailed family label is still shown on hover.
GROUP_OF = {
    "TMD (MX₂)": "TMD (MX₂)",
    "Other chalcogenide": "Other chalcogenide",
    "Oxide": "Oxide / hydroxide",
    "Hydroxide / oxyhydride": "Oxide / hydroxide",
    "Halide": "Halide",
    "Oxyhalide": "Mixed-anion halide",
    "Chalcohalide": "Mixed-anion halide",
    "Pnictide": "Pnictide",
    "Carbide / boride / silicide / germanide": "Carbide / boride / silicide",
    "Intermetallic": "Intermetallic",
    "Elemental": "Elemental & other",
    "MXene": "Elemental & other",
    "Hydride": "Elemental & other",
}

NUMERIC = [
    "bandgap_eV", "vbm_eV", "cbm_eV",
    "exfoliation_energy_meV_A2", "exfoliation_energy_per_atom_eV",
    "decomposition_energy_per_atom_eV", "energy_per_atom_eV",
    "magnetization_per_atom_muB", "thickness_A", "area_per_atom_A2",
    "atoms_per_nm2", "areal_mass_density_ug_cm2", "mean_atomic_mass",
    "a_A", "b_A", "nelements", "n_atoms_cell",
]
TEXT = ["material_id", "formula", "family", "discovery_process", "spacegroup", "relative_id",
        "crystal_system", "source_id"]


def sig(x, n=5):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return None
    if x == 0:
        return 0
    return float(f"{x:.{n}g}")


df = pd.read_csv(SRC)
df["group"] = df["family"].map(GROUP_OF).fillna("Elemental & other")
df["is_metal"] = df["is_metal"].fillna(False).astype(bool)

out = {"n": len(df), "columns": {}}
for c in TEXT + ["group"]:
    out["columns"][c] = df[c].fillna("").astype(str).tolist()
for c in NUMERIC:
    out["columns"][c] = [sig(v) for v in df[c].tolist()]
out["columns"]["is_metal"] = [int(v) for v in df["is_metal"]]

OUT.write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")))
print(f"wrote {len(df)} materials -> {OUT} ({OUT.stat().st_size/1e6:.2f} MB)")
print(df["group"].value_counts().to_string())
