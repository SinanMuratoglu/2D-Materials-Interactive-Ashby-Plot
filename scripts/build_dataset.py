"""
Build a clean, flat table from the 2DMatPedia dump (newline-delimited JSON).

Output: materials.csv — one row per 2D material, with
  - identifiers and metadata
  - DFT properties straight from 2DMatPedia
  - geometry-derived properties computed from the relaxed structure
  - a rule-based 'family' label for colour-coding / blobs
"""
import json
import sys

import numpy as np
import pandas as pd
import periodictable

SRC = sys.argv[1] if len(sys.argv) > 1 else "2dmatpedia.jsonl"
OUT = sys.argv[2] if len(sys.argv) > 2 else "materials.csv"

# ---------------------------------------------------------------- element data
MASS = {el.symbol: el.mass for el in periodictable.elements if el.number > 0}

# Pauling electronegativities (used to decide which element is the "anion")
EN = {
    "H": 2.20, "Li": 0.98, "Be": 1.57, "B": 2.04, "C": 2.55, "N": 3.04, "O": 3.44, "F": 3.98,
    "Na": 0.93, "Mg": 1.31, "Al": 1.61, "Si": 1.90, "P": 2.19, "S": 2.58, "Cl": 3.16,
    "K": 0.82, "Ca": 1.00, "Sc": 1.36, "Ti": 1.54, "V": 1.63, "Cr": 1.66, "Mn": 1.55,
    "Fe": 1.83, "Co": 1.88, "Ni": 1.91, "Cu": 1.90, "Zn": 1.65, "Ga": 1.81, "Ge": 2.01,
    "As": 2.18, "Se": 2.55, "Br": 2.96, "Rb": 0.82, "Sr": 0.95, "Y": 1.22, "Zr": 1.33,
    "Nb": 1.60, "Mo": 2.16, "Tc": 1.90, "Ru": 2.20, "Rh": 2.28, "Pd": 2.20, "Ag": 1.93,
    "Cd": 1.69, "In": 1.78, "Sn": 1.96, "Sb": 2.05, "Te": 2.10, "I": 2.66, "Cs": 0.79,
    "Ba": 0.89, "La": 1.10, "Ce": 1.12, "Pr": 1.13, "Nd": 1.14, "Pm": 1.13, "Sm": 1.17,
    "Eu": 1.20, "Gd": 1.20, "Tb": 1.10, "Dy": 1.22, "Ho": 1.23, "Er": 1.24, "Tm": 1.25,
    "Yb": 1.10, "Lu": 1.27, "Hf": 1.30, "Ta": 1.50, "W": 2.36, "Re": 1.90, "Os": 2.20,
    "Ir": 2.20, "Pt": 2.28, "Au": 2.54, "Hg": 2.00, "Tl": 1.62, "Pb": 2.33, "Bi": 2.02,
    "Th": 1.30, "Pa": 1.50, "U": 1.38, "Np": 1.36, "Pu": 1.28, "Ac": 1.10,
}

CHALCOGENS = {"S", "Se", "Te"}
HALOGENS = {"F", "Cl", "Br", "I"}
PNICTOGENS = {"N", "P", "As", "Sb", "Bi"}
TETRELS_ETC = {"C", "Si", "Ge", "B"}
NONMETALS = CHALCOGENS | HALOGENS | PNICTOGENS | TETRELS_ETC | {"O", "H", "Se", "Te"}
TRANSITION_METALS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
}
MXENE_METALS = {"Sc", "Ti", "V", "Cr", "Y", "Zr", "Nb", "Mo", "Hf", "Ta", "W"}


def anion_group(el):
    if el == "O":
        return "oxide"
    if el in HALOGENS:
        return "halide"
    if el in CHALCOGENS:
        return "chalcogenide"
    if el in {"N", "P", "As", "Sb"}:
        return "pnictide"
    if el in {"C", "Si", "Ge", "B"}:
        return "carbide/tetrelide"
    if el == "H":
        return "hydride"
    return None


ANION_PRIORITY = [  # first match wins
    ("halide", HALOGENS),
    ("oxide", {"O"}),
    ("chalcogenide", CHALCOGENS),
    ("pnictide", {"N", "P", "As"}),
    ("carbide/tetrelide", {"C", "Si", "B"}),
    ("heavy pnictide", {"Sb", "Bi"}),   # Zintl-type, only if nothing above present
    ("tetrelide", {"Ge"}),
    ("hydride", {"H"}),
]
NONMETAL_LIKE = set().union(*[s for _, s in ANION_PRIORITY])


def classify(comp):
    """comp: dict element -> amount (reduced). Returns a family label."""
    els = set(comp)
    if len(els) == 1:
        return "Elemental"
    if not (els & NONMETAL_LIKE):
        return "Intermetallic"
    # MXene-like: early TM + C/N, metal-rich, optionally O/F/H terminations
    core = els - {"O", "F", "H"}
    if (core & MXENE_METALS) and (core & {"C", "N"}) and core <= (MXENE_METALS | {"C", "N"}):
        m = sum(comp[e] for e in core & MXENE_METALS)
        x = sum(comp[e] for e in core & {"C", "N"})
        if m > x:
            return "MXene"

    # TMDs: binary transition metal + S/Se/Te in MX2 ratio
    if len(els) == 2 and (els & TRANSITION_METALS) and (els & CHALCOGENS):
        (tm,) = els & TRANSITION_METALS
        (ch,) = els & CHALCOGENS
        if abs(comp[ch] / comp[tm] - 2) < 1e-6:
            return "TMD (MX₂)"

    present = [g for g, s in ANION_PRIORITY if els & s]
    top = present[0]
    if top == "halide":
        if "oxide" in present:
            return "Oxyhalide"
        if "chalcogenide" in present:
            return "Chalcohalide"
        return "Halide"
    if top == "oxide":
        return "Hydroxide / oxyhydride" if "hydride" in present else "Oxide"
    if top == "chalcogenide":
        return "Other chalcogenide"
    if top in ("pnictide", "heavy pnictide"):
        return "Pnictide"
    if top in ("carbide/tetrelide", "tetrelide"):
        return "Carbide / boride / silicide / germanide"
    if top == "hydride":
        return "Hydride"
    return "Other"


def reduced_comp(formula_abc):
    # e.g. "F2 Ir1" -> {"F": 2.0, "Ir": 1.0}
    out = {}
    for tok in formula_abc.split():
        i = next(k for k, ch in enumerate(tok) if ch.isdigit() or ch == ".")
        out[tok[:i]] = float(tok[i:])
    return out


def geometry(struct):
    lat = np.array(struct["lattice"]["matrix"])
    a, b, c = lat
    area = float(np.linalg.norm(np.cross(a, b)))              # in-plane area, Å^2
    c_perp = float(abs(np.dot(c, np.cross(a, b))) / area)      # cell height ⟂ to layer, Å
    fz = np.array([s["abc"][2] % 1.0 for s in struct["sites"]])
    # Layer may wrap across the periodic boundary: thickness = cell height minus
    # the largest empty gap in fractional z (that gap is the vacuum).
    zs = np.sort(fz)
    gaps = np.diff(np.concatenate([zs, [zs[0] + 1.0]]))
    thickness = c_perp * (1.0 - gaps.max()) if len(zs) > 1 else 0.0
    mass = sum(MASS[sp["element"]] * sp["occu"] for s in struct["sites"] for sp in s["species"])
    nat = len(struct["sites"])
    return dict(
        n_atoms_cell=nat,
        area_A2=area,
        area_per_atom_A2=area / nat,
        thickness_A=thickness,
        areal_mass_density_ug_cm2=mass * 1.66054e-24 / (area * 1e-16) * 1e6,  # µg/cm²
        atoms_per_nm2=nat / (area / 100.0),
        mean_atomic_mass=mass / nat,
        a_A=struct["lattice"]["a"],
        b_A=struct["lattice"]["b"],
    )


rows = []
with open(SRC) as fh:
    for line in fh:
        r = json.loads(line)
        comp = reduced_comp(r["formula_reduced_abc"])
        g = geometry(r["structure"])
        bs = r.get("bandstructure", {}) or {}
        rows.append(dict(
            material_id=r["material_id"],
            formula=r["formula_pretty"],
            chemsys=r["chemsys"],
            nelements=r["nelements"],
            elements=" ".join(sorted(r["elements"])),
            formula_anonymous=r["formula_anonymous"],
            discovery_process=r["discovery_process"],
            source_id=r.get("source_id"),
            relative_id=r.get("relative_id"),
            spacegroup=r["sg_symbol"],
            spacegroup_number=r["sg_number"],
            crystal_system=r["spacegroup"]["crystal_system"],
            family=classify(comp),
            bandgap_eV=r.get("bandgap"),
            is_metal=bs.get("is_metal"),
            is_gap_direct=bs.get("is_gap_direct"),
            vbm_eV=bs.get("vbm"),
            cbm_eV=bs.get("cbm"),
            energy_per_atom_eV=r.get("energy_per_atom"),
            energy_vdw_per_atom_eV=r.get("energy_vdw_per_atom"),
            exfoliation_energy_per_atom_eV=r.get("exfoliation_energy_per_atom"),
            decomposition_energy_per_atom_eV=r.get("decomposition_energy"),
            total_magnetization_muB=r.get("total_magnetization"),
            magnetization_per_atom_muB=(abs(r["total_magnetization"]) / g["n_atoms_cell"]
                                        if r.get("total_magnetization") is not None else None),
            **g,
        ))

df = pd.DataFrame(rows)
# exfoliation energy per unit area (meV/Å²) is the more standard way to report it
df["exfoliation_energy_meV_A2"] = (df["exfoliation_energy_per_atom_eV"] * 1000
                                   / df["area_per_atom_A2"])
df.to_csv(OUT, index=False)
print(f"wrote {len(df)} rows -> {OUT}")
