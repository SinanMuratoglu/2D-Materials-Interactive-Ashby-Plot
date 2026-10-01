# 2D Materials Interactive Ashby Plot

An interactive, Ashby-style property map of ~6,350 two-dimensional materials from
[2DMatPedia](http://www.2dmatpedia.org/). Pick any two properties for the axes, colour
by material family, and wrap each family's densest region in a smooth envelope.

**Live site:** https://sinanmuratoglu.github.io/2D-Materials-Interactive-Ashby-Plot/

## Features

- 17 selectable axes, each with a linear or log scale: band gap, exfoliation energy, energy above hull,
  areal mass density, layer thickness, magnetization, band edges, lattice constants and more
- Every phase/polymorph is its own point (e.g. 2H and 1T′ MoS₂ are separate)
- Nine colour-coded families. Click a family to toggle it, double-click to isolate it, hover to highlight it
- Family envelopes: a kernel-density contour enclosing a chosen share (50–95%) of each
  family, drawn behind the dots or on its own
- Formula search: matching phases are ringed and labelled with their space group
- Filters: top-down vs bottom-up, metals, unphysical (negative) exfoliation energies, near-stable only
- Click a point to open its 2DMatPedia page; export the visible data as CSV or the plot as PNG
- The view is stored in the URL, so a link reproduces the same plot

## Repository layout

```
index.html                 the whole app (static; Plotly.js from CDN)
data/materials.json        compact data the page loads
data/materials.csv         full cleaned table (one row per material)
scripts/build_dataset.py   2DMatPedia JSON dump  -> data/materials.csv
scripts/export_web.py      data/materials.csv    -> data/materials.json
```

## Rebuilding the data

1. Download the 2DMatPedia database dump from http://www.2dmatpedia.org/download.
   It is newline-delimited JSON, about 29 MB uncompressed.
2. Run:
   ```bash
   pip install numpy pandas periodictable
   python scripts/build_dataset.py db.json data/materials.csv
   python scripts/export_web.py
   ```

## Notes on the data

- Values are DFT (PBE) from 2DMatPedia. PBE underestimates band gaps.
- **Derived properties** come from the relaxed structures: layer thickness is the atomic span normal to
  the layer (0 Å for graphene), and areal densities use the in-plane cell area.
- **Exfoliation energy per area** = exfoliation energy per atom ÷ area per atom.
- **Families** are assigned by composition rules in `scripts/build_dataset.py`. "TMD (MX₂)"
  means binary transition-metal + S/Se/Te in a 1:2 ratio. Other chalcogenides, such as Ta₂Te₃,
  fall under "Other chalcogenide". The finer family label shows on hover.
- 2DMatPedia contains no mechanical properties (modulus, strength).

## Citation

J. Zhou et al., "2DMatPedia, an open computational database of two-dimensional materials
from top-down and bottom-up approaches," *Sci. Data* **6**, 86 (2019).
https://doi.org/10.1038/s41597-019-0097-3
