# Tanglegram Visualization Toolkit

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Compare two phylogenetic trees built from the same genomes — for example an alignment-free
**VAMPhyRE** tree and a **GTDB-Tk** reference tree — and draw publication-ready figures of
either tree on its own.

The toolkit quantifies where the trees disagree (Robinson–Foulds distance, "rogue" genomes,
discordant clades, non-monophyletic taxa), colours everything by GTDB taxonomy, collapses
clades to keep large trees readable, and exports editable files for **iTOL**.

It runs as a command-line tool, as a local web app (Streamlit), or as a Python library.
All three share the same code path and produce identical files.

---

## Features

**Tanglegram (two trees)**
- Robinson–Foulds (RF) and normalised RF distance; split-level comparison independent of rooting.
- Greedy detection of rogue genomes (the leaves whose removal lowers RF the most).
- Per-node concordance, per-taxon monophyly and fragmentation tables.
- Common rooting and crossing minimisation (ladderize + untangle) for a readable layout.

**Single tree**
- Rectangular or **circular** layouts (full circle or partial arc).
- Branch-support display (number or dot, with a threshold) and support-based polytomy
  collapse (like `ape::di2multi`). GTDB-Tk 0–1 support is normalised to 0–100.
- Rooting: as in the file, midpoint, or on an outgroup taxon.
- **iTOL export**: tree + annotation datasets (colours, taxonomy strip with legend,
  readable labels, collapsed clades) bundled as a ZIP for iTOL *Batch upload*.

**Shared options**
- Colour by any GTDB rank (domain → species), highlight taxa, zoom into a subset,
  collapse pure clades, label leaves.
- Colour-blind-safe palette by default (8 colours); more colours are available on request,
  with explicit accessibility warnings (see [Colours and accessibility](#colours-and-accessibility)).
- PNG output plus optional SVG/PDF; CSV tables and a Markdown report for every run.
- Every run prints the **equivalent command line**, so any figure can be reproduced.

---

## Installation

Requires Python ≥ 3.12 (tested). For older Python versions, see `requirements-flexible.txt`.

```bash
git clone https://github.com/Emilianodz/tanglegram_visualization_toolkit.git
cd tanglegram_visualization_toolkit
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate

pip install -r requirements.txt        # command line + library
pip install -r requirements-app.txt    # + web interface (includes everything above)
pip install -r requirements-dev.txt    # + test suite
```

| File | Purpose |
|---|---|
| `requirements.txt` | Exact versions of the tested core dependencies |
| `requirements-app.txt` | Core + Streamlit for the web interface |
| `requirements-flexible.txt` | Minimum versions for older Python (untested) |
| `requirements-dev.txt` | Core + pytest and DendroPy (independent RF check) |

---

## Quick start

### Web interface

```bash
streamlit run app.py               # opens http://localhost:8501
```

Two pages: **Tanglegram** and **Single tree**. Upload the files, adjust the sidebar,
press *Generate*, and download the figure, tables and report (individually or as a ZIP).
Files are processed in a local temporary folder; nothing is sent to an external server.

### Command line

```bash
# Tanglegram: VAMPhyRE vs GTDB-Tk, coloured by family
python generate_tanglegram.py \
  --tree-vamphyre vamphyre.newick \
  --tree-gtdb gtdb_relabelled.newick \
  --tree-gtdb-original gtdbtk.unrooted.tree \
  --tabla gtdbtk.ar53.summary.tsv \
  --out-plot results/tanglegram.png

# Single circular tree, coloured by genus, with iTOL export
python generate_tree.py \
  --arbol gtdb_relabelled.newick \
  --tabla gtdbtk.ar53.summary.tsv \
  --mapa-ids results/tanglegram_mapa_ids.tsv \
  --color-rango genero --soporte punto --soporte-min 90 \
  --out-plot results/gtdb_circular.png
```

Run `python generate_tanglegram.py --help` or `python generate_tree.py --help` for every option.
The tanglegram is also available as `python -m tanglegram_utils`.

> **Note on option names.** Help text, messages, reports and the web interface are in English,
> but command-line flag names and output file suffixes keep their original Spanish names
> (e.g. `--tabla` = table, `--color-rango` = colour rank, `--resaltar` = highlight) for
> backward compatibility. The [User Guide](USER_GUIDE.md) includes a full English reference.

See the **[User Guide](USER_GUIDE.md)** for inputs, outputs, recipes and troubleshooting.

---

## Colours and accessibility

`--color-top N` sets how many taxa get their own colour; the rest are grey ("other").

| N | Palette | Recommendation |
|---|---|---|
| 1–8 | Validated categorical palette | **Safe for colour-vision deficiency.** Recommended for publication. |
| 9–20 | 12 extra hand-picked tones | Not validated for colour blindness. Pair with highlighting. |
| > 20 | Generated tones (golden-angle hue spacing) | Hard to tell apart even with normal vision and **not suitable for readers with colour blindness**. Use for exploration only. |

The tool prints a warning (and the web app shows one) whenever N exceeds 8 or 20.
For accessible figures with many groups, prefer: colouring at a higher rank,
`--resaltar` to highlight a few taxa, collapsing pure clades, or zooming with `--solo-taxa`.
The first 20 colours never change, so figures stay consistent across runs.

---

## iTOL export

`generate_tree.py` (and the *Single tree* page) writes, next to the figure:

| File | iTOL dataset | Content |
|---|---|---|
| `<prefix>_itol.tree` | Tree | Newick with branch lengths and support values |
| `<prefix>_itol_colores.txt` | `TREE_COLORS` | Branch colour of each pure clade and leaf label colours |
| `<prefix>_itol_franja.txt` | `DATASET_COLORSTRIP` | Taxonomy colour strip with legend |
| `<prefix>_itol_etiquetas.txt` | `LABELS` | Readable leaf names (GTDB species) instead of genome IDs |
| `<prefix>_itol_colapso.txt` | `COLLAPSE` | Same collapsed clades as the figure |
| `<prefix>_itol.zip` | — | All of the above, for *Batch upload* |

Upload `<prefix>_itol.zip` via **Upload → Batch upload** on [itol.embl.de](https://itol.embl.de),
or upload the `.tree` file and drag the `.txt` files onto the tree. Each layer stays a separate,
editable dataset. Internal nodes are referenced as `leafA|leafB` (last common ancestor), so
annotations do not depend on internal node names. Disable with `--sin-itol`.

---

## Project layout

```
├── generate_tanglegram.py   # CLI entry point: two-tree tanglegram
├── generate_tree.py         # CLI entry point: single tree (rectangular/circular, support, iTOL)
├── app.py                   # Streamlit web app (navigation)
├── paginas/                 # Web app pages: Tanglegram, Single tree
├── ui_comun.py              # Shared web UI components
├── tanglegram_utils/        # All logic, as reusable modules
└── tests/                   # pytest suite (synthetic data only)
```

| Module | Responsibility |
|---|---|
| `tree` | Tree structure, Newick I/O, pruning, bitmask indexing |
| `splits` | Splits, RF distance, rogue detection, per-node discordance |
| `taxonomy` | GTDB taxonomy parsing, taxon selection, colour scheme and palette |
| `idmap` | Table loading and name harmonisation between table and trees |
| `rooting`, `ordering`, `collapse` | Rooting, ladderize/untangle, clade collapsing |
| `support` | Support-value parsing, normalisation and low-support collapse |
| `layout`, `plotting`, `circular` | Geometry and matplotlib drawing |
| `itol` | iTOL tree and annotation export |
| `report`, `report_arbol` | CSV tables and Markdown reports |
| `pipeline`, `pipeline_arbol` | End-to-end orchestration |
| `cli`, `cli_arbol`, `models` | Argument parsing and data classes |
| `servicio` | Run-in-directory, log capture, ZIP packaging (used by the web app) |

---

## Library use

```python
from tanglegram_utils.models import Entradas, OpcionesAnalisis, Salida
from tanglegram_utils.pipeline import ejecutar

datos, an = ejecutar(
    Entradas("vamphyre.newick", "gtdb.newick", tabla="summary.tsv",
             tree_gtdb_original="gtdbtk.unrooted.tree"),
    OpcionesAnalisis(color_rango="genero"),
    salida=Salida(out_plot="results/tanglegram.png"),
)
print(an.rf, an.nrf)
```

---

## Testing

```bash
pip install -r requirements-dev.txt
pytest
```

The suite uses small synthetic trees only, so no research data is needed.

---

## Methodological notes

- RF and splits are rooting-independent; the common root is used only for drawing.
- The number of fragments per taxon depends on the chosen root; RF does not.
- `--contraer-soporte U` turns nodes with support < U into polytomies; leaf-to-leaf distances
  across the removed branch become shorter.
- GTDB-Tk support and bootstrap values are different confidence measures — do not compare
  them directly.
- VAMPhyRE internal labels such as `Inner123` are node names, not support, and are ignored.

---

## License

Released under the [MIT License](LICENSE).
