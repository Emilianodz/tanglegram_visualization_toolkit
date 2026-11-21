# User Guide

This guide covers inputs, outputs, common workflows and the full option reference for the
Tanglegram Toolkit. For installation and an overview, see the [README](README.md).

All commands assume you are in the project folder with the virtual environment active
(`source .venv/bin/activate`).

**Contents**

1. [Inputs](#1-inputs)
2. [Linking the taxonomy table to the trees](#2-linking-the-taxonomy-table-to-the-trees)
3. [Outputs](#3-outputs)
4. [Tanglegram recipes](#4-tanglegram-recipes)
5. [Single-tree recipes](#5-single-tree-recipes)
6. [Exporting to iTOL](#6-exporting-to-itol)
7. [Colours and accessibility](#7-colours-and-accessibility)
8. [Option reference](#8-option-reference)
9. [Troubleshooting](#9-troubleshooting)
10. [Using the toolkit as a library](#10-using-the-toolkit-as-a-library)

---

## 1. Inputs

| Input | Flag | Required | Notes |
|---|---|---|---|
| VAMPhyRE tree | `--tree-vamphyre` | Tanglegram | Newick |
| GTDB-Tk tree, relabelled to the VAMPhyRE leaf names | `--tree-gtdb` | Tanglegram | Newick |
| Single tree | `--arbol` | `generate_tree.py` | Any Newick: VAMPhyRE, GTDB-Tk, a consensus… |
| Taxonomy table | `--tabla` | Optional | `.xlsx`, `.tsv` or `.csv`, e.g. `gtdbtk.ar53.summary.tsv`. Enables everything taxonomic (colours, groups, highlighting, zoom). |

The table needs an ID column and a GTDB-style classification column
(`d__Archaea;p__…;c__…;o__…;f__…;g__…;s__…`). The defaults match the GTDB-Tk summary
(`--col-id user_genome`, `--col-tax classification`); change them only if your table uses
different headers.

To avoid repeating paths, define a shell variable:

```bash
BASE="--tree-vamphyre vamphyre.newick \
      --tree-gtdb gtdb_relabelled.newick \
      --tree-gtdb-original gtdbtk.unrooted.tree \
      --tabla gtdbtk.ar53.summary.tsv"
```

---

## 2. Linking the taxonomy table to the trees

GTDB-Tk tables use accession IDs (`GCF_…`), while VAMPhyRE trees often use sequence
names (`NZ_CP…_Species_…`). If the two cannot be matched, figures will have no colours
or group names. Provide **one** of the following:

1. `--tree-gtdb-original gtdbtk.unrooted.tree` (single tree: `--arbol-referencia`) — the
   original GTDB-Tk tree with table IDs and the same topology. The toolkit derives the mapping.
2. `--mapa-ids map.tsv` — two tab-separated columns: *table ID* ⇥ *leaf name in the tree*.
3. `--fasta-dir genomes/` — a folder of `.fna` files; maps file name → FASTA header.

The first run writes the mapping to `<prefix>_mapa_ids.tsv`; reuse it later with `--mapa-ids`
(this is the easiest way to feed taxonomy to `generate_tree.py`).

The web app reports coverage (e.g. *"With taxonomy: 0/537"*); the CLI prints a `[!]` warning.

---

## 3. Outputs

Output names derive from `--out-plot` (e.g. `results/x.png` → prefix `results/x`).

**Tanglegram (`generate_tanglegram.py`)**

| File | Content |
|---|---|
| `x.png` | Tanglegram (plus `.svg`/`.pdf` with `--formato-extra`) |
| `x_rebeldes.png` | Summary figure of rogue genomes and groups |
| `x_reporte.md` | Human-readable report: RF, rogues, least monophyletic groups |
| `x_hojas_rebeldes.csv` | Contribution of each genome to RF |
| `x_grupos.csv` | Monophyly / fragment count per family and genus |
| `x_nodos_vamphyre.csv`, `x_nodos_gtdb.csv` | Every node with ID, size and concordance |
| `x_mapa_ids.tsv` | ID mapping used |

**Single tree (`generate_tree.py`)**

| File | Content |
|---|---|
| `x.png` | Tree figure (plus `.svg`/`.pdf` with `--formato-extra`) |
| `x_arbol.nwk` | The tree as drawn (after pruning, rooting and support collapse) |
| `x_nodos.csv` | Internal nodes: ID, size, support (0–100), branch length, consensus taxon |
| `x_grupos.csv` | Monophyly of each taxon in this tree (requires a table) |
| `x_reporte.md` | Report |
| `x_mapa_ids.tsv` | ID mapping used (when a table is given) |
| `x_itol.*` | iTOL files — see [section 6](#6-exporting-to-itol) |

---

## 4. Tanglegram recipes

**Change the colour rank** — `domain`, `phylum`, `class`, `order`, `family`, `genus`,
`species` (Spanish names such as `genero` and the letters `d,p,c,o,f,g,s` also work):

```bash
python generate_tanglegram.py $BASE --color-rango order --out-plot results/order.png
python generate_tanglegram.py $BASE --color-rango genus --out-plot results/genus.png
```

**Highlight taxa** (everything else is dimmed):

```bash
python generate_tanglegram.py $BASE --color-rango genero \
    --resaltar g__Halorubrum,g__Haloquadratum --color-top 10 --out-plot results/highlight.png
```

A taxon outside the top-N is grey even when highlighted; raise `--color-top` if needed.
Taxa can be written with or without prefix (`g__Halorubrum` or `Halorubrum`), case-insensitive.

**Zoom into a group** (everything is recomputed on the subset; colours stay the same as in the
full figure):

```bash
python generate_tanglegram.py $BASE --color-rango genero --solo-taxa f__Haloferacaceae \
    --etiquetas-hojas --out-plot results/haloferacaceae.png
```

**Collapse clades**

```bash
# Genus-pure clades as triangles — best overview for hundreds of genomes
python generate_tanglegram.py $BASE --color-rango genero --colapsar taxonomia \
    --rango-colapso genero --out-plot results/collapsed_genus.png

# Collapse clades identical in both trees, to show only the conflict
python generate_tanglegram.py $BASE --colapsar concordantes --out-plot results/conflict.png

# The opposite: collapse the discordant part
python generate_tanglegram.py $BASE --colapsar discordantes --umbral 0.5 --out-plot results/x.png
```

Fine-tuning: `--min-colapso` / `--max-colapso` (clade size), `--colapsar-taxa` /
`--expandir-taxa`, and per node `--colapsar-nodos V12,G40` / `--expandir-nodos`
(node IDs are in the node CSVs, or drawn with `--mostrar-ids`).

**Reduce visual noise**

```bash
--ocultar-lineas g__Halorubrum   # do not draw these connecting lines
--grosor-linea 2.0               # thicker connecting lines
--escala longitud                # real branch lengths (default: cladogram)
--alto 40 --dpi 150              # height (inches) and resolution
--formato-extra svg,pdf          # vector output in addition to PNG
--max-rebeldes 0                 # skip the greedy rogue search (faster)
```

**One tree from the tanglegram run** — `--figura vamphyre|gtdb|individuales` draws a single
tree (or each tree in its own file) using the same leaf order as the tanglegram, so the
figures correspond. `concordantes`/`discordantes` collapse is still computed against the
other tree.

```bash
python generate_tanglegram.py $BASE --figura individuales --color-rango genero \
    --out-plot results/tree.png      # -> tree_vamphyre.png and tree_gtdb.png
```

---

## 5. Single-tree recipes

`generate_tree.py` needs only **one** Newick and accepts the same colour, collapse, zoom,
highlight and label options, plus layout, rooting and support options.

```bash
# GTDB-Tk, circular, coloured by genus; dots where support >= 90; nodes < 70 collapsed
python generate_tree.py --arbol gtdb.newick --tabla summary.tsv \
    --mapa-ids results/tanglegram_mapa_ids.tsv \
    --color-rango genero --raiz punto-medio --contraer-soporte 70 \
    --soporte punto --soporte-min 90 --out-plot results/gtdb_circular.png

# VAMPhyRE, rectangular, genus-pure clades as triangles
python generate_tree.py --arbol vamphyre.newick --tabla summary.tsv \
    --mapa-ids results/tanglegram_mapa_ids.tsv \
    --color-rango genero --colapsar taxonomia --rango-colapso genero \
    --forma rectangular --out-plot results/vamphyre.png
```

- `--forma circular|rectangular` (default `circular`) and `--arco 180..360` (the gap is at the top).
- `--raiz original|punto-medio|taxon` (+ `--raiz-taxon p__X`, which must form a clade). Changes the drawing only.
- `--soporte ninguno|numero|punto` and `--soporte-min` (0–100) control which support values are drawn.
- `--contraer-soporte U` turns nodes with support < U into polytomies.

In the web app, use the **Single tree** page; it can reuse the files uploaded on the
Tanglegram page.

---

## 6. Exporting to iTOL

Every `generate_tree.py` run also writes files for [iTOL](https://itol.embl.de), so the tree
can be refined interactively (fonts, layout, extra datasets) instead of editing an image.

| File | iTOL dataset type | What it contains |
|---|---|---|
| `x_itol.tree` | Tree | Newick with branch lengths and support (shown by iTOL as bootstrap) |
| `x_itol_colores.txt` | `TREE_COLORS` | Branch colour for each maximal pure clade; label colour for each leaf |
| `x_itol_franja.txt` | `DATASET_COLORSTRIP` | Colour strip by taxon, with legend |
| `x_itol_etiquetas.txt` | `LABELS` | GTDB species names as display labels |
| `x_itol_colapso.txt` | `COLLAPSE` | Clades collapsed in the figure |
| `x_itol.zip` | — | All of the above, for batch upload |

Colour, strip and label files are written only when a taxonomy table is given; the
collapse file only when something is collapsed.

**How to load it**

- *Option A (one step):* on itol.embl.de go to **Upload → Batch upload** and select `x_itol.zip`.
- *Option B:* upload `x_itol.tree`, then drag and drop each `.txt` file onto the open tree.

**What to expect**

- The tree reflects the final state of the run: pruning (`--solo-taxa`), rooting (`--raiz`)
  and support collapse (`--contraer-soporte`) are already applied.
- Highlighting (`--resaltar`) is preserved: non-highlighted taxa are written in the dimmed grey.
- Internal nodes are referenced as `leafA|leafB` (their last common ancestor), the identifier
  format documented by iTOL, so the files do not rely on internal node names.
- Each file is a separate dataset in iTOL and can be hidden, restyled or edited there.
- Use `--sin-itol` (web app: untick *Export for iTOL*) to skip these files.

---

## 7. Colours and accessibility

`--color-top N` (default 8) sets how many of the most frequent taxa at the chosen rank get
their own colour; the rest are grey ("other"). Genomes without a classification at that
rank are light grey ("unclassified").

| N | Palette | Colour-blind safe? | Use for |
|---|---|---|---|
| 1–8 | Validated categorical palette | **Yes** | Publication figures |
| 9–20 | + 12 extra hand-picked tones | No | Exploration; combine with `--resaltar` |
| > 20 | + generated tones (golden-angle hue spacing, varied lightness) | No — also hard to distinguish with normal vision | Exploration only |

When N is above 8 or above 20, the CLI logs a warning and the web app displays one.
The first 20 colours are fixed, so colours do not shift between runs or when zooming.

**Accessible alternatives to many colours**

- Colour at a higher rank (e.g. family instead of genus).
- Highlight a few taxa with `--resaltar`; the rest are dimmed.
- Collapse pure clades with `--colapsar taxonomia`; triangles are labelled with the taxon name.
- Zoom with `--solo-taxa` and colour within the subset.
- In iTOL, the taxonomy strip carries the taxon name for each leaf, which is readable
  regardless of colour perception.

---

## 8. Option reference

Flag names and their values are in Spanish (kept for backward compatibility); all help text and output are in English. *T* = `generate_tanglegram.py`, *S* = `generate_tree.py`.

### Inputs and output

| Flag | Tool | Meaning |
|---|---|---|
| `--tree-vamphyre`, `--tree-gtdb` | T | The two trees to compare |
| `--arbol` | S | The tree to draw |
| `--out-plot` | T, S | Output image path (PNG); defines the prefix for all other files |
| `--prefijo-salida` | T, S | Prefix for tables/report (default: `--out-plot` without extension) |
| `--tabla` | T, S | Taxonomy table (xlsx/tsv/csv) |
| `--col-id`, `--col-tax` | T, S | Table columns with genome ID and GTDB classification |
| `--tree-gtdb-original` | T | Original GTDB-Tk tree (table IDs), used to build the ID map |
| `--arbol-referencia` | S | Same purpose as above for a single tree |
| `--mapa-ids` | T, S | TSV: table ID ⇥ leaf name |
| `--fasta-dir` | T, S | Folder of `.fna` files for ID mapping |
| `--figura` | T | `tanglegram` (default), `vamphyre`, `gtdb`, `individuales` |
| `--titulo` | S | Figure title |
| `--formato-extra` | T, S | Extra figure formats: `svg,pdf` |
| `--sin-itol` | S | Do not write iTOL files |
| `--silencioso` | T, S | Show warnings only |

### Colour and selection

| Flag | Meaning |
|---|---|
| `--color-rango` | Rank used for colour: `domain`, `phylum`, `class`, `order`, `family`, `genus`, `species` (or Spanish `dominio`…`especie`, or `d,p,c,o,f,g,s`) |
| `--color-top` | Number of taxa with their own colour (see [section 7](#7-colours-and-accessibility)) |
| `--resaltar` | Comma-separated taxa to highlight (highlight) |
| `--solo-taxa` | Keep only these taxa (prune and recompute) |
| `--ocultar-lineas` | T only: taxa whose connecting lines are hidden |

### Collapse

| Flag | Meaning |
|---|---|
| `--colapsar` | `ninguno` (none), `taxonomia` (pure clades); T also: `concordantes` (concordant), `discordantes` (discordant) |
| `--umbral` | T: clade compatibility threshold (0–1) for concordant/discordant |
| `--min-colapso`, `--max-colapso` | Clade size range (leaves) eligible for collapse |
| `--rango-colapso` | Rank for `taxonomia` collapse (default: the colour rank) |
| `--colapsar-nodos`, `--expandir-nodos` | Node IDs always collapsed / kept open |
| `--colapsar-taxa`, `--expandir-taxa` | Taxa always collapsed / kept open |

### Tree (single tree)

| Flag | Meaning |
|---|---|
| `--raiz` | Root: `original`, `punto-medio` (midpoint), `taxon` (outgroup) |
| `--raiz-taxon` | Outgroup taxon when `--raiz taxon` (T: common root for both trees) |
| `--contraer-soporte` | Collapse nodes with support below this value (0–100) into polytomies |
| `--orden` | `ladderize` (by clade size) or `original` (file order) |

### Appearance

| Flag | Meaning |
|---|---|
| `--forma` | `circular` or `rectangular` |
| `--arco` | Circular: degrees spanned (180–360) |
| `--soporte`, `--soporte-min` | Draw support as `numero` (number) / `punto` (dot) / `ninguno`, above a threshold |
| `--escala` | `cladograma` (equal branches) or `longitud` (branch lengths) |
| `--etiquetas-hojas`, `--texto-hoja` | Write leaf names; text is `especie` (species) or `nombre` (raw name) |
| `--sin-etiquetas-grupos`, `--min-etiqueta` | Disable / size threshold for clade labels |
| `--mostrar-ids`, `--min-id` | Draw node IDs for nodes above a size |
| `--ancho-arbol`, `--ancho-centro`, `--alto` | Tree width (or radius), T centre width, figure height (inches) |
| `--dpi` | Resolution |
| `--grosor-rama`, `--grosor-linea`, `--escala-triangulo` | Branch width, T connecting-line width, triangle size |
| `--sin-desenredar` | T: do not reorder branches to reduce crossings |

### Analysis and tables

| Flag | Meaning |
|---|---|
| `--max-rebeldes` | T: number of greedy rogue-removal steps (0 = skip) |
| `--rango-grupos` | Ranks evaluated in the monophyly table (default `familia,genero`) |
| `--min-grupo` | Minimum group size evaluated |

---

## 9. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| No colours or group names | The table is not linked to the trees. Add `--tree-gtdb-original`, `--mapa-ids` or `--fasta-dir` ([section 2](#2-linking-the-taxonomy-table-to-the-trees)). |
| A highlighted taxon is grey | It is outside the top-N coloured taxa. Increase `--color-top`. |
| Leaf names are unreadable | Too many leaves. Combine `--etiquetas-hojas` with `--solo-taxa`. |
| `--raiz-taxon` ignored | The taxon does not form a clade in that tree; the original root is kept (a warning is printed). |
| Colours look alike | More than 8 (or 20) colours. See [section 7](#7-colours-and-accessibility). |
| iTOL shows no colours | Make sure the annotation `.txt` files were loaded with the matching `_itol.tree`; leaf names must match exactly. |

---

## 10. Using the toolkit as a library

```python
from tanglegram_utils.models import Entradas, OpcionesAnalisis, Salida
from tanglegram_utils.pipeline import ejecutar

datos, an = ejecutar(
    Entradas("vamphyre.newick", "gtdb.newick", tabla="summary.tsv",
             tree_gtdb_original="gtdbtk.unrooted.tree"),
    OpcionesAnalisis(color_rango="genero", solo_taxa="f__Haloferacaceae"),
    salida=Salida(out_plot="results/haloferacaceae.png"),
)
print(an.rf, an.nrf)
```

Single tree:

```python
from tanglegram_utils.models import EntradaArbol, OpcionesArbol, Salida
from tanglegram_utils.pipeline_arbol import ejecutar_arbol
from tanglegram_utils.plotting import EstiloFigura

an = ejecutar_arbol(
    EntradaArbol("gtdb.newick", tabla="summary.tsv", mapa_ids="results/tanglegram_mapa_ids.tsv"),
    OpcionesArbol(color_rango="genero", color_top=8),
    estilo=EstiloFigura(forma="circular", soporte="punto", soporte_min=90),
    salida=Salida(out_plot="results/gtdb.png"),
)
```

`tanglegram_utils.servicio` offers helpers to run in a directory, capture the log and
package results as a ZIP (used by the web app).
