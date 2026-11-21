"""Export to iTOL (https://itol.embl.de): tree + editable annotation files.

iTOL does not import images: it reads a Newick tree and, separately, plain-text
annotation files that are dragged onto the tree (or uploaded together in a ZIP via
"Batch upload"). Each layer (colours, taxonomy strip, labels, collapsed clades) thus
becomes a dataset that can be edited, hidden or restyled inside iTOL.

Files written (with `<prefijo>`):

- `_itol.tree`            Newick with branch lengths and support (shown by iTOL as bootstrap).
- `_itol_colores.txt`     TREE_COLORS: branch colour of each pure clade and label colour of each leaf.
- `_itol_franja.txt`      DATASET_COLORSTRIP: colour strip by taxon, with legend.
- `_itol_etiquetas.txt`   LABELS: readable name (GTDB species) instead of the genome ID.
- `_itol_colapso.txt`     COLLAPSE: the same collapsed clades as in the figure.
- `_itol.zip`             all of the above, ready for iTOL "Batch upload".

Internal nodes are identified as `leafA|leafB` (last common ancestor), the format
documented by iTOL, which does not depend on internal node names in the Newick."""
import logging
import os
import zipfile

from .constants import GRIS_ATENUADO, GRIS_OTROS, GRIS_SIN, SIN
from .taxonomy import nombre_rango, sin_prefijo
from .tree import a_newick, hojas_de, preorden

log = logging.getLogger(__name__)

SUFIJOS_ITOL = ("_itol.tree", "_itol_colores.txt", "_itol_franja.txt", "_itol_etiquetas.txt",
                "_itol_colapso.txt", "_itol.zip")


def _limpio(txt):
    """Text that is safe inside a tab-separated field."""
    return str(txt).replace("\t", " ").replace("\n", " ").strip()


def _primera_hoja(n):
    while not n.es_hoja:
        n = n.hijos[0]
    return n


def _ultima_hoja(n):
    while not n.es_hoja:
        n = n.hijos[-1]
    return n


def id_nodo_itol(n):
    """iTOL ID for a node: the leaf name for a leaf; 'A|B' (common ancestor) for an internal node."""
    if n.es_hoja:
        return n.nombre
    return f"{_primera_hoja(n.hijos[0]).nombre}|{_ultima_hoja(n.hijos[-1]).nombre}"


def _cabecera(tipo, etiqueta=None, color=None):
    lineas = [tipo, "SEPARATOR TAB"]
    if etiqueta:
        lineas.append(f"DATASET_LABEL\t{_limpio(etiqueta)}")
    if color:
        lineas.append(f"COLOR\t{color}")
    return lineas


def _escribir(ruta, lineas):
    with open(ruta, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lineas) + "\n")
    return ruta


def _color_hoja(esquema, i, hl, hay_res):
    if hay_res and not hl[i]:
        return GRIS_ATENUADO
    return esquema.col_cat[esquema.cat[i]]


def _clados_puros_maximos(raiz, esquema):
    """Internal nodes whose subtree is a single taxon and whose parent is not."""
    out, pila = [], [raiz]
    while pila:
        n = pila.pop()
        if n.es_hoja:
            continue
        if n.n >= 2 and esquema.categoria_pura(n.hojas) is not None:
            out.append(n)
        else:
            pila.extend(n.hijos)
    return out


def exportar_itol(raiz, datos, prefijo, colapsados=(), hl=None, hay_res=False):
    """Write the tree and its annotations for iTOL. Returns the list of files created."""
    esquema, idx = datos.esquema, datos.idx
    hojas = hojas_de(raiz)
    creados = []

    ruta_arbol = prefijo + "_itol.tree"
    with open(ruta_arbol, "w", encoding="utf-8", newline="\n") as f:
        f.write(a_newick(raiz) + "\n")
    creados.append(ruta_arbol)

    if esquema.hay_tax:
        rango = nombre_rango(esquema.rango)
        # 1) branch colours (pure clades) and leaf label colours
        lineas = _cabecera("TREE_COLORS") + ["DATA"]
        for n in _clados_puros_maximos(raiz, esquema):
            cat = esquema.categoria_pura(n.hojas)
            col = esquema.col_cat[cat]
            if col in (GRIS_OTROS, GRIS_SIN):
                continue
            if hay_res and not hl[n.hojas].any():
                col = GRIS_ATENUADO
            lineas.append(f"{id_nodo_itol(n)}\tclade\t{col}\tnormal\t1")
        for h in hojas:
            i = idx[h.nombre]
            lineas.append(f"{h.nombre}\tlabel\t{_color_hoja(esquema, i, hl, hay_res)}\tnormal\t1")
        creados.append(_escribir(prefijo + "_itol_colores.txt", lineas))

        # 2) colour strip by taxon, with legend (only taxa present in this tree)
        presentes = {esquema.cat[idx[h.nombre]] for h in hojas}
        leyenda = [c for c in esquema.orden if c in presentes and esquema.col_cat[c] != GRIS_OTROS]
        etiquetas = [sin_prefijo(c) for c in leyenda]
        colores = [esquema.col_cat[c] for c in leyenda]
        if any(esquema.col_cat.get(c) == GRIS_OTROS for c in presentes):
            etiquetas.append("other")
            colores.append(GRIS_OTROS)
        if SIN in presentes:
            etiquetas.append("unclassified")
            colores.append(GRIS_SIN)
        lineas = _cabecera("DATASET_COLORSTRIP", f"Taxonomy ({rango})", colores[0] if colores else "#000000")
        lineas += [f"LEGEND_TITLE\t{_limpio(rango.capitalize())}",
                   "LEGEND_SHAPES\t" + "\t".join(["1"] * len(colores)),
                   "LEGEND_COLORS\t" + "\t".join(colores),
                   "LEGEND_LABELS\t" + "\t".join(_limpio(e) for e in etiquetas),
                   "STRIP_WIDTH\t25", "MARGIN\t5", "SHOW_INTERNAL\t0", "DATA"]
        for h in hojas:
            i = idx[h.nombre]
            cat = esquema.cat[i]
            txt = "unclassified" if cat == SIN else sin_prefijo(cat)
            lineas.append(f"{h.nombre}\t{_color_hoja(esquema, i, hl, hay_res)}\t{_limpio(txt)}")
        creados.append(_escribir(prefijo + "_itol_franja.txt", lineas))

        # 3) readable names (GTDB species) instead of the genome ID
        lineas = _cabecera("LABELS") + ["DATA"]
        for h in hojas:
            lineas.append(f"{h.nombre}\t{_limpio(datos.corto(idx[h.nombre]))}")
        creados.append(_escribir(prefijo + "_itol_etiquetas.txt", lineas))

    # 4) collapsed clades (same as in the figure)
    vivos = {id(n) for n in preorden(raiz)}
    nodos = [n for n in colapsados if id(n) in vivos and not n.es_hoja]
    if nodos:
        lineas = ["COLLAPSE", "DATA"] + [id_nodo_itol(n) for n in nodos]
        creados.append(_escribir(prefijo + "_itol_colapso.txt", lineas))

    # ZIP for "Batch upload" (.tree file + .txt annotations)
    ruta_zip = prefijo + "_itol.zip"
    with zipfile.ZipFile(ruta_zip, "w", zipfile.ZIP_DEFLATED) as z:
        for ruta in creados:
            z.write(ruta, arcname=os.path.basename(ruta))
    creados.append(ruta_zip)
    log.info(f"[+] iTOL: {prefijo}_itol.tree + {len(creados) - 2} annotation file(s) and {os.path.basename(ruta_zip)}")
    return creados
