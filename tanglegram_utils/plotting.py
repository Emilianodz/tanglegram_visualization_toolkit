"""Tanglegram drawing (matplotlib): trees, triangles, colour strips,
group labels, connecting lines and legend."""
import collections
import logging
import math
import os
from dataclasses import dataclass

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle

from .constants import (GRIS_ATENUADO, GRIS_MIXTO, GRIS_OTROS, GRIS_RAMA, GRIS_SIN,
                        LETRA_A_NOMBRE, SIN, TINTA, TINTA_SUAVE)
from .layout import escalar_y, layout_arbol
from .ordering import posiciones
from .support import escala_soporte
from .taxonomy import etiqueta_consenso, sin_prefijo

log = logging.getLogger(__name__)

NOMBRES_PANEL = {"vamphyre": "VAMPhyRE", "gtdb": "GTDB-Tk"}
# Display names for the collapse modes (CLI values are in Spanish)
NOMBRES_COLAPSO = {"ninguno": "none", "taxonomia": "taxonomy", "concordantes": "concordant",
                   "discordantes": "discordant"}


@dataclass
class EstiloFigura:
    """Visual parameters (everything that is not data)."""
    escala: str = "cladograma"           # cladograma (cladogram) | longitud (branch lengths)
    etiquetas_hojas: bool = False
    texto_hoja: str = "especie"          # especie (species) | nombre (raw name)
    etiquetas_grupos: bool = True
    min_etiqueta: int = 4
    mostrar_ids: bool = False
    min_id: int = 15
    ancho_arbol: float = None            # inches; None = automatic (4.5 in a tanglegram, 8 for a single tree)
    ancho_centro: float = 5.0
    alto: float = None
    dpi: int = 150
    grosor_linea: float = 1.3
    grosor_rama: float = 1.1
    escala_triangulo: float = 1.0
    forma: str = "rectangular"           # rectangular | circular
    arco: float = 360.0                  # circular: degrees spanned by the tree (the gap is at the top)
    soporte: str = "ninguno"             # ninguno (none) | numero (number) | punto (dot): node support display
    soporte_min: float = 70.0            # only support >= this value is shown (0-100 scale)


@dataclass
class ContenidoFigura:
    """What is drawn: trees, collapses, colours and the subtitle summary."""
    rV: object
    rG: object
    colV: set
    colG: set
    esquema: object
    tax: list
    nombres: list
    hl: np.ndarray
    hay_res: bool
    vis_linea: np.ndarray
    N: int
    N0: int
    rf: int
    nrf: float
    cr0: int
    cr1: int
    modo_colapso: str = "ninguno"
    conservar: set = None       # indices of displayed leaves (None = all); filters the legend
    titulo: str = None          # single tree: custom title (default: the panel name)
    subtitulo: str = None       # single tree: text below the title


class DibujanteTanglegram:
    """Build the figure. Usage: DibujanteTanglegram(contenido, estilo).dibujar().

    `panel`: 'ambos' (tanglegram), 'vamphyre' or 'gtdb' (a single tree, with the
    same colour, collapse, label and highlight options)."""

    # fixed layout margins and widths (inches)
    MARGEN, STRIP_W, PAD = 0.35, 0.14, 0.05
    W_HOJA, W_GRUPO = 2.3, 2.4
    ANCHO_PAR, ANCHO_SOLO = 4.5, 8.0        # automatic width of each tree (inches)
    SUP, INF = 4.5, 1.5

    def __init__(self, contenido, estilo, panel="ambos"):
        if panel not in ("ambos", "vamphyre", "gtdb"):
            raise ValueError(f"unknown panel: {panel}")
        self.panel = panel
        self.solo = panel != "ambos"
        self.c = contenido
        self.e = estilo
        self.es = contenido.esquema
        self._cache_rama = {}
        self._calcular_geometria()

    # ------------------------------------------------------------------ geometry
    def _calcular_geometria(self):
        c, e = self.c, self.e
        w = e.ancho_arbol or (self.ANCHO_SOLO if self.solo else self.ANCHO_PAR)
        self.w_hoja = self.W_HOJA if e.etiquetas_hojas else 0.0
        self.usar_grupos = self.es.hay_tax and e.etiquetas_grupos
        self.w_grupo = self.W_GRUPO if self.usar_grupos else 0.0
        gutter = 0.06 + self.STRIP_W + self.PAD + ((self.w_hoja + self.PAD) if self.w_hoja else 0.0)
        gutter += (0.07 + 0.06 + self.w_grupo) if self.usar_grupos else 0.0
        gutter += 0.1
        if self.solo:
            self.W = 2 * self.MARGEN + w + gutter
            self.xL_raiz = self.MARGEN
            self.xL_tip = self.xL_raiz + w
            self.bL = self._bandas(self.xL_tip, +1)
            raiz, col = self._arbol_solo()
            self.lay1 = layout_arbol(raiz, col, e.escala, "izq", self.xL_raiz, w, e.escala_triangulo)
            self.H = self.lay1["filas"]
        else:
            self.W = 2 * (self.MARGEN + w + gutter) + e.ancho_centro
            self.xL_raiz, self.xR_raiz = self.MARGEN, self.W - self.MARGEN
            self.xL_tip, self.xR_tip = self.xL_raiz + w, self.xR_raiz - w
            self.bL = self._bandas(self.xL_tip, +1)
            self.bR = self._bandas(self.xR_tip, -1)

            self.layV = layout_arbol(c.rV, c.colV, e.escala, "izq", self.xL_raiz, w, e.escala_triangulo)
            self.layG = layout_arbol(c.rG, c.colG, e.escala, "der", self.xR_raiz, w, e.escala_triangulo)
            self.H = max(self.layV["filas"], self.layG["filas"])
            escalar_y(self.layV, self.H / self.layV["filas"])
            escalar_y(self.layG, self.H / self.layG["filas"])

        self.alto = e.alto or min(110.0, max(9.0, self.H * 0.11 + 4.0))
        self.bottom_frac = min(0.2, 1.7 / self.alto)
        self.frac_ax = 1.0 - self.bottom_frac - 0.004
        self.pitch_pt = self.alto * self.frac_ax * 72.0 / (self.H + self.SUP + self.INF)
        self.fs_hoja = float(np.clip(0.78 * self.pitch_pt, 2.5, 9))
        self.fs_grupo = float(np.clip(1.05 * self.pitch_pt, 5, 11))

    def _arbol_solo(self):
        """(root, collapsed) of the single tree being drawn."""
        c = self.c
        return (c.rV, c.colV) if self.panel == "vamphyre" else (c.rG, c.colG)

    def _bandas(self, tip, s):
        """x positions of the colour strip, labels and line edge (s=+1 left, -1 right)."""
        d = {}
        c = tip + s * 0.06
        d["strip"] = (c, c + s * self.STRIP_W)
        c += s * (self.STRIP_W + self.PAD)
        d["hoja_x"] = c
        if self.w_hoja:
            c += s * (self.w_hoja + self.PAD)
        if self.usar_grupos:
            d["gbar"] = (c, c + s * 0.07)
            d["gtxt_x"] = c + s * 0.13
            c = d["gtxt_x"] + s * self.w_grupo
        d["link"] = c + s * 0.1
        return d

    # ------------------------------------------------------------------ colours
    def _color_hoja(self, i):
        if not self.es.hay_tax:
            # no taxonomy: in the tanglegram, colour links each leaf to its counterpart;
            # a single tree has no counterpart, so it is drawn in neutral grey
            return GRIS_RAMA if self.solo else self._cmap(self._pos_v[i] / max(self.c.N, 1))
        return self.es.col_cat[self.es.cat[i]]

    def _color_rama(self, h):
        if h in self._cache_rama:
            return self._cache_rama[h]
        c, es = self.c, self.es
        if not es.hay_tax:
            col = GRIS_RAMA
        else:
            cat = es.categoria_pura(h.hojas)
            if cat is not None:
                col = es.col_cat[cat]
                if c.hay_res and not c.hl[h.hojas].all():
                    col = GRIS_ATENUADO
            else:
                col = GRIS_RAMA if (not c.hay_res or c.hl[h.hojas].any()) else GRIS_ATENUADO
        self._cache_rama[h] = col
        return col

    # ------------------------------------------------------------------ drawing
    def dibujar(self):
        """Return the matplotlib figure, ready to save."""
        log.info("[-] Drawing tanglegram..." if not self.solo else f"[-] Drawing tree {self.c.titulo or NOMBRES_PANEL[self.panel]}...")
        self.fig = plt.figure(figsize=(self.W, self.alto), facecolor="white")
        self.ax = self.fig.add_axes([0.0, self.bottom_frac, 1.0, self.frac_ax])
        self.ax.set_xlim(0, self.W)
        self.ax.set_ylim(self.H + self.INF, -self.SUP)
        self.ax.axis("off")
        self._cmap = mpl.colormaps["turbo"]
        self._pos_v = posiciones(self.c.rV, self.c.N0)

        if self.solo:
            raiz, col = self._arbol_solo()
            self._dibujar_arbol(self.lay1, raiz, col, self.bL, +1)
        else:
            self._dibujar_arbol(self.layV, self.c.rV, self.c.colV, self.bL, +1)
            self._dibujar_arbol(self.layG, self.c.rG, self.c.colG, self.bR, -1)
            self._dibujar_enlaces()
        self._titulos()
        self._leyenda()
        return self.fig

    def _dibujar_arbol(self, lay, raiz, colapsados, bandas, s):
        c, e, es, ax = self.c, self.e, self.es, self.ax
        segs, cols = [], []
        for x1, y1, x2, y2, h in lay["segmentos"]:
            segs.append([(x1, y1), (x2, y2)])
            cols.append(self._color_rama(h))
        ax.add_collection(LineCollection(segs, colors=cols, linewidths=e.grosor_rama,
                                         capstyle="round", zorder=2))
        # triangles + colour strip at the tip
        a0, a1 = bandas["strip"]
        for u in lay["unidades"]:
            n = u["nodo"]
            if u["tri"]:
                cat = es.categoria_pura(n.hojas)
                col = es.col_cat[cat] if cat is not None else GRIS_RAMA
                if c.hay_res and not c.hl[n.hojas].any():
                    col = GRIS_ATENUADO
                hh = u["h"] * 0.46
                ax.add_patch(Polygon([(u["x_apice"], u["y"]), (u["x_base"], u["y"] - hh), (u["x_base"], u["y"] + hh)],
                                     closed=True, facecolor=col, edgecolor=col, alpha=0.9, lw=0.6, zorder=3))
                fc, alt = col, u["h"]
            else:
                i = int(n.hojas[0])
                fc = self._color_hoja(i) if (not c.hay_res or c.hl[i]) else GRIS_ATENUADO
                alt = 1.0
            ax.add_patch(Rectangle((min(a0, a1), u["y"] - alt / 2 + 0.06), abs(a1 - a0), alt - 0.12,
                                   facecolor=fc, edgecolor="none", zorder=3))
            if e.etiquetas_hojas and not u["tri"]:
                i = int(n.hojas[0])
                sp_ = c.tax[i].get("s")
                txt = sin_prefijo(sp_) if (e.texto_hoja == "especie" and sp_) else c.nombres[i]
                ax.text(bandas["hoja_x"], u["y"], txt[:48], fontsize=self.fs_hoja, ha="left" if s > 0 else "right",
                        va="center", color=TINTA if (not c.hay_res or c.hl[i]) else GRIS_ATENUADO, zorder=4)
        if self.usar_grupos:
            self._rotular_grupos(lay, raiz, colapsados, bandas, s)
        if e.soporte != "ninguno":
            self._dibujar_soporte(lay, raiz, colapsados, s)
        if e.mostrar_ids:
            for n in lay["y"]:
                if (not n.es_hoja) and n.id and n.n >= e.min_id and n not in colapsados:
                    ax.text(lay["x"][n] - s * 0.03, lay["y"][n] - 0.35, n.id, fontsize=max(3.0, self.pitch_pt * 0.5),
                            color=GRIS_MIXTO, ha="right" if s > 0 else "left", va="bottom", zorder=4)

    def _items_grupos(self, lay, raiz, colapsados):
        """Group labels: [(priority, n, y0, y1, text, colour)] for each triangle and pure clade."""
        c, e, es = self.c, self.e, self.es
        items = []
        for u in lay["unidades"]:
            if u["tri"]:
                n = u["nodo"]
                cat = es.categoria_pura(n.hojas)
                if cat is not None:
                    txt, col = sin_prefijo(cat), es.col_cat[cat]
                else:
                    lab, _ = etiqueta_consenso(n.hojas, c.tax)
                    txt, col = "≈ " + sin_prefijo(lab), GRIS_RAMA
                items.append((1, n.n, lay["y0"][n], lay["y1"][n], f"{txt} ({n.n})", col))
        pila = [raiz]
        while pila:
            n = pila.pop()
            if n.es_hoja or n in colapsados or n not in lay["y"]:
                continue
            if n.padre is not None and n.n >= e.min_etiqueta:
                cat = es.categoria_pura(n.hojas)
                if cat is not None:
                    items.append((0, n.n, lay["y0"][n], lay["y1"][n], f"{sin_prefijo(cat)} ({n.n})", es.col_cat[cat]))
                    continue
            pila.extend(n.hijos)
        items.sort(key=lambda t: (-t[0], -t[1]))
        return items

    def _rotular_grupos(self, lay, raiz, colapsados, bandas, s):
        """Bar + text for each triangle and pure clade (without overlapping labels)."""
        ax = self.ax
        items = self._items_grupos(lay, raiz, colapsados)
        colocados, alto_txt = [], self.fs_grupo / self.pitch_pt * 1.15
        g0, g1 = bandas["gbar"]
        for _, _nn, ya, yb, txt, col in items:
            ax.add_patch(Rectangle((min(g0, g1), ya + 0.1), abs(g1 - g0), max(yb - ya - 0.2, 0.2),
                                   facecolor=col, edgecolor="none", zorder=3))
            yc = (ya + yb) / 2
            if all(abs(yc - yp) >= alto_txt for yp in colocados):
                colocados.append(yc)
                ax.text(bandas["gtxt_x"], yc, txt, fontsize=self.fs_grupo, ha="left" if s > 0 else "right",
                        va="center", color=TINTA, zorder=4)

    def _pt(self, x, y):
        """Map layout space to canvas space (identity here; polar in the circular tree)."""
        return x, y

    def _nodos_con_soporte(self, lay, raiz, colapsados):
        """[(node, value 0-100)] for visible internal nodes with support >= soporte_min."""
        f = escala_soporte(raiz)
        out = []
        for n in lay["y"]:
            if n.es_hoja or n.padre is None or n in colapsados or n.soporte is None:
                continue
            v = n.soporte * f
            if v >= self.e.soporte_min:
                out.append((n, v))
        return out

    def _dibujar_soporte(self, lay, raiz, colapsados, s):
        """Support values: a number next to the node or a dot on it."""
        ax, e = self.ax, self.e
        nodos = self._nodos_con_soporte(lay, raiz, colapsados)
        if not nodos:
            return
        if e.soporte == "punto":
            xy = np.array([self._pt(lay["x"][n], lay["y"][n]) for n, _ in nodos])
            tam = float(np.clip(0.5 * self.pitch_pt, 1.5, 7.0)) ** 2
            ax.scatter(xy[:, 0], xy[:, 1], s=tam, c=TINTA, linewidths=0, zorder=5)
        else:
            for n, v in nodos:
                x, y = self._pt(lay["x"][n] - s * 0.03, lay["y"][n] - 0.18)
                ax.text(x, y, f"{v:.0f}", fontsize=max(3.0, self.pitch_pt * 0.5), color=GRIS_MIXTO,
                        ha="right" if s > 0 else "left", va="bottom", zorder=5)

    def _dibujar_enlaces(self):
        """Connecting lines: one per pair of units (leaf/triangle) sharing genomes."""
        c, e, es = self.c, self.e, self.es
        destino = np.full(c.N0, -1)
        for j, u in enumerate(self.layG["unidades"]):
            destino[u["hojas"]] = j
        enlaces = []
        for u in self.layV["unidades"]:
            grupos = collections.defaultdict(list)
            for i in u["hojas"]:
                if c.vis_linea[i] and destino[i] >= 0:
                    grupos[destino[i]].append(int(i))
            for j, lst in grupos.items():
                enlaces.append((u["y"], self.layG["unidades"][j]["y"], lst))
        buckets = {1: [], 2: [], 3: []}
        for yv, yg, lst in enlaces:
            k = len(lst)
            cnt = collections.Counter(es.cat[i] for i in lst)
            top, c_top = cnt.most_common(1)[0]
            if not es.hay_tax:
                col = self._color_hoja(lst[0])
            elif c_top / k >= 0.6:
                col = es.col_cat[top]
            else:
                col = GRIS_MIXTO
            lw = e.grosor_linea * (1.0 if k == 1 else min(1.0 + 0.9 * math.sqrt(k - 1), 11.0))
            es_hl = bool(c.hl[lst].mean() > 0.5)
            if c.hay_res:
                if es_hl:
                    lw, alpha, z = lw * 1.8, 0.92, 3
                else:
                    col, lw, alpha, z = GRIS_ATENUADO, lw * 0.6, 0.45, 1
            else:
                alpha = 0.78 if col not in (GRIS_OTROS, GRIS_SIN) else 0.5
                z = 2
            buckets[z].append(([(self.bL["link"], yv), (self.bR["link"], yg)], col, lw, alpha))
        for z, lst in buckets.items():
            if not lst:
                continue
            segs = [a for a, _, _, _ in lst]
            cols = [to_rgba(col, al) for _, col, _, al in lst]
            self.ax.add_collection(LineCollection(segs, colors=cols, linewidths=[lw for _, _, lw, _ in lst],
                                                  capstyle="round", zorder=z))

    def _titulos(self):
        c, ax = self.c, self.ax
        if self.solo:
            ax.text((self.xL_raiz + self.xL_tip) / 2, -2.4, c.titulo or NOMBRES_PANEL[self.panel], fontsize=22,
                    ha="center", weight="bold", color=TINTA)
            sub = c.subtitulo or f"{c.N} genomes"
            if (c.colV if self.panel == "vamphyre" else c.colG):
                sub += f" · collapsed clades: {NOMBRES_COLAPSO.get(c.modo_colapso, c.modo_colapso)}"
            ax.text(self.W / 2, -0.7, sub, fontsize=11, ha="center", va="center", color=TINTA_SUAVE)
            return
        ax.text((self.xL_raiz + self.xL_tip) / 2, -2.4, "VAMPhyRE", fontsize=22, ha="center", weight="bold", color=TINTA)
        ax.text((self.xR_raiz + self.xR_tip) / 2, -2.4, "GTDB-Tk", fontsize=22, ha="center", weight="bold", color=TINTA)
        sub = f"{c.N} genomes · RF = {c.rf} (normalised {c.nrf:.2f}) · line crossings {c.cr0} → {c.cr1}"
        if c.colV or c.colG:
            sub += f" · collapsed clades: {NOMBRES_COLAPSO.get(c.modo_colapso, c.modo_colapso)}"
        ax.text(self.W / 2, -2.4, sub, fontsize=11, ha="center", va="center", color=TINTA_SUAVE)

    def _leyenda(self):
        es = self.es
        if not es.hay_tax:
            return
        # The legend lists only what is visible (with --solo-taxa: the subset); colours
        # remain those of the full dataset.
        vistos = range(len(es.cat)) if self.c.conservar is None else sorted(self.c.conservar)
        cont = collections.Counter(es.cat[i] for i in vistos if es.cat[i] != SIN)
        orden = sorted(cont, key=lambda c: (-cont[c], c))
        manos = [Line2D([0], [0], color=es.col_cat[c], lw=6, solid_capstyle="butt",
                        label=f"{sin_prefijo(c)} ({cont[c]})") for c in orden if es.col_cat[c] != GRIS_OTROS]
        otros = [c for c in orden if es.col_cat[c] == GRIS_OTROS]
        if otros:
            manos.append(Line2D([0], [0], color=GRIS_OTROS, lw=6, solid_capstyle="butt",
                                label=f"other: {len(otros)} taxa ({sum(cont[c] for c in otros)})"))
        n_sin = sum(1 for i in vistos if es.cat[i] == SIN)
        if n_sin:
            manos.append(Line2D([0], [0], color=GRIS_SIN, lw=6, solid_capstyle="butt",
                                label=f"unclassified ({n_sin})"))
        self.fig.legend(handles=manos, loc="lower center", ncol=min(6, len(manos)), frameon=False,
                        fontsize=12, title=f"Colour by {LETRA_A_NOMBRE[es.rango]}", title_fontsize=13,
                        bbox_to_anchor=(0.5, 0.0))


def guardar_figura(fig, ruta_principal, dpi, formato_extra=None):
    """Save in the format of `ruta_principal` and, optionally, in others (svg,pdf)."""
    formatos = [os.path.splitext(ruta_principal)[1].lstrip(".") or "png"]
    if formato_extra:
        formatos += [f.strip() for f in formato_extra.split(",") if f.strip()]
    base = os.path.splitext(ruta_principal)[0]
    for fmt in formatos:
        ruta = ruta_principal if fmt == formatos[0] else f"{base}.{fmt}"
        fig.savefig(ruta, dpi=dpi, bbox_inches="tight", facecolor="white")
        log.info(f"[+] Figure: {ruta}")
    plt.close(fig)
