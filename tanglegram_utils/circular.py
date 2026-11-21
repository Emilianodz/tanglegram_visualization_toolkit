"""Single tree in circular layout (matplotlib).

Reuses the rectangular tree geometry (`layout_arbol`) and maps it to polar
coordinates: distance from the root becomes the radius and vertical position
becomes the angle. Colours, legend and group logic are inherited from
`DibujanteTanglegram`, so colour, highlighting, collapse and support behave
exactly as in the rectangular figure."""
import logging
import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon, Wedge

from .constants import GRIS_ATENUADO, GRIS_MIXTO, GRIS_RAMA, TINTA, TINTA_SUAVE
from .layout import layout_arbol
from .plotting import NOMBRES_COLAPSO, NOMBRES_PANEL, DibujanteTanglegram
from .taxonomy import sin_prefijo
from .tree import hojas_de as raiz_hojas

log = logging.getLogger(__name__)


class DibujanteCircular(DibujanteTanglegram):
    """Usage: DibujanteCircular(contenido, estilo, panel).dibujar() (panel: 'vamphyre' or 'gtdb')."""

    ANILLO_W, PAD, BARRA_W = 0.16, 0.05, 0.07
    TITULO_H, LEYENDA_H = 1.0, 1.7
    PASO_ARCO, PASO_ARCO_ETIQ = 0.07, 0.11      # inches of circumference per leaf (without / with names)

    def __init__(self, contenido, estilo, panel="vamphyre"):
        if panel == "ambos":
            raise ValueError("The circular layout draws a single tree (panel 'vamphyre' or 'gtdb').")
        super().__init__(contenido, estilo, panel)

    # ------------------------------------------------------------------ geometry
    def _calcular_geometria(self):
        e = self.e
        self.arco_rad = math.radians(min(max(e.arco, 30.0), 360.0))
        raiz, col = self._arbol_solo()
        # first with radius 1 to get the number of rows (leaves + triangles)
        filas = layout_arbol(raiz, col, e.escala, "izq", 0.0, 1.0, e.escala_triangulo)["filas"]
        paso = self.PASO_ARCO_ETIQ if e.etiquetas_hojas else self.PASO_ARCO
        self.R = e.ancho_arbol or float(np.clip(filas * paso / self.arco_rad, 3.0, 16.0))
        self.lay1 = layout_arbol(raiz, col, e.escala, "izq", 0.0, self.R, e.escala_triangulo)
        self.H = self.lay1["filas"]

        # points per row at the tree tips (font sizes)
        self.pitch_pt = max(self.R * self.arco_rad / self.H * 72.0, 1.0)
        self.fs_hoja = float(np.clip(0.78 * self.pitch_pt, 2.5, 9))
        self.fs_grupo = float(np.clip(1.05 * self.pitch_pt, 5, 11))

        # text widths (estimated from length) to size the figure
        def ancho_txt(nchar, fs):
            return nchar * fs * 0.62 / 72.0 + 0.1

        self.usar_grupos = self.es.hay_tax and e.etiquetas_grupos
        self.w_hoja = 0.0
        if e.etiquetas_hojas:
            self.w_hoja = ancho_txt(max(len(self._texto_hoja(int(h.hojas[0]))) for h in raiz_hojas(raiz)), self.fs_hoja)
        r = self.R + 0.06
        self.r_strip = (r, r + self.ANILLO_W)
        r += self.ANILLO_W + self.PAD
        self.r_hoja = r
        if self.w_hoja:
            r += self.w_hoja + self.PAD
        self.r_gbar = (r, r + self.BARRA_W)
        self.r_gtxt = r + self.BARRA_W + 0.06
        self.w_grupo = 0.0
        if self.usar_grupos:
            items = self._items_grupos(self.lay1, raiz, col)
            self.w_grupo = ancho_txt(max((len(t) for _, _, _, _, t, _ in items), default=8), self.fs_grupo)
        self.R_tot = (self.r_gtxt + self.w_grupo if self.usar_grupos else r) + 0.15

        self.W = max(2 * self.R_tot, 10.0)
        self.alto = e.alto or (2 * self.R_tot + self.TITULO_H + self.LEYENDA_H)
        self.bottom_frac = self.LEYENDA_H / self.alto

    def _texto_hoja(self, i):
        sp_ = self.c.tax[i].get("s")
        txt = sin_prefijo(sp_) if (self.e.texto_hoja == "especie" and sp_) else self.c.nombres[i]
        return txt[:48]

    # ------------------------------------------------------------------ polar
    def _ang(self, y):
        """Angle (mathematical radians) of row y: clockwise from the top, with the gap centred at the top."""
        gap = 2 * math.pi - self.arco_rad
        return math.pi / 2 - (gap / 2 + y / self.H * self.arco_rad)

    def _pt(self, x, y):
        th = self._ang(y)
        return x * math.cos(th), x * math.sin(th)

    def _arco(self, r, y1, y2):
        n = max(2, int(abs(y2 - y1) / self.H * self.arco_rad * 57.3) + 2)
        ys = np.linspace(y1, y2, n)
        th = np.array([self._ang(v) for v in ys])
        return list(zip(r * np.cos(th), r * np.sin(th)))

    def _texto_radial(self, r, y, txt, fs, color, zorder=4):
        th = self._ang(y)
        deg = math.degrees(th)
        izq = math.cos(th) < 0
        x, yy = r * math.cos(th), r * math.sin(th)
        self.ax.text(x, yy, txt, fontsize=fs, color=color, ha="right" if izq else "left", va="center",
                     rotation=deg + 180 if izq else deg, rotation_mode="anchor", zorder=zorder)

    def _anillo(self, r0, r1, y0, y1, color, zorder=3):
        a1, a2 = math.degrees(self._ang(y1)), math.degrees(self._ang(y0))
        self.ax.add_patch(Wedge((0, 0), r1, a1, a2, width=r1 - r0, facecolor=color, edgecolor="none", zorder=zorder))

    # ------------------------------------------------------------------ drawing
    def dibujar(self):
        log.info(f"[-] Drawing circular tree {self.c.titulo or NOMBRES_PANEL[self.panel]}...")
        self.fig = plt.figure(figsize=(self.W, self.alto), facecolor="white")
        self.ax = self.fig.add_axes([0.0, self.bottom_frac, 1.0, 1.0 - self.bottom_frac])
        self.ax.set_xlim(-self.W / 2, self.W / 2)
        alto_ax = self.alto - self.LEYENDA_H
        self.ax.set_ylim(-(alto_ax - self.TITULO_H) / 2, (alto_ax + self.TITULO_H) / 2)
        self.ax.set_aspect("equal")
        self.ax.axis("off")
        raiz, col = self._arbol_solo()
        self._dibujar_arbol_circular(raiz, col)
        self._titulos()
        self._leyenda()
        return self.fig

    def _dibujar_arbol_circular(self, raiz, colapsados):
        c, e, es, ax, lay = self.c, self.e, self.es, self.ax, self.lay1
        segs, cols = [], []
        for x1, y1, x2, y2, h in lay["segmentos"]:
            if abs(y1 - y2) < 1e-12:                      # horizontal segment -> radial
                segs.append([self._pt(x1, y1), self._pt(x2, y2)])
            else:                                         # vertical segment -> arc
                segs.append(self._arco(x1, y1, y2))
            cols.append(self._color_rama(h))
        ax.add_collection(LineCollection(segs, colors=cols, linewidths=e.grosor_rama, capstyle="round", zorder=2))

        s0, s1 = self.r_strip
        # with real branch lengths tips end at different radii: dotted guides out to the ring
        guias = [[self._pt(u["x_apice"], u["y"]), self._pt(s0, u["y"])]
                 for u in lay["unidades"] if not u["tri"] and u["x_apice"] < s0 - 0.1]
        if guias:
            ax.add_collection(LineCollection(guias, colors=GRIS_ATENUADO, linewidths=0.4, linestyles=":", zorder=1))
        for u in lay["unidades"]:
            n = u["nodo"]
            y, alto = u["y"], u["h"]
            if u["tri"]:
                cat = es.categoria_pura(n.hojas)
                color = es.col_cat[cat] if cat is not None else GRIS_RAMA
                if c.hay_res and not c.hl[n.hojas].any():
                    color = GRIS_ATENUADO
                hh = alto * 0.46
                pts = [self._pt(u["x_apice"], y)] + self._arco(u["x_base"], y - hh, y + hh)
                ax.add_patch(Polygon(pts, closed=True, facecolor=color, edgecolor=color, alpha=0.9, lw=0.6, zorder=3))
                fc = color
            else:
                i = int(n.hojas[0])
                fc = self._color_hoja(i) if (not c.hay_res or c.hl[i]) else GRIS_ATENUADO
                alto = 1.0
            self._anillo(s0, s1, y - alto / 2 + 0.06, y + alto / 2 - 0.06, fc)
            if e.etiquetas_hojas and not u["tri"]:
                i = int(n.hojas[0])
                self._texto_radial(self.r_hoja, y, self._texto_hoja(i), self.fs_hoja,
                                   TINTA if (not c.hay_res or c.hl[i]) else GRIS_ATENUADO)
        if self.usar_grupos:
            self._rotular_grupos_circular(lay, raiz, colapsados)
        if e.soporte != "ninguno":
            self._dibujar_soporte(lay, raiz, colapsados, +1)
        if e.mostrar_ids:
            for n in lay["y"]:
                if (not n.es_hoja) and n.id and n.n >= e.min_id and n not in colapsados:
                    x, yy = self._pt(lay["x"][n], lay["y"][n])
                    ax.text(x, yy, n.id, fontsize=max(3.0, self.pitch_pt * 0.5), color=GRIS_MIXTO,
                            ha="center", va="bottom", zorder=4)

    def _rotular_grupos_circular(self, lay, raiz, colapsados):
        items = self._items_grupos(lay, raiz, colapsados)
        g0, g1 = self.r_gbar
        paso_txt = self.r_gtxt * self.arco_rad / self.H * 72.0          # pt per row at the text radius
        alto_txt = self.fs_grupo / max(paso_txt, 1.0) * 1.15
        colocados = []
        for _, _nn, ya, yb, txt, color in items:
            self._anillo(g0, g1, ya + 0.1, max(yb - 0.1, ya + 0.2), color)
            yc = (ya + yb) / 2
            if all(abs(yc - yp) >= alto_txt for yp in colocados):
                colocados.append(yc)
                self._texto_radial(self.r_gtxt, yc, txt, self.fs_grupo, TINTA)

    def _titulos(self):
        c, ax = self.c, self.ax
        y_top = ax.get_ylim()[1]
        ax.text(0, y_top - 0.35, c.titulo or NOMBRES_PANEL[self.panel], fontsize=22, ha="center", va="center",
                weight="bold", color=TINTA)
        sub = c.subtitulo or f"{c.N} genomes"
        if (c.colV if self.panel == "vamphyre" else c.colG):
            sub += f" · collapsed clades: {NOMBRES_COLAPSO.get(c.modo_colapso, c.modo_colapso)}"
        ax.text(0, y_top - 0.8, sub, fontsize=11, ha="center", va="center", color=TINTA_SUAVE)
