"""GTDB taxonomy (d__;p__;...), taxon selection and colour scheme."""
import collections
import re
from dataclasses import dataclass, field

import numpy as np

import colorsys

from .constants import (GRIS_OTROS, GRIS_SIN, LETRA_A_NOMBRE, MAX_DISTINGUIBLE, MAX_SEGURO_DALTONISMO,
                        NOMBRE_RANGO, PALETA, PALETA_EXTRA, RANGOS, SIN)
from .errors import ErrorDatos
import logging

log = logging.getLogger(__name__)


def parsear_tax(cadena):
    d = {}
    for tok in str(cadena).split(";"):
        tok = tok.strip()
        if len(tok) > 3 and tok[1:3] == "__":
            d[tok[0]] = tok
    return d


def sin_prefijo(tok):
    return tok[3:] if tok and len(tok) > 3 and tok[1:3] == "__" else tok


def tokens_de_hoja(tax_d):
    s = set()
    for t in tax_d.values():
        s.add(t.lower())
        s.add(sin_prefijo(t).lower())
    return s


def etiqueta_consenso(idxs, tax, umbral=0.8):
    n = len(idxs)
    for r in reversed(RANGOS):
        c = collections.Counter(tax[i].get(r) for i in idxs)
        t, k = c.most_common(1)[0]
        if t is not None and k / n >= umbral:
            return t, k / n
    return "mixed", 0.0


def genero_base(tok):
    return re.sub(r"_[A-Z]{1,2}$", "", sin_prefijo(tok)) if tok else ""


def rango_letra(txt):
    """'family' | 'familia' | 'f' -> 'f'. Raises ErrorDatos if the rank does not exist."""
    t = str(txt).strip().lower()
    if t in NOMBRE_RANGO:
        return NOMBRE_RANGO[t]
    if t in RANGOS:
        return t
    raise ErrorDatos(f"Unknown rank: {txt}")


def conj_taxa(toks, txt):
    """Indices of leaves that belong to any of the taxa in `txt` ('a,b,c').

    Each taxon may be written with prefix (g__Halorubrum) or without it (Halorubrum);
    case-insensitive. Returns None if `txt` is empty."""
    if not txt:
        return None
    qs = {q.strip().lower() for q in txt.split(",") if q.strip()}
    return {i for i, t in toks.items() if t & qs}


@dataclass
class EsquemaColor:
    """Colour assignment per taxon for one rank. Computed with ALL genomes, so
    colours do not change when filtering with --solo-taxa."""
    rango: str                  # rank letter (e.g. 'f')
    hay_tax: bool               # False if there is no taxonomy table
    cat: list                   # category of each leaf (or SIN)
    arr: np.ndarray             # numeric category id per leaf
    col_cat: dict               # category -> hex colour
    orden: list                 # categories by descending frequency
    conteo: collections.Counter = field(default_factory=collections.Counter)

    def categoria_pura(self, hojas):
        """Category shared by all `hojas`, or None if mixed/unclassified."""
        v = self.arr[hojas]
        if bool((v == v[0]).all()) and self.cat[int(hojas[0])] != SIN:
            return self.cat[int(hojas[0])]
        return None


def paleta(n):
    """List of `n` colours: the 8 validated ones, then the 12 extra ones and, if more are
    needed, generated tones (golden-angle hue steps in HLS, alternating lightness and
    saturation so neighbours in the list do not look alike). The first 20 never change."""
    base = PALETA + PALETA_EXTRA
    if n <= len(base):
        return base[:n]
    out = list(base)
    k = 0
    while len(out) < n:
        h = (0.13 + k * 0.618033988749895) % 1.0
        luz = (0.38, 0.55, 0.30, 0.47)[k % 4]
        sat = (0.75, 0.60, 0.85)[k % 3]
        r, g, b = colorsys.hls_to_rgb(h, luz, sat)
        col = "#{:02x}{:02x}{:02x}".format(round(r * 255), round(g * 255), round(b * 255))
        if col not in out and col not in (GRIS_OTROS, GRIS_SIN):
            out.append(col)
        k += 1
    return out


def aviso_colores(top):
    """Warning message for `top` colours, or None if the palette is safe."""
    if top > MAX_DISTINGUIBLE:
        return (f"{top} colours: above {MAX_DISTINGUIBLE}, colours are hard to tell apart even with normal "
                f"vision and are NOT suitable for readers with colour blindness. For publication figures "
                f"use ≤ {MAX_SEGURO_DALTONISMO} (or combine with highlighting, collapsing or a higher rank).")
    if top > MAX_SEGURO_DALTONISMO:
        return (f"{top} colours: colours {MAX_SEGURO_DALTONISMO + 1}+ are not validated for colour blindness. "
                f"For accessible figures use ≤ {MAX_SEGURO_DALTONISMO}, or combine with highlighting.")
    return None


def construir_esquema(tax, rc, top):
    """Build the EsquemaColor: the `top` most frequent taxa get their own colour
    (fixed palette order); the rest are grey ('other')."""
    hay_tax = any(tax)
    cat = [(t.get(rc) or SIN) for t in tax] if hay_tax else [SIN] * len(tax)
    aviso = aviso_colores(top)
    if aviso:
        log.warning(f"    [!] {aviso}")
    cont = collections.Counter(c for c in cat if c != SIN)
    orden = [c for c, _ in sorted(cont.items(), key=lambda kv: (-kv[1], kv[0]))]
    colores = paleta(min(top, len(orden)))
    col_cat = {c: (colores[i] if i < len(colores) else GRIS_OTROS)
               for i, c in enumerate(orden)}
    col_cat[SIN] = GRIS_SIN
    ids = {c: k for k, c in enumerate(sorted(set(cat)))}
    arr = np.array([ids[c] for c in cat])
    return EsquemaColor(rango=rc, hay_tax=hay_tax, cat=cat, arr=arr, col_cat=col_cat,
                        orden=orden, conteo=cont)


def nombre_rango(letra):
    return LETRA_A_NOMBRE[letra]
