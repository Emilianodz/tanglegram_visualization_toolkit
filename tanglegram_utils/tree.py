"""Lightweight tree structure (Nodo) and basic operations: Newick reading,
traversals, pruning and bitmask indexing."""
import re

import numpy as np
from Bio import Phylo

from .names import limpiar_etiqueta


class Nodo:
    __slots__ = ("nombre", "largo", "soporte", "hijos", "padre", "mask", "n",
                 "hojas", "altura", "prof", "profmax", "id")

    def __init__(self, nombre=None, largo=0.0, soporte=None):
        self.nombre = nombre
        self.largo = largo
        self.soporte = soporte
        self.hijos = []
        self.padre = None
        self.mask = 0
        self.n = 0
        self.hojas = None
        self.altura = 0
        self.prof = 0.0
        self.profmax = 0.0
        self.id = None

    @property
    def es_hoja(self):
        return self.nombre is not None


_RE_SOPORTE = re.compile(r"^\s*([0-9]*\.?[0-9]+)\s*(?:[:/|,]|$)")


def soporte_de_clado(c):
    """Support (bootstrap, SH, probability...) of an internal node, or None.

    Taken from the node's numeric label ('95', '0.98', '1.00') or from the first
    number of composite labels ('95.0:g__Halorubrum', '87/100'). Non-numeric names
    (e.g. VAMPhyRE's 'Inner299') are not support."""
    if c.confidence is not None:
        return float(c.confidence)
    m = _RE_SOPORTE.match(c.name or "")
    return float(m.group(1)) if m else None


def desde_phylo(c):
    nodo = Nodo(None, float(c.branch_length or 0.0), soporte_de_clado(c) if c.clades else c.confidence)
    if c.clades:
        for h in c.clades:
            hn = desde_phylo(h)
            hn.padre = nodo
            nodo.hijos.append(hn)
    else:
        nodo.nombre = limpiar_etiqueta(c.name or "")
    return nodo


def cargar_arbol(ruta):
    return desde_phylo(Phylo.read(ruta, "newick").root)


def postorden(raiz):
    out, pila = [], [(raiz, False)]
    while pila:
        n, visto = pila.pop()
        if visto:
            out.append(n)
        else:
            pila.append((n, True))
            for h in reversed(n.hijos):
                pila.append((h, False))
    return out


def preorden(raiz):
    out, pila = [], [raiz]
    while pila:
        n = pila.pop()
        out.append(n)
        pila.extend(reversed(n.hijos))
    return out


def hojas_de(raiz):
    return [n for n in preorden(raiz) if n.es_hoja]


def limpiar_unarios(raiz):
    for n in postorden(raiz):
        if n.es_hoja or len(n.hijos) != 1:
            continue
        h = n.hijos[0]
        h.largo += n.largo
        if n.padre is None:
            h.padre = None
            raiz = h
        else:
            p = n.padre
            p.hijos[p.hijos.index(n)] = h
            h.padre = p
    return raiz


def podar(raiz, conservar):
    for n in postorden(raiz):
        if n.es_hoja:
            continue
        n.hijos = [h for h in n.hijos
                   if (h.nombre in conservar if h.es_hoja else bool(h.hijos))]
        for h in n.hijos:
            h.padre = n
    return limpiar_unarios(raiz)


def indexar(raiz, idx):
    """Assign bitmask, size, heights and depths."""
    for n in postorden(raiz):
        if n.es_hoja:
            n.mask = 1 << idx[n.nombre]
            n.n = 1
            n.altura = 0
            n.hojas = np.array([idx[n.nombre]], dtype=np.int64)
        else:
            n.mask, n.n, arrs = 0, 0, []
            for h in n.hijos:
                n.mask |= h.mask
                n.n += h.n
                arrs.append(h.hojas)
            n.altura = 1 + max(h.altura for h in n.hijos)
            n.hojas = np.concatenate(arrs)
    raiz.prof = 0.0
    for n in preorden(raiz):
        for h in n.hijos:
            h.prof = n.prof + h.largo
    for n in postorden(raiz):
        n.profmax = n.prof if n.es_hoja else max(h.profmax for h in n.hijos)


def desde_cadena(texto):
    """Build a tree from a Newick string (handy in tests and notebooks)."""
    from io import StringIO
    return desde_phylo(Phylo.read(StringIO(texto), "newick").root)


_RE_NOMBRE_SEGURO = re.compile(r"^[A-Za-z0-9_.\-+|]+$")


def _etiqueta_newick(texto):
    return texto if _RE_NOMBRE_SEGURO.match(texto) else "'" + texto.replace("'", "''") + "'"


def a_newick(raiz, soporte=True):
    """Tree as a Newick string: branch lengths and, if `soporte`, each internal node's support."""
    def rec(n):
        if n.es_hoja:
            txt = _etiqueta_newick(n.nombre)
        else:
            txt = "(" + ",".join(rec(h) for h in n.hijos) + ")"
            if soporte and n.soporte is not None and n.padre is not None:
                txt += f"{n.soporte:g}"
        if n.padre is not None:
            txt += f":{n.largo:.10g}"
        return txt
    return rec(raiz) + ";"
