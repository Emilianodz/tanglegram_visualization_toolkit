"""Collapse clades into triangles by concordance, taxonomy or manual selection."""
import logging
from dataclasses import dataclass

import numpy as np

from .constants import SIN
from .errors import ErrorDatos
from .taxonomy import conj_taxa
from .tree import preorden

log = logging.getLogger(__name__)


@dataclass
class OpcionesColapso:
    modo: str = "ninguno"            # ninguno (none) | concordantes | discordantes | taxonomia
    umbral: float = None             # clade compatibility (see seleccionar_colapsos)
    min_colapso: int = 3
    max_colapso: int = None
    rango: str = "f"                 # rank letter for 'taxonomia' mode
    colapsar_nodos: str = None       # "V12,G40"
    expandir_nodos: str = None
    colapsar_taxa: str = None
    expandir_taxa: str = None


def seleccionar_colapsos(raiz, info, modo, umbral, min_colapso, max_colapso, pura, bloqueado, forzados):
    """Return the set of nodes to collapse (the outermost ones that qualify).

    concordantes: subtree compatibility >= umbral (default 1.0 = identical topology)
    discordantes: compatibility < umbral (default 0.5)
    taxonomia   : clades made of a single taxon (`pura(node)`)
    `forzados` nodes are always collapsed; `bloqueado(node)` nodes never are."""
    if umbral is None:
        umbral = 1.0 if modo == "concordantes" else 0.5

    def cumple(n):
        if n in forzados:
            return True
        if modo == "concordantes":
            return info[n]["comp"] >= umbral - 1e-12
        if modo == "discordantes":
            return info[n]["comp"] < umbral
        if modo == "taxonomia":
            return pura(n)
        return False

    sel = set()
    if modo == "ninguno" and not forzados:
        return sel
    pila = [raiz]
    while pila:
        n = pila.pop()
        if n.es_hoja:
            continue
        ok_tam = n.n >= min_colapso and (max_colapso is None or n.n <= max_colapso)
        if n.padre is not None and not bloqueado(n) and ((n in forzados) or (ok_tam and cumple(n))):
            sel.add(n)
            continue
        pila.extend(n.hijos)
    return sel


def _pureza_por_rango(tax, rango):
    """Return pura(node): do all leaves share one taxon at `rango` (other than SIN)?"""
    catcol = [(t.get(rango) or SIN) for t in tax]
    ids_cat = {c: k for k, c in enumerate(sorted(set(catcol)))}
    arr = np.array([ids_cat[c] for c in catcol])
    sin_id = ids_cat.get(SIN, -1)

    def pura(n):
        v = arr[n.hojas]
        return n.n >= 2 and v[0] != sin_id and bool((v == v[0]).all())

    return pura


def planificar_colapsos(rV, rG, infoV, infoG, idsV, idsG, tax, toks, todo, op):
    """Decide which nodes of each tree are collapsed. Returns (colV, colG)."""
    pura = _pureza_por_rango(tax, op.rango)

    # leaves that must stay visible (they block collapsing any ancestor)
    Mexp = 0
    for i in (conj_taxa(toks, op.expandir_taxa) or ()):
        Mexp |= 1 << i
    forz = {"V": set(), "G": set()}
    for lista, accion in ((op.colapsar_nodos, "col"), (op.expandir_nodos, "exp")):
        for tok in [t.strip().upper() for t in (lista or "").split(",") if t.strip()]:
            nodo = (idsV if tok.startswith("V") else idsG).get(tok)
            if nodo is None:
                log.warning(f"    [!] Unknown node ID: {tok}")
            elif accion == "col":
                forz[tok[0]].add(nodo)
            else:
                Mexp |= nodo.mask
    Mfz = 0
    for i in (conj_taxa(toks, op.colapsar_taxa) or ()):
        Mfz |= 1 << i
    if Mfz:
        for r, k in ((rV, "V"), (rG, "G")):
            for n in preorden(r):
                if not n.es_hoja and n.padre is not None and (n.mask & ~Mfz) & todo == 0 and n.n >= 2:
                    forz[k].add(n)

    def bloqueado(n):
        return bool(n.mask & Mexp)

    colV = seleccionar_colapsos(rV, infoV, op.modo, op.umbral, op.min_colapso, op.max_colapso, pura, bloqueado, forz["V"])
    colG = seleccionar_colapsos(rG, infoG, op.modo, op.umbral, op.min_colapso, op.max_colapso, pura, bloqueado, forz["G"])
    if colV or colG:
        log.info(f"[+] Collapsed clades: VAMPhyRE {len(colV)} | GTDB {len(colG)}")
    return colV, colG


def planificar_colapsos_arbol(raiz, ids, tax, toks, todo, op):
    """Collapse for a single tree: by taxonomy or by chosen nodes/taxa
    (the 'concordantes'/'discordantes' modes compare two trees and do not apply here).
    `ids`: {'N12': node}. Returns the set of collapsed nodes."""
    if op.modo in ("concordantes", "discordantes"):
        raise ErrorDatos(f"--colapsar {op.modo} compares two trees: for a single tree use 'taxonomia' or 'ninguno'.")
    pura = _pureza_por_rango(tax, op.rango)
    Mexp, forz = 0, set()
    for i in (conj_taxa(toks, op.expandir_taxa) or ()):
        Mexp |= 1 << i
    for lista, accion in ((op.colapsar_nodos, "col"), (op.expandir_nodos, "exp")):
        for tok in [t.strip().upper() for t in (lista or "").split(",") if t.strip()]:
            nodo = ids.get(tok)
            if nodo is None:
                log.warning(f"    [!] Unknown node ID: {tok}")
            elif accion == "col":
                forz.add(nodo)
            else:
                Mexp |= nodo.mask
    Mfz = 0
    for i in (conj_taxa(toks, op.colapsar_taxa) or ()):
        Mfz |= 1 << i
    if Mfz:
        for n in preorden(raiz):
            if not n.es_hoja and n.padre is not None and (n.mask & ~Mfz) & todo == 0 and n.n >= 2:
                forz.add(n)
    sel = seleccionar_colapsos(raiz, None, op.modo, op.umbral, op.min_colapso, op.max_colapso, pura,
                               lambda n: bool(n.mask & Mexp), forz)
    if sel:
        log.info(f"[+] Collapsed clades: {len(sel)}")
    return sel
