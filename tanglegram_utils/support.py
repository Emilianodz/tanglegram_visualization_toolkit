"""Node support values (bootstrap, SH-like, posterior probabilities...).

GTDB-Tk trees carry 0-1 support; bootstrap usually ranges 0-100. Everything is
normalised to 0-100 here, and poorly supported nodes can be collapsed into
polytomies, as is usually done before presenting a consensus tree."""
import statistics

from .tree import postorden, preorden


def _internos_con_soporte(raiz):
    return [n for n in preorden(raiz) if not n.es_hoja and n.padre is not None and n.soporte is not None]


def escala_soporte(raiz):
    """Factor that brings support to 0-100: 100 if everything is within 0-1, otherwise 1."""
    vals = [n.soporte for n in _internos_con_soporte(raiz)]
    return 100.0 if vals and max(vals) <= 1.0 + 1e-9 else 1.0


def normalizar_soporte(raiz):
    """Rescale every node's support to 0-100 (in place). Returns the factor applied."""
    f = escala_soporte(raiz)
    if f != 1.0:
        for n in _internos_con_soporte(raiz):
            n.soporte *= f
    return f


def resumen_soporte(raiz):
    """Support statistics (assumes 0-100 scale); None if the tree has no support values."""
    internos = [n for n in preorden(raiz) if not n.es_hoja and n.padre is not None]
    vals = [n.soporte for n in internos if n.soporte is not None]
    if not vals:
        return None
    n = len(vals)
    return {"nodos_internos": len(internos), "con_soporte": n, "minimo": min(vals),
            "mediana": statistics.median(vals),
            "pct_ge_70": 100.0 * sum(v >= 70 for v in vals) / n,
            "pct_ge_90": 100.0 * sum(v >= 90 for v in vals) / n}


def contraer_soporte_bajo(raiz, umbral):
    """Turn every internal node with support < umbral (0-100 scale) into a polytomy.

    The node's children move to its parent and the node's branch is removed (like ape's
    `di2multi`): child branches keep their length, so distances between leaves that
    crossed that branch become shorter by its length.
    Returns (root, number of collapsed nodes); nodes without support are left untouched."""
    k = 0
    for n in postorden(raiz):
        if n.es_hoja or n.padre is None or n.soporte is None or n.soporte >= umbral:
            continue
        p = n.padre
        i = p.hijos.index(n)
        for h in n.hijos:
            h.padre = p
        p.hijos[i:i + 1] = n.hijos
        k += 1
    return raiz, k
