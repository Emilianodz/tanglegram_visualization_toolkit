"""Branch ordering: ladderize and untangle (minimise line crossings)."""
import numpy as np

from .tree import postorden, preorden


def ladderize(raiz):
    for n in preorden(raiz):
        if n.hijos:
            n.hijos.sort(key=lambda h: (h.n, int(h.hojas.min())))


def posiciones(raiz, N):
    pos = np.full(N, -1, dtype=np.int64)
    k = 0
    for n in preorden(raiz):
        if n.es_hoja:
            pos[n.hojas[0]] = k
            k += 1
    return pos


def contar_cruces(rA, rB, N):
    pa, pb = posiciones(rA, N), posiciones(rB, N)
    ok = pa >= 0
    orden = np.argsort(pa[ok])
    seq = pb[ok][orden]
    return int(np.triu(seq[:, None] > seq[None, :], k=1).sum())


def orientar(raiz, pos_otro):
    """Optimal orientation of each binary node with respect to the other tree's leaf order."""
    cambios = 0
    for n in postorden(raiz):
        if n.es_hoja or len(n.hijos) < 2:
            continue
        if len(n.hijos) == 2:
            a, b = n.hijos
            A = pos_otro[a.hojas]
            B = np.sort(pos_otro[b.hojas])
            c_ab = int(np.searchsorted(B, A).sum())
            if c_ab > len(A) * len(B) - c_ab:
                n.hijos = [b, a]
                cambios += 1
        else:
            nuevo = sorted(n.hijos, key=lambda h: float(pos_otro[h.hojas].mean()))
            if nuevo != n.hijos:
                n.hijos = nuevo
                cambios += 1
    return cambios


def desenredar(rI, rD, N, iteraciones=12):
    """Reorder the branches of both trees to minimise line crossings.

    Returns (crossings_before, crossings_after). "Before" is measured after a ladderize
    (children sorted by size), not on the raw Newick order."""
    ladderize(rI)
    ladderize(rD)
    antes = contar_cruces(rI, rD, N)
    for _ in range(iteraciones):
        c = orientar(rD, posiciones(rI, N))
        c += orientar(rI, posiciones(rD, N))
        if c == 0:
            break
    return antes, contar_cruces(rI, rD, N)
