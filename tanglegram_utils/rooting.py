"""Tree rooting (affects drawing only; RF does not depend on the root)."""
import logging

from .bitutils import pc
from .splits import canon, conj_splits
from .tree import Nodo, indexar, limpiar_unarios, preorden

log = logging.getLogger(__name__)


def reenraizar_en(raiz, v, frac=0.5):
    """Place a new root on the edge above v. Returns the new root.

    `frac` is the fraction of the edge left on v's side (0.5 = midway)."""
    if v is raiz:
        return raiz
    if v.padre is raiz and len(raiz.hijos) == 2:
        return raiz
    cadena, x = [], v.padre
    while x is not None:
        cadena.append(x)
        x = x.padre
    attrs = [(c.largo, c.soporte) for c in cadena]
    p0 = cadena[0]
    p0.hijos.remove(v)
    for i in range(len(cadena) - 1):
        a, b = cadena[i], cadena[i + 1]
        b.hijos.remove(a)
        a.hijos.append(b)
        b.padre = a
        b.largo, b.soporte = attrs[i]
    L, S = v.largo, v.soporte
    nueva = Nodo(None, 0.0, None)
    v.padre = nueva
    v.largo = L * frac
    p0.padre = nueva
    p0.largo = L * (1.0 - frac)
    p0.soporte = S
    nueva.hijos = [v, p0]
    return limpiar_unarios(nueva)


def raiz_comun(rV, rG, todo, idx, tokens_hoja, raiz_taxon):
    """Pick a split present in both trees to root them the same way (visual only)."""
    sV, sG = conj_splits(rV, todo), conj_splits(rG, todo)
    comunes = sV & sG
    N = pc(todo)
    if raiz_taxon:
        M = 0
        for i, toks in tokens_hoja.items():
            if raiz_taxon.lower() in toks:
                M |= 1 << i
        c = canon(M, todo) if M else None
        if c is not None and c in comunes:
            return c, "taxon " + raiz_taxon
        log.warning(f"    [!] --raiz-taxon {raiz_taxon}: not monophyletic in both trees; using the automatic root.")
    if not comunes:
        return None, "no common split"
    mejor = max(comunes, key=lambda s: (min(pc(s), N - pc(s)), -s))
    return mejor, "most balanced common split"


def enraizar_en_split(raiz, split, todo, idx):
    for n in preorden(raiz):
        if n.padre is None or n.es_hoja:
            continue
        if canon(n.mask, todo) == split:
            nueva = reenraizar_en(raiz, n)
            indexar(nueva, idx)
            return nueva
    return raiz


def enraizar_balanceado(raiz, todo, idx):
    """Fallback: root on the most balanced split of the tree itself."""
    N, mejor, mejor_n = pc(todo), None, -1
    for n in preorden(raiz):
        if n.padre is None or n.es_hoja:
            continue
        b = min(n.n, N - n.n)
        if b > mejor_n:
            mejor, mejor_n = n, b
    nueva = reenraizar_en(raiz, mejor)
    indexar(nueva, idx)
    return nueva


def _desenraizar(raiz):
    """If the root is bifurcating, merge it with its internal child (making it trifurcating),
    so the root edge becomes an ordinary edge that can be rerooted on."""
    if len(raiz.hijos) != 2:
        return raiz
    a, b = raiz.hijos
    interno, otro = (a, b) if not a.es_hoja else (b, a)
    if interno.es_hoja:
        return raiz
    otro.largo += interno.largo
    raiz.hijos = [otro] + interno.hijos
    for h in raiz.hijos:
        h.padre = raiz
    return raiz


def enraizar_punto_medio(raiz):
    """Root at the midpoint of the longest leaf-to-leaf path (using branch lengths).

    Without branch lengths (all 0) there is no midpoint: the tree is returned unchanged."""
    raiz = _desenraizar(raiz)
    nodos = preorden(raiz)
    vecinos = {n: [] for n in nodos}
    for n in nodos:
        if n.padre is not None:
            vecinos[n].append((n.padre, n.largo))
            vecinos[n.padre].append((n, n.largo))

    def distancias(origen):
        dist, previo, pila = {origen: 0.0}, {origen: None}, [origen]
        while pila:
            u = pila.pop()
            for w, l in vecinos[u]:
                if w not in dist:
                    dist[w], previo[w] = dist[u] + l, (u, l)
                    pila.append(w)
        return dist, previo

    hojas = [n for n in nodos if n.es_hoja]
    d0, _ = distancias(hojas[0])
    A = max(hojas, key=lambda h: d0[h])
    dA, previo = distancias(A)
    B = max(hojas, key=lambda h: dA[h])
    D = dA[B]
    if D <= 0:
        return raiz
    # path B -> A; the midpoint lies at D/2 from B
    x, acum, mitad = B, 0.0, D / 2.0
    while previo[x] is not None:
        u, l = previo[x]
        if acum + l >= mitad:
            t = mitad - acum                      # distance from x to the midpoint
            if x.padre is u:                      # x is the child end of the edge
                v, frac = x, (t / l if l else 0.5)
            else:                                 # u is the child end
                v, frac = u, ((l - t) / l if l else 0.5)
            return reenraizar_en(raiz, v, frac)
        acum += l
        x = u
    return raiz


def enraizar_en_grupo(raiz, mascara_grupo, todo, idx):
    """Root on the edge separating the outgroup from the rest, if it is a clade
    (or its complement is). Returns (root, ok); ok=False if it is not monophyletic."""
    for n in preorden(raiz):
        if n.padre is None:
            continue
        if n.mask == mascara_grupo or n.mask == (todo & ~mascara_grupo):
            nueva = reenraizar_en(raiz, n)
            indexar(nueva, idx)
            return nueva, True
    return raiz, False
