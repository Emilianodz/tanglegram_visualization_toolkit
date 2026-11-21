"""Bipartitions (splits), Robinson-Foulds distance and discordance analysis.

Each leaf is one bit; a split is stored as the integer of its canonical side
(the side that does NOT contain the lowest-index leaf), so it is independent of
the root and comparable between trees."""
import logging

from .bitutils import bits, pc
from .tree import hojas_de, postorden, preorden

log = logging.getLogger(__name__)


def canon(m, todo):
    m &= todo
    c = pc(m)
    if c < 2 or pc(todo) - c < 2:
        return None
    return todo ^ m if (m & (todo & -todo)) else m


def conj_splits(raiz, todo):
    s = set()
    for n in preorden(raiz):
        if n.padre is None or n.es_hoja:
            continue
        c = canon(n.mask, todo)
        if c is not None:
            s.add(c)
    return s


def _restringir(lista, todo):
    lb = todo & -todo
    t = pc(todo)
    out = set()
    for m in lista:
        m &= todo
        c = pc(m)
        if c < 2 or t - c < 2:
            continue
        out.add(todo ^ m if m & lb else m)
    return out


def rebeldes_voraz(LA, LB, todo, max_pasos):
    """Greedily remove, step by step, the genome whose removal lowers RF the most."""
    pasos, actual = [], todo
    rf_act = len(_restringir(LA, actual) ^ _restringir(LB, actual))
    loo = {}
    for paso in range(max_pasos):
        if rf_act == 0:
            break
        mejor_i, mejor_rf = None, None
        for i in bits(actual):
            t2 = actual & ~(1 << i)
            rf2 = len(_restringir(LA, t2) ^ _restringir(LB, t2))
            if paso == 0:
                loo[i] = rf_act - rf2
            if mejor_rf is None or rf2 < mejor_rf:
                mejor_i, mejor_rf = i, rf2
        if mejor_rf >= rf_act:
            break
        pasos.append((mejor_i, rf_act, mejor_rf))
        actual &= ~(1 << mejor_i)
        rf_act = mejor_rf
    return pasos, loo


def asignar_ids(raiz, pref):
    k, pila = 0, [raiz]
    por_id = {}
    while pila:
        n = pila.pop()
        if n.es_hoja:
            continue
        if n.padre is not None:
            k += 1
            n.id = f"{pref}{k}"
            por_id[n.id] = n
        hs = sorted(n.hijos, key=lambda h: (-h.n, int(h.hojas.min())))
        pila.extend(reversed(hs))
    return por_id


def fragmentos(raiz, S):
    """Minimum number of disjoint clades that exactly cover group S (1 = monophyletic)."""
    k, pila = 0, [raiz]
    while pila:
        n = pila.pop()
        if n.mask & S == 0:
            continue
        if n.mask & ~S == 0:
            k += 1
            continue
        pila.extend(n.hijos)
    return k


def info_arbol(raiz, otra_raiz, otra_splits, todo):
    """Per node: concordance (does its split exist in the other tree?), subtree
    compatibility, and best matching clade in the other tree (Jaccard)."""
    clados = [(n.mask, n) for n in preorden(otra_raiz) if n.padre is not None and not n.es_hoja]
    info = {}
    for n in postorden(raiz):
        d = {"conc": None, "cc": 0, "ct": 0, "comp": 1.0, "jmax": None, "pareja": None, "dif": 0}
        if not n.es_hoja:
            cc = sum(info[h]["cc"] for h in n.hijos)
            ct = sum(info[h]["ct"] for h in n.hijos)
            if n.padre is not None:
                c = canon(n.mask, todo)
                conc = True if c is None else (c in otra_splits)
                d["conc"] = conc
                cc += int(conc)
                ct += 1
            d["cc"], d["ct"] = cc, ct
            d["comp"] = cc / ct if ct else 1.0
        info[n] = d
    for n, d in info.items():
        if n.es_hoja or n.padre is None:
            continue
        if d["conc"]:
            d["jmax"], d["dif"] = 1.0, 0
            continue
        mejor_j, mejor_o, mejor_i = -1.0, None, 0
        for m, o in clados:
            i = pc(n.mask & m)
            if i == 0:
                continue
            j = i / (n.n + o.n - i)
            if j > mejor_j:
                mejor_j, mejor_o, mejor_i = j, o, i
        d["jmax"], d["pareja"] = mejor_j, mejor_o
        d["dif"] = n.n + (mejor_o.n if mejor_o else 0) - 2 * mejor_i
    return info


def ancestros_concordantes(raiz, info):
    r = {}
    for h in hojas_de(raiz):
        tot = c = 0
        x = h.padre
        while x is not None and x.padre is not None:
            tot += 1
            c += 1 if info[x]["conc"] else 0
            x = x.padre
        r[int(h.hojas[0])] = (c, tot)
    return r
