import pytest

from conftest import NWK_G, NWK_V, preparar
from tanglegram_utils.bitutils import pc
from tanglegram_utils.rooting import reenraizar_en
from tanglegram_utils.splits import (_restringir, canon, conj_splits, fragmentos,
                                     rebeldes_voraz)
from tanglegram_utils.tree import indexar, podar, preorden


def rf(a, b, todo):
    return len(conj_splits(a, todo) ^ conj_splits(b, todo))


def test_rf_cero_si_misma_topologia_rotada():
    rA, rB, _, todo = preparar("((a,b),(c,d),e);", "(e,(d,c),(b,a));")
    assert rf(rA, rB, todo) == 0


def test_rf_conocido():
    # {ab|cde, cd|abe} vs {ac|bde, bd|ace}: none in common -> RF = 4
    rA, rB, _, todo = preparar("((a,b),(c,d),e);", "((a,c),(b,d),e);")
    assert rf(rA, rB, todo) == 4


def test_rf_coincide_con_dendropy(arboles):
    dendropy = pytest.importorskip("dendropy")
    from dendropy.calculate import treecompare
    ns = dendropy.TaxonNamespace()
    kw = dict(schema="newick", taxon_namespace=ns, rooting="force-unrooted")
    a, b = dendropy.Tree.get(data=NWK_V, **kw), dendropy.Tree.get(data=NWK_G, **kw)
    a.encode_bipartitions()
    b.encode_bipartitions()
    rA, rB, _, todo = arboles
    assert rf(rA, rB, todo) == treecompare.symmetric_difference(a, b) > 0


def test_splits_no_dependen_de_la_raiz(arboles):
    rA, _, idx, todo = arboles
    esperado = conj_splits(rA, todo)
    for nodo in [n for n in preorden(rA) if n.padre is not None]:
        nueva = reenraizar_en(rA, nodo)
        indexar(nueva, idx)
        assert conj_splits(nueva, todo) == esperado
        rA = nueva


def test_podar_equivale_a_restringir_splits(arboles):
    rA, _, idx, todo = arboles
    quitar = {"a1", "c3"}
    keep_mask = todo
    for q in quitar:
        keep_mask &= ~(1 << idx[q])
    esperado = _restringir(list(conj_splits(rA, todo)), keep_mask)
    podado = podar(rA, set(idx) - quitar)
    indexar(podado, idx)
    assert conj_splits(podado, keep_mask) == esperado


def test_canon_descarta_triviales():
    todo = 0b1111
    assert canon(0b0001, todo) is None          # a single leaf
    assert canon(0b0111, todo) is None          # the complement is a single leaf
    assert canon(0b0011, todo) == canon(0b1100, todo)


def test_rebeldes_voraz_encuentra_las_hojas_movidas():
    # f and g swapped: removing both leaves RF = 0
    rA, rB, _, todo = preparar("((a,b),(c,d),(e,f),g);", "((a,b),(c,d),(e,g),f);")
    LA, LB = list(conj_splits(rA, todo)), list(conj_splits(rB, todo))
    pasos, loo = rebeldes_voraz(LA, LB, todo, 5)
    assert pasos[0][1] > pasos[0][2]            # the first step lowers RF
    assert pasos[-1][2] == 0 and len(pasos) <= 2
    assert all(v >= 0 for v in loo.values())


def test_fragmentos(arboles):
    rA, rB, idx, _ = arboles
    S = sum(1 << idx[x] for x in ("a1", "a2", "a3", "a4"))
    assert fragmentos(rA, S) == 1               # monophyletic in V
    assert fragmentos(rB, S) > 1                # split in G (a3 moved away)
    assert pc(S) == 4
