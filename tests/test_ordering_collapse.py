from conftest import NWK_V, preparar
from tanglegram_utils.collapse import OpcionesColapso, planificar_colapsos
from tanglegram_utils.ordering import contar_cruces, desenredar
from tanglegram_utils.splits import asignar_ids, conj_splits, info_arbol

NWK_V_ROTADO = ("(((c3:1,c4:1):1,(c1:1,c2:1):1):1,((b3:1,b4:1):1,(b1:1,b2:1):1):1,"
                "((a3:1,a4:1):1,(a1:1,a2:1):1):1);")


def test_desenredar_deja_cero_cruces_si_mismo_arbol():
    rA, rB, idx, _ = preparar(NWK_V, NWK_V_ROTADO)
    assert contar_cruces(rA, rB, len(idx)) > 0          # raw file order: tangled
    _, despues = desenredar(rA, rB, len(idx))
    assert despues == 0


def test_desenredar_no_aumenta_cruces(arboles):
    rA, rB, idx, _ = arboles
    antes, despues = desenredar(rA, rB, len(idx))     # "before" = after ladderize
    assert despues <= antes
    assert despues == contar_cruces(rA, rB, len(idx))


def _info(arboles):
    rA, rB, idx, todo = arboles
    sA, sB = conj_splits(rA, todo), conj_splits(rB, todo)
    return (rA, rB, info_arbol(rA, rB, sB, todo), info_arbol(rB, rA, sA, todo),
            asignar_ids(rA, "V"), asignar_ids(rB, "G"), todo)


def test_colapsar_concordantes_solo_toma_subarboles_identicos(arboles):
    rA, rB, iA, iB, idsA, idsB, todo = _info(arboles)
    tax = [{} for _ in range(12)]
    op = OpcionesColapso(modo="concordantes", min_colapso=3)
    colA, colB = planificar_colapsos(rA, rB, iA, iB, idsA, idsB, tax, {i: set() for i in range(12)}, todo, op)
    # the c* clade is identical in both trees and must be collapsed in both
    assert any(n.n == 4 for n in colA) and any(n.n == 4 for n in colB)
    # no clade containing a3/b1 (the swapped leaves) is identical
    assert all(iA[n]["comp"] == 1.0 for n in colA) and all(iB[n]["comp"] == 1.0 for n in colB)


def test_modo_ninguno_no_colapsa(arboles):
    rA, rB, iA, iB, idsA, idsB, todo = _info(arboles)
    colA, colB = planificar_colapsos(rA, rB, iA, iB, idsA, idsB, [{}] * 12, {i: set() for i in range(12)},
                                     todo, OpcionesColapso())
    assert not colA and not colB


def test_colapsar_nodos_forzado_y_expandir(arboles):
    rA, rB, iA, iB, idsA, idsB, todo = _info(arboles)
    tox = {i: set() for i in range(12)}
    nodo = next(n for n in idsA.values() if n.n == 4)
    colA, _ = planificar_colapsos(rA, rB, iA, iB, idsA, idsB, [{}] * 12, tox, todo,
                                  OpcionesColapso(colapsar_nodos=nodo.id))
    assert nodo in colA
    colA, _ = planificar_colapsos(rA, rB, iA, iB, idsA, idsB, [{}] * 12, tox, todo,
                                  OpcionesColapso(colapsar_nodos=nodo.id, expandir_nodos=nodo.id))
    assert nodo not in colA
