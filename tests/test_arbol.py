import itertools
import os
import shlex

import pytest

from tanglegram_utils.cli_arbol import comando_equivalente, construir_parser, opciones_desde_args
from tanglegram_utils.collapse import OpcionesColapso
from tanglegram_utils.errors import ErrorDatos
from tanglegram_utils.models import EntradaArbol, OpcionesAnalisis, OpcionesArbol, Salida
from tanglegram_utils.plotting import EstiloFigura
from tanglegram_utils.rooting import enraizar_en_grupo, enraizar_punto_medio
from tanglegram_utils.servicio import ejecutar_arbol_en_directorio, ejecutar_en_directorio
from tanglegram_utils.support import contraer_soporte_bajo, normalizar_soporte, resumen_soporte
from tanglegram_utils.tree import a_newick, desde_cadena, hojas_de, indexar, preorden

# 12 leaves with support (0-1) on internal nodes
NWK_SOP = ("(((a1:1,a2:1)0.99:1,(a3:1,a4:1)0.5:1)0.95:1,((b1:1,b2:1)1.0:1,(b3:1,b4:1)0.4:1)0.9:1,"
           "((c1:1,c2:1)1.0:1,(c3:1,c4:1)0.97:1)0.88:1);")


def _distancias(raiz):
    def anc(n):
        out = []
        while n is not None:
            out.append(n)
            n = n.padre
        return out
    d = {}
    for a, b in itertools.combinations(hojas_de(raiz), 2):
        A, B = anc(a), anc(b)
        s = set(B)
        c = next(x for x in A if x in s)
        d[frozenset((a.nombre, b.nombre))] = round(
            sum(x.largo for x in A[:A.index(c)]) + sum(x.largo for x in B[:B.index(c)]), 9)
    return d


def _altura(n):
    return 0 if n.es_hoja else max(h.largo + _altura(h) for h in n.hijos)


# ---------------------------------------------------------------- support
def test_soporte_se_lee_de_etiquetas_numericas_y_compuestas():
    r = desde_cadena("(((a:1,b:1)90:1,c:1)'55.0:g__X':1,((d:1,e:1)Inner3:1,f:1)0.3:1,g:1);")
    vals = [n.soporte for n in preorden(r) if not n.es_hoja]
    assert vals == [None, 55.0, 90.0, 0.3, None]       # 'Inner3' (VAMPhyRE) is not support


def test_normalizar_y_contraer_soporte():
    r = desde_cadena(NWK_SOP)
    d0 = _distancias(r)
    assert normalizar_soporte(r) == 100.0
    res = resumen_soporte(r)
    assert res["minimo"] == pytest.approx(40.0) and res["con_soporte"] == 9
    r, k = contraer_soporte_bajo(r, 60)
    assert k == 2                                      # 0.5 and 0.4 -> polytomies
    d1 = _distancias(r)
    assert len(hojas_de(r)) == 12
    assert d1[frozenset(("a1", "a2"))] == d0[frozenset(("a1", "a2"))]            # does not cross the removed branch
    assert d1[frozenset(("a1", "a3"))] == d0[frozenset(("a1", "a3"))] - 1.0      # does: loses that branch (1.0)
    assert sorted(len(h.hijos) for h in r.hijos[0].hijos) == [0, 0, 2]     # (a1,a2) plus leaves a3, a4 hanging together from clade A


def test_arbol_sin_soporte_no_falla():
    r = desde_cadena("((a:1,b:1):1,(c:1,d:1):1);")
    assert resumen_soporte(r) is None and normalizar_soporte(r) == 1.0


# ---------------------------------------------------------------- root
@pytest.mark.parametrize("nwk", ["((a:1,b:2):0.5,(c:8,(d:1,e:1):1):2,f:0.2);",
                                 "((a:1,b:2):0.5,(c:8,(d:1,e:1):1):2);",
                                 "(a:1,b:1,(c:1,d:5):1);"])
def test_punto_medio_conserva_distancias_y_equilibra(nwk):
    r = desde_cadena(nwk)
    d0 = _distancias(r)
    r2 = enraizar_punto_medio(r)
    assert _distancias(r2) == d0
    alturas = [h.largo + _altura(h) for h in r2.hijos]
    assert alturas[0] == pytest.approx(alturas[1])


def test_enraizar_en_grupo_externo():
    r = desde_cadena("((a1:1,a2:1):1,(b1:1,b2:1):1,(c1:1,c2:1):1);")
    idx = {n: i for i, n in enumerate(sorted(h.nombre for h in hojas_de(r)))}
    indexar(r, idx)
    todo = (1 << len(idx)) - 1
    mask = (1 << idx["c1"]) | (1 << idx["c2"])
    r2, ok = enraizar_en_grupo(r, mask, todo, idx)
    assert ok
    # the outgroup ends up as one of the two sides of the root
    lados = [{hojas.nombre for hojas in hojas_de(h)} for h in r2.hijos]
    assert {"c1", "c2"} in lados
    # group that is not a clade -> no rerooting
    m2 = (1 << idx["a1"]) | (1 << idx["b1"])
    _, ok2 = enraizar_en_grupo(r2, m2, todo, idx)
    assert not ok2


def test_newick_ida_y_vuelta():
    r = desde_cadena(NWK_SOP)
    normalizar_soporte(r)
    txt = a_newick(r)
    r2 = desde_cadena(txt)
    assert a_newick(r2) == txt
    assert _distancias(r2) == _distancias(r)


# ---------------------------------------------------------------- pipeline and CLI
@pytest.fixture
def arbol_con_tabla(datos_sinteticos):
    d = datos_sinteticos
    (d / "sop.newick").write_text(NWK_SOP + "\n")
    return d


def _entrada(d, **kw):
    return EntradaArbol(arbol=str(d / "sop.newick"), tabla=str(d / "tabla.csv"), **kw)


@pytest.mark.parametrize("forma", ["circular", "rectangular"])
def test_arbol_individual_completo(arbol_con_tabla, tmp_path, forma):
    op = OpcionesArbol(color_rango="familia", raiz="punto-medio", contraer_soporte=60.0, rango_grupos="familia,genero")
    est = EstiloFigura(forma=forma, dpi=40, soporte="punto", soporte_min=80.0)
    res = ejecutar_arbol_en_directorio(_entrada(arbol_con_tabla), op, OpcionesColapso(), est, str(tmp_path / "o"),
                                       formato_extra="svg")
    assert res.n_genomas == 12 and res.n_con_taxonomia == 12
    assert res.n_contraidos == 2 and res.motivo_raiz == "midpoint"
    assert res.soporte["con_soporte"] == 9
    for suf in (".png", ".svg", "_arbol.nwk", "_nodos.csv", "_grupos.csv", "_reporte.md"):
        assert res.con_sufijo(suf), suf
    assert os.path.getsize(res.imagen) > 2000
    assert f"--forma {forma}" in res.comando or forma == "circular"        # circular is this CLI's default
    assert "--arbol sop.newick" in res.comando


def test_arbol_individual_color_colapso_zoom_y_resaltado(arbol_con_tabla, tmp_path):
    op = OpcionesArbol(color_rango="genero", resaltar="f__Fa", solo_taxa="f__Fa,f__Fb")
    col = OpcionesColapso(modo="taxonomia", rango="g", min_colapso=2)
    res = ejecutar_arbol_en_directorio(_entrada(arbol_con_tabla), op, col, EstiloFigura(forma="circular", dpi=40),
                                       str(tmp_path / "o"))
    assert res.n_genomas == 8
    assert "Collapsed clades" in res.log


def test_grupo_externo_y_orden(arbol_con_tabla, tmp_path):
    op = OpcionesArbol(raiz="taxon", raiz_taxon="f__Fc", orden="original")
    res = ejecutar_arbol_en_directorio(_entrada(arbol_con_tabla), op, OpcionesColapso(), EstiloFigura(dpi=40),
                                       str(tmp_path / "o"))
    assert res.motivo_raiz == "outgroup f__Fc"


def test_errores_de_datos_claros(arbol_con_tabla, tmp_path):
    est = EstiloFigura(dpi=40)
    for op, col, esperado in (
            (OpcionesArbol(raiz="taxon"), OpcionesColapso(), "--raiz taxon"),
            (OpcionesArbol(raiz="otra"), OpcionesColapso(), "Unknown root option"),
            (OpcionesArbol(), OpcionesColapso(modo="concordantes"), "compares two trees"),
            (OpcionesArbol(solo_taxa="f__Zz"), OpcionesColapso(), "fewer than 4")):
        with pytest.raises(ErrorDatos, match=esperado):
            ejecutar_arbol_en_directorio(_entrada(arbol_con_tabla), op, col, est, str(tmp_path / "x"))


def test_arbol_sin_tabla_usa_gris_y_no_escribe_grupos(arbol_con_tabla, tmp_path):
    ent = EntradaArbol(arbol=str(arbol_con_tabla / "sop.newick"))
    res = ejecutar_arbol_en_directorio(ent, OpcionesArbol(), OpcionesColapso(), EstiloFigura(dpi=40), str(tmp_path / "o"))
    assert res.n_con_taxonomia == 0 and not res.con_sufijo("_grupos.csv")


def test_comando_arbol_ida_y_vuelta():
    ent = EntradaArbol("a b.nwk", tabla="t.xlsx", mapa_ids="m.tsv", titulo="My tree")
    op = OpcionesArbol(color_rango="genero", raiz="punto-medio", contraer_soporte=70.0, orden="original")
    col = OpcionesColapso(modo="taxonomia", rango="f", min_colapso=5)
    est = EstiloFigura(forma="rectangular", arco=300.0, soporte="numero", soporte_min=90.0, etiquetas_hojas=True)
    sal = Salida(out_plot="salida/x.png", formato_extra="svg,pdf")
    cmd = comando_equivalente(ent, op, col, est, sal)
    args = construir_parser().parse_args(shlex.split(cmd.replace("\\\n", " "))[2:])
    assert opciones_desde_args(args) == (ent, op, col, est, sal)


# ---------------------------------------------------------------- single tree from the tanglegram
def test_figura_circular_desde_el_tanglegram(datos_sinteticos, tmp_path):
    d = datos_sinteticos
    from tanglegram_utils.models import Entradas
    ent = Entradas(str(d / "v.newick"), str(d / "g.newick"), tabla=str(d / "tabla.csv"))
    est = EstiloFigura(forma="circular", dpi=40)
    res = ejecutar_en_directorio(ent, OpcionesAnalisis(max_rebeldes=0), OpcionesColapso(), est,
                                 str(tmp_path / "ind"), figura="individuales")
    assert [os.path.basename(i) for i in res.imagenes] == ["tanglegram_gtdb.png", "tanglegram_vamphyre.png"]
    with pytest.raises(ErrorDatos, match="circular"):
        ejecutar_en_directorio(ent, OpcionesAnalisis(max_rebeldes=0), OpcionesColapso(), est,
                               str(tmp_path / "tg"), figura="tanglegram")
