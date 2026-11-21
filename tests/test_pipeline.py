import pytest

from tanglegram_utils.cli import main
from tanglegram_utils.models import Entradas, OpcionesAnalisis, Salida
from tanglegram_utils.pipeline import ejecutar


def _entradas(d, **kw):
    return Entradas(tree_vamphyre=str(d / "v.newick"), tree_gtdb=str(d / "g.newick"),
                    tabla=str(d / "tabla.csv"), **kw)


def test_pipeline_completo_genera_todas_las_salidas(datos_sinteticos):
    d = datos_sinteticos
    datos, an = ejecutar(_entradas(d), OpcionesAnalisis(color_rango="familia"),
                         salida=Salida(out_plot=str(d / "out" / "t.png")))
    assert an.N == 12 and an.rf > 0 and 0 < an.nrf <= 1
    for suf in ("t.png", "t_reporte.md", "t_hojas_rebeldes.csv", "t_nodos_vamphyre.csv",
                "t_nodos_gtdb.csv", "t_grupos.csv", "t_mapa_ids.tsv"):
        assert (d / "out" / suf).exists(), suf
    assert "RF" in (d / "out" / "t_reporte.md").read_text(encoding="utf-8")


def test_rf_pipeline_coincide_con_dendropy(datos_sinteticos):
    dendropy = pytest.importorskip("dendropy")
    from dendropy.calculate import treecompare
    d = datos_sinteticos
    ns = dendropy.TaxonNamespace()
    kw = dict(schema="newick", taxon_namespace=ns, rooting="force-unrooted")
    a, b = (dendropy.Tree.get(path=str(d / n), **kw) for n in ("v.newick", "g.newick"))
    a.encode_bipartitions()
    b.encode_bipartitions()
    _, an = ejecutar(_entradas(d), OpcionesAnalisis(max_rebeldes=0),
                     salida=Salida(out_plot=str(d / "o" / "t.png")))
    assert an.rf == treecompare.symmetric_difference(a, b)


def test_solo_taxa_recalcula_sobre_el_subconjunto(datos_sinteticos):
    d = datos_sinteticos
    _, an = ejecutar(_entradas(d), OpcionesAnalisis(solo_taxa="f__Fa,f__Fb", max_rebeldes=0),
                     salida=Salida(out_plot=str(d / "o2" / "t.png")))
    assert an.N == 8


def test_sin_tabla_funciona_con_gradiente(datos_sinteticos):
    d = datos_sinteticos
    ent = Entradas(tree_vamphyre=str(d / "v.newick"), tree_gtdb=str(d / "g.newick"))
    _, an = ejecutar(ent, OpcionesAnalisis(max_rebeldes=2), salida=Salida(out_plot=str(d / "o3" / "t.png")))
    assert (d / "o3" / "t.png").exists() and an.N == 12


def test_cli_ok_y_error_de_datos(datos_sinteticos, capsys):
    d = datos_sinteticos
    base = ["--tree-vamphyre", str(d / "v.newick"), "--tree-gtdb", str(d / "g.newick"),
            "--tabla", str(d / "tabla.csv"), "--out-plot", str(d / "c" / "t.png"), "--max-rebeldes", "0"]
    assert main(base) == 0
    assert main(base + ["--col-id", "no_existe"]) == 1
    assert "no_existe" in capsys.readouterr().out
