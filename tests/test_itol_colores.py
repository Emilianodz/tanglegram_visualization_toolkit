"""iTOL export and extended palette (more than 8 colours)."""
import logging
import os
import zipfile

from tanglegram_utils.collapse import OpcionesColapso
from tanglegram_utils.constants import GRIS_OTROS, PALETA, PALETA_EXTRA
from tanglegram_utils.itol import id_nodo_itol
from tanglegram_utils.models import EntradaArbol, OpcionesArbol
from tanglegram_utils.plotting import EstiloFigura
from tanglegram_utils.servicio import ejecutar_arbol_en_directorio
from tanglegram_utils.taxonomy import aviso_colores, construir_esquema, paleta, parsear_tax
from tanglegram_utils.tree import desde_cadena, hojas_de

NWK = ("(((a1:1,a2:1)0.99:1,(a3:1,a4:1)0.5:1)0.95:1,((b1:1,b2:1)1.0:1,(b3:1,b4:1)0.4:1)0.9:1,"
       "((c1:1,c2:1)1.0:1,(c3:1,c4:1)0.97:1)0.88:1);")


# ---------------------------------------------------------------- palette
def test_paleta_conserva_los_20_primeros_y_crece_sin_repetir():
    assert paleta(8) == PALETA
    assert paleta(20) == PALETA + PALETA_EXTRA
    p = paleta(80)
    assert p[:20] == PALETA + PALETA_EXTRA
    assert len(p) == 80 and len(set(p)) == 80 and GRIS_OTROS not in p


def test_aviso_colores_por_niveles():
    assert aviso_colores(8) is None
    assert "colour blindness" in aviso_colores(12)
    assert "hard to tell apart" in aviso_colores(25)


def test_esquema_con_mas_de_20_taxones(caplog):
    tax = [parsear_tax(f"d__A;g__G{k:02d}") for k in range(30) for _ in range(2)]
    with caplog.at_level(logging.WARNING):
        e = construir_esquema(tax, "g", 30)
    assert len({e.col_cat[c] for c in e.orden}) == 30 and GRIS_OTROS not in e.col_cat.values()
    assert "hard to tell apart" in caplog.text
    e = construir_esquema(tax, "g", 25)
    assert sum(e.col_cat[c] == GRIS_OTROS for c in e.orden) == 5


# ---------------------------------------------------------------- iTOL
def test_id_nodo_itol_usa_ancestro_comun():
    r = desde_cadena(NWK)
    assert id_nodo_itol(r.hijos[0]) == "a1|a4"
    assert id_nodo_itol(hojas_de(r)[0]) == "a1"


def test_exportacion_itol(datos_sinteticos, tmp_path):
    (datos_sinteticos / "sop.newick").write_text(NWK + "\n")
    ent = EntradaArbol(arbol=str(datos_sinteticos / "sop.newick"), tabla=str(datos_sinteticos / "tabla.csv"))
    col = OpcionesColapso(modo="taxonomia", rango="f", min_colapso=2)
    res = ejecutar_arbol_en_directorio(ent, OpcionesArbol(color_rango="familia"), col,
                                       EstiloFigura(forma="rectangular", dpi=40), str(tmp_path / "o"))
    for suf in ("_itol.tree", "_itol_colores.txt", "_itol_franja.txt", "_itol_etiquetas.txt",
                "_itol_colapso.txt", "_itol.zip"):
        assert res.con_sufijo(suf), suf
    arbol = open(res.con_sufijo("_itol.tree")[0]).read()
    assert hojas_de(desde_cadena(arbol)) and "99" in arbol             # 0-1 support arrives as 0-100
    colores = open(res.con_sufijo("_itol_colores.txt")[0]).read().splitlines()
    assert colores[:3] == ["TREE_COLORS", "SEPARATOR TAB", "DATA"]
    assert f"a1\tlabel\t{PALETA[0]}\tnormal\t1" in colores or f"a1\tlabel\t{PALETA[1]}\tnormal\t1" in colores
    assert any("\tclade\t" in l and l.split("\t")[0].count("|") == 1 for l in colores)
    franja = open(res.con_sufijo("_itol_franja.txt")[0]).read()
    assert franja.startswith("DATASET_COLORSTRIP\nSEPARATOR TAB\n") and "LEGEND_LABELS\tF" in franja
    colapso = open(res.con_sufijo("_itol_colapso.txt")[0]).read().splitlines()
    assert colapso[:2] == ["COLLAPSE", "DATA"] and len(colapso) == 5     # 3 families + header
    etiquetas = open(res.con_sufijo("_itol_etiquetas.txt")[0]).read()
    assert "a1\tGa sp1" in etiquetas
    with zipfile.ZipFile(res.con_sufijo("_itol.zip")[0]) as z:
        nombres = z.namelist()
    assert "arbol_itol.tree" in nombres and all(n.endswith((".tree", ".txt")) for n in nombres)


def test_sin_itol(datos_sinteticos, tmp_path):
    (datos_sinteticos / "sop.newick").write_text(NWK + "\n")
    ent = EntradaArbol(arbol=str(datos_sinteticos / "sop.newick"))
    res = ejecutar_arbol_en_directorio(ent, OpcionesArbol(itol=False), OpcionesColapso(),
                                       EstiloFigura(dpi=40), str(tmp_path / "o"))
    assert not [a for a in res.archivos if "_itol" in os.path.basename(a)]
    assert "--sin-itol" in res.comando
