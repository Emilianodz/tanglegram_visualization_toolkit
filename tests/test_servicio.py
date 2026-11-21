import os
import shlex

from tanglegram_utils.cli import comando_equivalente, construir_parser, opciones_desde_args
from tanglegram_utils.collapse import OpcionesColapso
from tanglegram_utils.errors import ErrorDatos
from tanglegram_utils.models import Entradas, OpcionesAnalisis, Salida
from tanglegram_utils.plotting import EstiloFigura
from tanglegram_utils.servicio import (adivinar_columnas, columnas_tabla, diagnosticar_ids,
                                       ejecutar_en_directorio, empaquetar_zip, listar_taxa)

import pytest


def _entradas(d):
    return Entradas(str(d / "v.newick"), str(d / "g.newick"), tabla=str(d / "tabla.csv"))


def test_comando_por_defecto_solo_lleva_obligatorios():
    cmd = comando_equivalente(Entradas("v.nwk", "g.nwk"), OpcionesAnalisis(), OpcionesColapso(),
                              EstiloFigura(), Salida(out_plot="o.png"))
    assert cmd.count("--") == 3


def test_comando_ida_y_vuelta():
    """The generated command, fed through the real parser, rebuilds the same options."""
    ent = Entradas("v.nwk", "g.nwk", tabla="mi tabla.xlsx", col_id="id", tree_gtdb_original="o.tree")
    op = OpcionesAnalisis(color_rango="genero", color_top=10, resaltar="g__A,g__B", solo_taxa="f__X",
                          max_rebeldes=0, sin_desenredar=True)
    col = OpcionesColapso(modo="taxonomia", rango="f", min_colapso=5, expandir_taxa="g__A")
    est = EstiloFigura(etiquetas_hojas=True, etiquetas_grupos=False, dpi=300, alto=40.0, escala="longitud", ancho_arbol=6.0)
    sal = Salida(out_plot="salida/x.png", formato_extra="svg,pdf", figura="gtdb")
    cmd = comando_equivalente(ent, op, col, est, sal)
    args = construir_parser().parse_args(shlex.split(cmd.replace("\\\n", " "))[2:])
    ent2, op2, col2, est2, sal2 = opciones_desde_args(args)
    assert (ent2, op2, col2, est2, sal2) == (ent, op, col, est, sal)


def test_listar_taxa_y_columnas(datos_sinteticos):
    tabla = str(datos_sinteticos / "tabla.csv")
    assert columnas_tabla(tabla) == ["user_genome", "classification"]
    taxa = listar_taxa(tabla)
    assert taxa["f"] == [("f__Fa", 4), ("f__Fb", 4), ("f__Fc", 4)]
    assert ("g__Ga", 4) in taxa["g"]


def test_ejecutar_en_directorio(datos_sinteticos, tmp_path):
    out = tmp_path / "salida"
    res = ejecutar_en_directorio(_entradas(datos_sinteticos), OpcionesAnalisis(max_rebeldes=2),
                                 OpcionesColapso(), EstiloFigura(dpi=40), str(out), formato_extra="svg")
    assert res.n_genomas == 12 and res.rf > 0 and res.n_con_taxonomia == 12
    assert res.imagen.endswith("tanglegram.png")
    assert res.con_sufijo(".svg") and res.con_sufijo("_reporte.md")
    assert "RF" in res.log                     # the pipeline log was captured
    assert "--out-plot salida/tanglegram.png" in res.comando
    assert "--tree-vamphyre v.newick" in res.comando      # without the temporary path
    assert len(empaquetar_zip(res.archivos)) > 1000


def test_error_de_datos_lleva_el_log(datos_sinteticos, tmp_path):
    ent = _entradas(datos_sinteticos)
    ent.col_id = "no_existe"
    with pytest.raises(ErrorDatos) as exc:
        ejecutar_en_directorio(ent, OpcionesAnalisis(), OpcionesColapso(), EstiloFigura(dpi=40),
                               str(tmp_path / "o"))
    assert hasattr(exc.value, "log")


def test_diagnostico_de_ids_y_columnas(datos_sinteticos, tmp_path):
    v, tabla = str(datos_sinteticos / "v.newick"), str(datos_sinteticos / "tabla.csv")
    assert adivinar_columnas(tabla) == ("user_genome", "classification")
    assert diagnosticar_ids(v, tabla) == 1.0
    # table with IDs that do not match the leaves -> coverage 0
    otra = tmp_path / "otra.csv"
    otra.write_text("id,tax\nGCF_1.1_x,d__A;p__B\nGCF_2.1_y,d__A;p__C\n")
    assert adivinar_columnas(str(otra)) == ("id", "tax")
    assert diagnosticar_ids(v, str(otra), "id", "tax") == 0.0


def test_tabla_sin_ligar_avisa_y_cuenta_cero(datos_sinteticos, tmp_path):
    otra = tmp_path / "otra.csv"
    filas = ["user_genome,classification"] + [f"GCF_{i}.1_x,d__A;p__B;c__C;o__D;f__E;g__F;s__G" for i in range(12)]
    otra.write_text("\n".join(filas) + "\n")
    ent = _entradas(datos_sinteticos)
    ent.tabla = str(otra)
    res = ejecutar_en_directorio(ent, OpcionesAnalisis(max_rebeldes=0), OpcionesColapso(),
                                 EstiloFigura(dpi=40), str(tmp_path / "o"))
    assert res.n_con_taxonomia == 0
    assert "[!] Only 0 of 12 genomes have taxonomy" in res.log


def test_figura_un_solo_arbol_e_individuales(datos_sinteticos, tmp_path):
    from PIL import Image
    ent, est = _entradas(datos_sinteticos), EstiloFigura(dpi=40)
    anchos = {}
    for fig in ("tanglegram", "vamphyre", "gtdb", "individuales"):
        res = ejecutar_en_directorio(ent, OpcionesAnalisis(max_rebeldes=0), OpcionesColapso(), est,
                                     str(tmp_path / fig), figura=fig)
        nombres = [os.path.basename(i) for i in res.imagenes]
        esperado = {"tanglegram": ["tanglegram.png"], "vamphyre": ["tanglegram.png"], "gtdb": ["tanglegram.png"],
                    "individuales": ["tanglegram_gtdb.png", "tanglegram_vamphyre.png"]}[fig]
        assert nombres == esperado
        assert f"--figura {fig}" in res.comando or fig == "tanglegram"
        anchos[fig] = Image.open(res.imagenes[0]).size[0]
    assert anchos["vamphyre"] < anchos["tanglegram"] and anchos["gtdb"] < anchos["tanglegram"]


def test_figura_invalida(datos_sinteticos, tmp_path):
    with pytest.raises(ErrorDatos):
        ejecutar_en_directorio(_entradas(datos_sinteticos), OpcionesAnalisis(), OpcionesColapso(),
                               EstiloFigura(dpi=40), str(tmp_path / "x"), figura="circular")
