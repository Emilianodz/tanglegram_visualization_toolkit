import pandas as pd
import pytest

from conftest import NWK_V
from tanglegram_utils.errors import ErrorDatos
from tanglegram_utils.idmap import alias_por_topologia, cargar_tabla, homogeneizar
from tanglegram_utils.tree import desde_cadena, hojas_de

# same tree as V but with GCF_*-style names
RENOMBRE = {f"{p}{k}": f"GCF_{ord(p) * 1000 + k:09d}.1_ASM_genomic" for p in "abc" for k in range(1, 5)}


def renombrar(nwk):
    for viejo, nuevo in sorted(RENOMBRE.items(), key=lambda kv: -len(kv[0])):
        nwk = nwk.replace(viejo + ":", nuevo + ":")
    return nwk


def test_alias_por_topologia_recupera_el_mapa():
    # distinct branch lengths per leaf so the matching is unambiguous
    orig = desde_cadena("(((x1:0.1,x2:0.2):1,(x3:0.3,x4:0.4):1):1,(x5:0.5,x6:0.6):1,(x7:0.7,x8:0.8):1);")
    canon = desde_cadena("(((n1:0.1,n2:0.2):1,(n3:0.3,n4:0.4):1):1,(n5:0.5,n6:0.6):1,(n7:0.7,n8:0.8):1);")
    alias, ambiguos = alias_por_topologia(orig, canon)
    assert alias == {f"x{i}": f"n{i}" for i in range(1, 9)} and not ambiguos


def test_alias_marca_ambiguos_si_hermanas_identicas():
    orig = desde_cadena("((x1:0,x2:0):1,x3:1,x4:1);")
    canon = desde_cadena("((n1:0,n2:0):1,n3:1,n4:1);")
    _, ambiguos = alias_por_topologia(orig, canon)
    assert {"x1", "x2"} <= ambiguos


def test_cargar_tabla_columna_ausente(tmp_path):
    f = tmp_path / "t.tsv"
    pd.DataFrame({"id": ["a"], "tax": ["d__x"]}).to_csv(f, sep="\t", index=False)
    with pytest.raises(ErrorDatos):
        cargar_tabla(str(f))
    assert list(cargar_tabla(str(f), "id", "tax").columns) == ["id", "tax"]


def test_homogeneizar_con_gtdb_en_gcf_y_arbol_original(tmp_path):
    # V (VAMPhyRE leaf names) and an already "clean" GTDB tree (same names); the table uses GCF names
    rV, rG = desde_cadena(NWK_V), desde_cadena(NWK_V)
    original = tmp_path / "orig.tree"
    original.write_text(renombrar(NWK_V))
    tabla = pd.DataFrame({"user_genome": [v.replace("_genomic", "_genomic") for v in RENOMBRE.values()],
                          "classification": [f"d__A;g__G{k[0]}" for k in RENOMBRE]})
    rV2, rG2, clasif, n = homogeneizar(rV, rG, tabla, str(tmp_path / "salida"),
                                       tree_gtdb_original=str(original))
    assert n == 12 and len(clasif) == 12
    assert clasif["a1"] == "d__A;g__Ga"
    assert (tmp_path / "salida_mapa_ids.tsv").exists()
    assert len(hojas_de(rV2)) == len(hojas_de(rG2)) == 12


def test_homogeneizar_pocas_hojas_comunes():
    rV, rG = desde_cadena("((a,b),(c,d),e);"), desde_cadena("((x,y),(z,w),v);")
    with pytest.raises(ErrorDatos):
        homogeneizar(rV, rG, None, "/tmp/no_importa")
