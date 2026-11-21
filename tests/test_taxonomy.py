import numpy as np
import pytest

from tanglegram_utils.constants import GRIS_OTROS, GRIS_SIN, PALETA, SIN
from tanglegram_utils.errors import ErrorDatos
from tanglegram_utils.taxonomy import (conj_taxa, construir_esquema, etiqueta_consenso,
                                       genero_base, parsear_tax, rango_letra, sin_prefijo,
                                       tokens_de_hoja)


def test_parsear_tax_ignora_rangos_vacios():
    d = parsear_tax("d__Archaea;p__P1;c__;g__Halorubrum;s__Halorubrum sp. x")
    assert "c" not in d and d["g"] == "g__Halorubrum" and d["s"] == "s__Halorubrum sp. x"


def test_tokens_y_conj_taxa():
    tax = [parsear_tax("p__P1;g__Hal"), parsear_tax("p__P2;g__Hal"), parsear_tax("p__P2;g__Nat")]
    toks = {i: tokens_de_hoja(t) for i, t in enumerate(tax)}
    assert conj_taxa(toks, "g__hal") == {0, 1}
    assert conj_taxa(toks, "Nat, P1") == {0, 2}
    assert conj_taxa(toks, None) is None


def test_rango_letra():
    assert rango_letra("familia") == "f" and rango_letra("G") == "g" and rango_letra("género") == "g"
    assert rango_letra("family") == "f" and rango_letra("Genus") == "g"
    with pytest.raises(ErrorDatos):
        rango_letra("tribu")


def test_sin_prefijo_y_genero_base():
    assert sin_prefijo("g__Halorubrum") == "Halorubrum"
    assert genero_base("g__Thermococcus_A") == "Thermococcus"


def test_etiqueta_consenso():
    tax = [parsear_tax("f__F1;g__A"), parsear_tax("f__F1;g__B"), parsear_tax("f__F1;g__B"), parsear_tax("f__F2;g__C")]
    assert etiqueta_consenso(np.array([1, 2]), tax)[0] == "g__B"
    assert etiqueta_consenso(np.array([0, 1, 2]), tax)[0] == "f__F1"


def test_esquema_color_orden_fijo_y_otros():
    tax = [parsear_tax(f"f__{f}") for f in ["A"] * 5 + ["B"] * 3 + ["C"] * 2 + [""]]
    e = construir_esquema(tax, "f", top=2)
    assert e.orden == ["f__A", "f__B", "f__C"]            # by frequency
    assert e.col_cat["f__A"] == PALETA[0] and e.col_cat["f__B"] == PALETA[1]
    assert e.col_cat["f__C"] == GRIS_OTROS and e.col_cat[SIN] == GRIS_SIN
    assert e.categoria_pura(np.array([0, 1, 2])) == "f__A"
    assert e.categoria_pura(np.array([0, 5])) is None      # mixed
    assert e.categoria_pura(np.array([10])) is None        # unclassified


def test_esquema_sin_tabla():
    e = construir_esquema([{}, {}, {}], "f", top=8)
    assert not e.hay_tax and set(e.cat) == {SIN}
