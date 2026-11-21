from tanglegram_utils.names import acc_gcf, limpiar_etiqueta


def test_limpiar_etiqueta():
    assert limpiar_etiqueta(" >Hal rub=x ") == "Hal_rub_x"
    assert limpiar_etiqueta("[ 12 ] foo__bar") == "foo_bar"
    assert limpiar_etiqueta("a#b") == "ab"


def test_acc_gcf():
    assert acc_gcf("GCF_000007185.1_ASM718v1_genomic") == "GCF_000007185.1"
    assert acc_gcf("NZ_CP014855.1_Thermococcus") is None
