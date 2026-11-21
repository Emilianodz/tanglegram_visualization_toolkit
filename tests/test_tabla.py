"""Reading taxonomy tables: separator detection and clear column errors."""
import pytest

from tanglegram_utils.errors import ErrorDatos
from tanglegram_utils.idmap import cargar_tabla, detectar_separador
from tanglegram_utils.servicio import adivinar_columnas, columnas_tabla

# GTDB-Tk style: tab-separated, but one header contains commas
CABECERA = ["names", "classification", "other_related_references(genome_id,species_name,radius,ANI,AF)", "note"]
FILAS = [["g1", "d__Archaea;p__P;f__Fa;g__Ga;s__Ga sp1", "GCF_1, s__X, 95.0", "x"],
         ["g2", "d__Archaea;p__P;f__Fb;g__Gb;s__Gb sp2", "N/A", "y"]]


@pytest.fixture
def tsv_con_comas(tmp_path):
    ruta = tmp_path / "gtdbtk.ar53.summary.tsv"
    ruta.write_text("\n".join("\t".join(f) for f in [CABECERA] + FILAS) + "\n")
    return ruta


def test_tab_wins_over_commas_in_header(tsv_con_comas):
    assert detectar_separador(tsv_con_comas) == "\t"
    assert columnas_tabla(tsv_con_comas) == CABECERA
    assert adivinar_columnas(tsv_con_comas) == ("names", "classification")


def test_load_with_custom_id_column(tsv_con_comas):
    df = cargar_tabla(tsv_con_comas, col_id="names")
    assert list(df.iloc[:, 0]) == ["g1", "g2"]


def test_missing_column_error_lists_columns_and_fix(tsv_con_comas):
    with pytest.raises(ErrorDatos) as e:
        cargar_tabla(tsv_con_comas)                       # default --col-id user_genome is absent
    msg = str(e.value)
    assert "'user_genome'" in msg and "--col-id" in msg
    assert "names, classification" in msg and "Likely fix: --col-id names" in msg


def test_plain_csv_and_semicolon(tmp_path):
    for sep in (",", ";"):
        ruta = tmp_path / f"t{ord(sep)}.csv"
        ruta.write_text(f"user_genome{sep}classification\ng1{sep}d__A;g__G\n" if sep == "," else
                        f"user_genome{sep}classification\ng1{sep}d__A\n")
        assert detectar_separador(ruta) == sep
        assert len(cargar_tabla(ruta)) == 1
