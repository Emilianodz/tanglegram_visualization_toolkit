"""Service layer: runs the pipeline inside a working directory and returns everything
an interface (web, notebook) needs, without depending on it.

It is the only thing `app.py` (Streamlit) uses besides the option dataclasses, so
the logic can be tested with pytest without starting any interface."""
import io
import logging
import os
import threading
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass, field, replace


from .cli import comando_equivalente
from .cli_arbol import comando_equivalente as comando_equivalente_arbol
from .constants import RANGOS
from .errors import ErrorDatos
from .idmap import cargar_tabla, leer_tabla
from .names import limpiar_etiqueta
from .models import Salida
from .pipeline import ejecutar
from .pipeline_arbol import ejecutar_arbol, n_con_taxonomia
from .support import normalizar_soporte, resumen_soporte
from .taxonomy import parsear_tax
from .tree import cargar_arbol, hojas_de


@dataclass
class Resultado:
    """Files and figures of a run."""
    directorio: str
    imagen: str
    archivos: list = field(default_factory=list)   # all outputs (absolute paths)
    log: str = ""
    comando: str = ""
    n_genomas: int = 0
    rf: int = 0
    rf_norm: float = 0.0
    cruces_antes: int = 0
    cruces_despues: int = 0
    n_con_taxonomia: int = 0     # displayed genomes that have taxonomy

    @property
    def imagenes(self):
        """Figure PNGs (tanglegram and/or single trees), excluding the rogue summary figure."""
        return [a for a in self.archivos if a.endswith(".png") and not a.endswith("_rebeldes.png")]

    def con_sufijo(self, sufijo):
        """Output paths whose name ends with `sufijo`."""
        return [a for a in self.archivos if a.endswith(sufijo)]


@dataclass
class ResultadoArbol(Resultado):
    """Result of a single-tree run (inherits files, log and command)."""
    motivo_raiz: str = ""
    n_contraidos: int = 0        # nodes collapsed into polytomies due to low support
    soporte: dict = None         # support summary (0-100 scale), or None if the tree has none


@contextmanager
def capturar_log():
    """Capture the package messages produced by THIS thread (several web sessions
    may run at the same time)."""
    buf = io.StringIO()
    hilo = threading.get_ident()
    handler = logging.StreamHandler(buf)
    handler.setFormatter(logging.Formatter("%(message)s"))
    handler.addFilter(lambda r: r.thread == hilo)
    logger = logging.getLogger("tanglegram_utils")
    nivel, propaga = logger.level, logger.propagate
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    try:
        yield buf
    finally:
        logger.removeHandler(handler)
        logger.setLevel(nivel)
        logger.propagate = propaga


def ejecutar_en_directorio(entradas, op, colapso, estilo, directorio, formato_extra=None, nombre="tanglegram",
                           figura="tanglegram"):
    """Run the pipeline writing every output into `directorio`.

    Returns a `Resultado`. Data errors (`ErrorDatos`) propagate; the log captured up
    to that point is stored in `e.log` so it can be shown to the user."""
    os.makedirs(directorio, exist_ok=True)
    salida = Salida(out_plot=os.path.join(directorio, f"{nombre}.png"), formato_extra=formato_extra, figura=figura)
    with capturar_log() as buf:
        try:
            datos, an = ejecutar(entradas, op, colapso, estilo, salida)
        except ErrorDatos as e:
            e.log = buf.getvalue()
            raise
    archivos = sorted(os.path.join(directorio, f) for f in os.listdir(directorio)
                      if f.startswith(nombre) and os.path.isfile(os.path.join(directorio, f)))
    figuras = [r for r, _ in salida.figuras()]
    return Resultado(directorio=directorio, imagen=figuras[0], archivos=archivos, log=buf.getvalue(),
                     comando=comando_equivalente(_con_nombres_base(entradas), op, colapso, estilo,
                                                 Salida(out_plot=f"salida/{nombre}.png",
                                                        formato_extra=formato_extra, figura=figura)),
                     n_genomas=an.N, rf=an.rf, rf_norm=an.nrf,
                     cruces_antes=an.cr0, cruces_despues=an.cr1,
                     n_con_taxonomia=sum(1 for i in an.conservar if datos.tax[i]))


def _con_nombres_base(entradas, campos=("tree_vamphyre", "tree_gtdb", "tabla", "tree_gtdb_original", "mapa_ids")):
    """Copy of `entradas` with only each file's name (not the temporary path), so the
    equivalent command is readable and easy to copy."""
    return replace(entradas, **{c: os.path.basename(getattr(entradas, c))
                                for c in campos if getattr(entradas, c)})


def ejecutar_arbol_en_directorio(entrada, op, colapso, estilo, directorio, formato_extra=None, nombre="arbol"):
    """Like `ejecutar_en_directorio`, but for a single tree (generate_tree.py)."""
    os.makedirs(directorio, exist_ok=True)
    salida = Salida(out_plot=os.path.join(directorio, f"{nombre}.png"), formato_extra=formato_extra)
    with capturar_log() as buf:
        try:
            an = ejecutar_arbol(entrada, op, colapso, estilo, salida)
        except ErrorDatos as e:
            e.log = buf.getvalue()
            raise
    archivos = sorted(os.path.join(directorio, f) for f in os.listdir(directorio)
                      if f.startswith(nombre) and os.path.isfile(os.path.join(directorio, f)))
    campos = ("arbol", "tabla", "mapa_ids", "arbol_referencia")
    comando = comando_equivalente_arbol(_con_nombres_base(entrada, campos), op, colapso, estilo,
                                        Salida(out_plot=f"salida/{nombre}.png", formato_extra=formato_extra))
    return ResultadoArbol(directorio=directorio, imagen=salida.out_plot, archivos=archivos, log=buf.getvalue(),
                          comando=comando, n_genomas=an.N, n_con_taxonomia=n_con_taxonomia(an),
                          motivo_raiz=an.motivo_raiz, n_contraidos=an.n_contraidos, soporte=an.resumen_soporte)


def empaquetar_zip(archivos):
    """In-memory ZIP of `archivos` (base names only inside the archive)."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for a in archivos:
            z.write(a, arcname=os.path.basename(a))
    return buf.getvalue()


def columnas_tabla(ruta):
    """Column names of an xlsx/csv/tsv table (to choose the ID and taxonomy columns)."""
    return [str(c).strip() for c in leer_tabla(ruta, nrows=5).columns]


def listar_taxa(ruta, col_id="user_genome", col_tax="classification"):
    """Table taxa by rank: {'g': [('g__Halorubrum', 31), ...], ...} (key = rank letter,
    in order d,p,c,o,f,g,s), from most to least frequent. Used to populate the
    interface selectors."""
    df = cargar_tabla(ruta, col_id, col_tax)
    cuentas = {}
    for cadena in df[col_tax]:
        for letra, tok in parsear_tax(cadena).items():
            cuentas.setdefault(letra, {}).setdefault(tok, 0)
            cuentas[letra][tok] += 1
    return {l: sorted(cuentas[l].items(), key=lambda kv: (-kv[1], kv[0]))
            for l in RANGOS if l in cuentas}


def adivinar_columnas(ruta):
    """Likely (col_id, col_tax): 'user_genome'/'classification' if present; otherwise
    the first column and the column that contains 'd__...' strings."""
    df = leer_tabla(ruta, nrows=20)
    df.columns = [str(c).strip() for c in df.columns]
    cols = list(df.columns)
    col_id = "user_genome" if "user_genome" in cols else cols[0]
    if "classification" in cols:
        return col_id, "classification"
    con_tax = [c for c in cols if df[c].astype(str).str.contains("d__", regex=False).any()]
    return col_id, (con_tax[0] if con_tax else cols[-1])


def diagnosticar_ids(tree_vamphyre, tabla, col_id="user_genome", col_tax="classification"):
    """Fraction of VAMPhyRE leaves whose name appears verbatim in the table (the case
    the pipeline links without help). If low, the original GTDB tree, an ID map or
    the FASTA folder is needed to link the taxonomy."""
    hojas = [n.nombre for n in hojas_de(cargar_arbol(tree_vamphyre))]
    ids = {limpiar_etiqueta(g) for g in cargar_tabla(tabla, col_id, col_tax)[col_id]}
    return sum(1 for h in hojas if h in ids) / len(hojas) if hojas else 0.0


def inspeccionar_arbol(ruta):
    """(leaf count, 0-100 support summary or None): what the interface shows when a tree is loaded."""
    raiz = cargar_arbol(ruta)
    normalizar_soporte(raiz)
    return len(hojas_de(raiz)), resumen_soporte(raiz)
