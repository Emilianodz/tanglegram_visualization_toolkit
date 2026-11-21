"""Data containers passed between pipeline stages."""
import os
from dataclasses import dataclass

import numpy as np

from .errors import ErrorDatos

FIGURAS = ("tanglegram", "vamphyre", "gtdb", "individuales")


@dataclass
class Entradas:
    tree_vamphyre: str
    tree_gtdb: str
    tabla: str = None
    col_id: str = "user_genome"
    col_tax: str = "classification"
    tree_gtdb_original: str = None
    fasta_dir: str = None
    mapa_ids: str = None


@dataclass
class OpcionesAnalisis:
    color_rango: str = "familia"
    color_top: int = 8
    resaltar: str = None
    ocultar_lineas: str = None
    solo_taxa: str = None
    raiz_taxon: str = None
    sin_desenredar: bool = False
    max_rebeldes: int = 25
    rango_grupos: str = "familia,genero"
    min_grupo: int = 3


@dataclass
class EntradaArbol:
    """A single tree (e.g. VAMPhyRE, GTDB-Tk or a future consensus) and, for colouring
    by taxonomy, the table that describes it."""
    arbol: str
    tabla: str = None
    col_id: str = "user_genome"
    col_tax: str = "classification"
    mapa_ids: str = None            # TSV table_id <tab> leaf_name (e.g. the tanglegram *_mapa_ids.tsv)
    fasta_dir: str = None
    arbol_referencia: str = None    # tree with the SAME topology and the table names (the original GTDB tree)
    titulo: str = None


@dataclass
class OpcionesArbol:
    color_rango: str = "familia"
    color_top: int = 8
    resaltar: str = None
    solo_taxa: str = None
    raiz: str = "original"          # original | punto-medio | taxon
    raiz_taxon: str = None          # outgroup when raiz == 'taxon'
    contraer_soporte: float = None  # nodes with support < value (0-100) become polytomies
    orden: str = "ladderize"        # ladderize | original
    rango_grupos: str = "familia,genero"
    min_grupo: int = 3
    itol: bool = True               # export tree + editable annotations for iTOL (_itol.*)


@dataclass
class Salida:
    out_plot: str
    prefijo: str = None
    formato_extra: str = None
    figura: str = "tanglegram"      # tanglegram | vamphyre | gtdb | individuales (one PNG per tree)

    @property
    def prefijo_efectivo(self):
        return self.prefijo or os.path.splitext(self.out_plot)[0]

    def figuras(self):
        """[(path, panel)] to draw according to `figura`."""
        if self.figura not in FIGURAS:
            raise ErrorDatos(f"Unknown figure type: {self.figura} (use {', '.join(FIGURAS)})")
        if self.figura == "individuales":
            base, ext = os.path.splitext(self.out_plot)
            return [(f"{base}_vamphyre{ext}", "vamphyre"), (f"{base}_gtdb{ext}", "gtdb")]
        return [(self.out_plot, {"tanglegram": "ambos"}.get(self.figura, self.figura))]


@dataclass
class Datos:
    """Harmonised trees + taxonomy per leaf."""
    rV: object
    rG: object
    nombres: list              # name of each leaf (index = bit)
    idx: dict                  # name -> index
    tax: list                  # dict rank->taxon per leaf
    toks: dict                 # index -> set of search tokens
    esquema: object            # EsquemaColor

    @property
    def N0(self):
        return len(self.nombres)

    def corto(self, i):
        """Readable leaf name: GTDB species or truncated name."""
        s = self.tax[i].get("s")
        return s[3:] if s else self.nombres[i][:40]


@dataclass
class Analisis:
    """Results of the tree comparison."""
    rV: object
    rG: object
    conservar: set
    todo: int
    N: int
    motivo_raiz: str
    sV: set
    sG: set
    rf: int
    nrf: float
    infoV: dict
    infoG: dict
    idsV: dict
    idsG: dict
    ancV: dict
    ancG: dict
    pasos: list
    loo: dict
    cr0: int
    cr1: int
    sin_desenredar: bool = False


@dataclass
class Visibilidad:
    hl: np.ndarray             # highlighted genomes
    hay_res: bool              # was --resaltar requested?
    vis_linea: np.ndarray      # genomes whose connecting lines are drawn
