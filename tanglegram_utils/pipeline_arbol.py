"""Single-tree pipeline: one Newick (+ optional taxonomy table) -> rectangular or
circular figure with colour, collapse, highlighting and support, plus its tables.

Does not compare against another tree: works for VAMPhyRE, GTDB-Tk or a future consensus."""
import logging
import os
from dataclasses import dataclass

import numpy as np

from .circular import DibujanteCircular
from .collapse import OpcionesColapso, planificar_colapsos_arbol
from .errors import ErrorDatos
from .itol import exportar_itol
from .idmap import cargar_tabla, homogeneizar_arbol
from .models import Datos, OpcionesArbol
from .ordering import ladderize
from .pipeline import calcular_visibilidad
from .plotting import ContenidoFigura, DibujanteTanglegram, EstiloFigura, guardar_figura
from .report_arbol import (escribir_reporte_arbol, rangos_de, tabla_grupos_arbol, tabla_nodos_arbol)
from .rooting import enraizar_en_grupo, enraizar_punto_medio
from .splits import asignar_ids
from .support import contraer_soporte_bajo, normalizar_soporte, resumen_soporte
from .taxonomy import conj_taxa, construir_esquema, parsear_tax, rango_letra, tokens_de_hoja
from .tree import a_newick, cargar_arbol, hojas_de, indexar, podar

log = logging.getLogger(__name__)


@dataclass
class AnalisisArbol:
    datos: Datos
    raiz: object
    N: int
    motivo_raiz: str
    resumen_soporte: dict
    n_contraidos: int
    conservar: set


def ejecutar_arbol(entrada, op=None, colapso=None, estilo=None, salida=None):
    """Draw and describe a single tree. Returns an `AnalisisArbol`."""
    op = op or OpcionesArbol()
    colapso = colapso or OpcionesColapso()
    estilo = estilo or EstiloFigura(forma="circular")
    if op.raiz not in ("original", "punto-medio", "taxon"):
        raise ErrorDatos(f"Unknown root option: {op.raiz} (use original, punto-medio or taxon)")
    if op.orden not in ("ladderize", "original"):
        raise ErrorDatos(f"Unknown order option: {op.orden} (use ladderize or original)")
    prefijo = salida.prefijo_efectivo
    os.makedirs(os.path.dirname(os.path.abspath(prefijo)) or ".", exist_ok=True)

    log.info("[-] Loading tree...")
    raiz = cargar_arbol(entrada.arbol)
    if len(hojas_de(raiz)) < 4:
        raise ErrorDatos("The tree has fewer than 4 leaves.")
    factor = normalizar_soporte(raiz)
    resumen = resumen_soporte(raiz)
    if resumen:
        log.info(f"    [i] Support: {resumen['con_soporte']} of {resumen['nodos_internos']} internal nodes"
                 + (" (0-1 rescaled to 0-100)" if factor != 1.0 else "")
                 + f"; median {resumen['mediana']:.1f}, ≥70: {resumen['pct_ge_70']:.0f}%, ≥90: {resumen['pct_ge_90']:.0f}%")
    else:
        log.info("    [i] The tree has no support values on internal nodes")

    df = cargar_tabla(entrada.tabla, entrada.col_id, entrada.col_tax) if entrada.tabla else None
    clasif = {}
    if df is not None:
        clasif = homogeneizar_arbol(raiz, df, prefijo, mapa_ids=entrada.mapa_ids, fasta_dir=entrada.fasta_dir,
                                    arbol_referencia=entrada.arbol_referencia)
    nombres = sorted(n.nombre for n in hojas_de(raiz))
    idx = {nm: i for i, nm in enumerate(nombres)}
    tax = [parsear_tax(clasif.get(nm, "")) for nm in nombres]
    toks = {i: tokens_de_hoja(tax[i]) for i in range(len(nombres))}
    esquema = construir_esquema(tax, rango_letra(op.color_rango), op.color_top)   # colours computed on the WHOLE tree
    datos = Datos(rV=raiz, rG=None, nombres=nombres, idx=idx, tax=tax, toks=toks, esquema=esquema)

    conservar = set(range(len(nombres)))
    if op.solo_taxa:
        conservar = conj_taxa(toks, op.solo_taxa)
        if len(conservar) < 4:
            raise ErrorDatos("--solo-taxa leaves fewer than 4 genomes.")
        log.info(f"    [i] --solo-taxa: drawing {len(conservar)} of {len(nombres)} genomes")
        raiz = podar(raiz, {nombres[i] for i in conservar})
    todo = 0
    for i in conservar:
        todo |= 1 << i
    indexar(raiz, idx)

    # root (changes only the drawing and which groups appear as clades)
    motivo = "as in file"
    if op.raiz == "punto-medio":
        raiz = enraizar_punto_medio(raiz)
        indexar(raiz, idx)
        motivo = "midpoint"
    elif op.raiz == "taxon":
        sel = conj_taxa(toks, op.raiz_taxon) & conservar if op.raiz_taxon else set()
        if not sel:
            raise ErrorDatos("--raiz taxon requires --raiz-taxon with a group present in the tree.")
        mask = 0
        for i in sel:
            mask |= 1 << i
        raiz, ok = enraizar_en_grupo(raiz, mask, todo, idx)
        if ok:
            motivo = f"outgroup {op.raiz_taxon}"
        else:
            log.warning(f"    [!] {op.raiz_taxon} does not form a clade in this tree: keeping the root from the file.")

    n_contraidos = 0
    if op.contraer_soporte is not None:
        raiz, n_contraidos = contraer_soporte_bajo(raiz, op.contraer_soporte)
        log.info(f"[+] Nodes with support < {op.contraer_soporte:g} collapsed into polytomies: {n_contraidos}")
        indexar(raiz, idx)
    if op.orden == "ladderize":
        ladderize(raiz)
    ids = asignar_ids(raiz, "N")

    col = planificar_colapsos_arbol(raiz, ids, tax, toks, todo, colapso)
    vis = calcular_visibilidad(datos, conservar, op.resaltar, None)
    titulo = entrada.titulo or os.path.splitext(os.path.basename(entrada.arbol))[0]
    subtitulo = f"{len(conservar)} genomes · root: {motivo}"
    contenido = ContenidoFigura(rV=raiz, rG=None, colV=col, colG=set(), esquema=esquema, tax=tax, nombres=nombres,
                                hl=vis.hl, hay_res=vis.hay_res, vis_linea=vis.vis_linea, N=len(conservar),
                                N0=len(nombres), rf=0, nrf=0.0, cr0=0, cr1=0, modo_colapso=colapso.modo,
                                conservar=conservar, titulo=titulo, subtitulo=subtitulo)
    dibujante = DibujanteCircular if estilo.forma == "circular" else DibujanteTanglegram
    fig = dibujante(contenido, estilo, "vamphyre").dibujar()
    guardar_figura(fig, salida.out_plot, estilo.dpi, salida.formato_extra)

    # tables and supporting files
    with open(prefijo + "_arbol.nwk", "w", encoding="utf-8") as f:
        f.write(a_newick(raiz) + "\n")
    tabla_nodos_arbol(raiz, tax, esquema.hay_tax).to_csv(prefijo + "_nodos.csv", index=False)
    dfg = None
    if esquema.hay_tax:
        dfg = tabla_grupos_arbol(raiz, tax, conservar, rangos_de(op.rango_grupos), op.min_grupo)
        dfg.to_csv(prefijo + "_grupos.csv", index=False)
    escribir_reporte_arbol(prefijo + "_reporte.md", titulo, len(conservar), motivo, resumen, n_contraidos,
                           op.contraer_soporte, dfg, esquema.hay_tax)
    log.info(f"[+] Tables and report: {prefijo}_nodos.csv, _grupos.csv, _reporte.md, _arbol.nwk")
    if op.itol:
        exportar_itol(raiz, datos, prefijo, colapsados=col, hl=vis.hl, hay_res=vis.hay_res)
    return AnalisisArbol(datos=datos, raiz=raiz, N=len(conservar), motivo_raiz=motivo, resumen_soporte=resumen,
                         n_contraidos=n_contraidos, conservar=conservar)


def n_con_taxonomia(analisis):
    """Number of drawn genomes that have taxonomy."""
    return int(np.sum([bool(analisis.datos.tax[i]) for i in analisis.conservar]))
