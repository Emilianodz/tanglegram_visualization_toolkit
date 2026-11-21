"""Orchestration: loading/harmonisation -> tree comparison -> figure and reports.

Each stage is an independent function with explicit inputs/outputs, so it can be
reused (or tested) separately from other scripts."""
import logging
import os

import numpy as np

from .circular import DibujanteCircular
from .collapse import OpcionesColapso, planificar_colapsos
from .errors import ErrorDatos
from .idmap import cargar_tabla, homogeneizar
from .models import Analisis, Datos, OpcionesAnalisis, Visibilidad
from .ordering import contar_cruces, desenredar
from .plotting import ContenidoFigura, DibujanteTanglegram, EstiloFigura, guardar_figura
from .report import escribir_tablas_y_reporte
from .rooting import enraizar_balanceado, enraizar_en_split, raiz_comun
from .splits import (ancestros_concordantes, asignar_ids, conj_splits, info_arbol,
                     rebeldes_voraz)
from .taxonomy import (conj_taxa, construir_esquema, parsear_tax, rango_letra,
                       tokens_de_hoja)
from .tree import cargar_arbol, hojas_de, indexar, podar

log = logging.getLogger(__name__)


def preparar_datos(entradas, prefijo, color_rango="familia", color_top=8):
    """Read trees and table, harmonise names and build the colour scheme."""
    log.info("[-] Loading and cleaning trees...")
    rV = cargar_arbol(entradas.tree_vamphyre)
    rG = cargar_arbol(entradas.tree_gtdb)
    df = cargar_tabla(entradas.tabla, entradas.col_id, entradas.col_tax) if entradas.tabla else None
    rV, rG, clasif, _ = homogeneizar(rV, rG, df, prefijo, mapa_ids=entradas.mapa_ids,
                                     fasta_dir=entradas.fasta_dir,
                                     tree_gtdb_original=entradas.tree_gtdb_original)
    nombres = sorted(n.nombre for n in hojas_de(rV))
    idx = {nm: i for i, nm in enumerate(nombres)}
    tax = [parsear_tax(clasif.get(nm, "")) for nm in nombres]
    toks = {i: tokens_de_hoja(tax[i]) for i in range(len(nombres))}
    esquema = construir_esquema(tax, rango_letra(color_rango), color_top)
    return Datos(rV=rV, rG=rG, nombres=nombres, idx=idx, tax=tax, toks=toks, esquema=esquema)


def analizar(datos, op):
    """Compare the trees: RF, per-node discordance, rogue genomes and branch order."""
    rV, rG, N0 = datos.rV, datos.rG, datos.N0
    conservar = set(range(N0))
    if op.solo_taxa:
        conservar = conj_taxa(datos.toks, op.solo_taxa)
        if len(conservar) < 4:
            raise ErrorDatos("--solo-taxa leaves fewer than 4 genomes.")
        log.info(f"    [i] --solo-taxa: analysing {len(conservar)} of {N0} genomes")
        nom_cons = {datos.nombres[i] for i in conservar}
        rV, rG = podar(rV, nom_cons), podar(rG, nom_cons)
    todo = 0
    for i in conservar:
        todo |= 1 << i
    N = len(conservar)

    # common rooting (visual only; RF does not depend on the root)
    indexar(rV, datos.idx)
    indexar(rG, datos.idx)
    sp, motivo = raiz_comun(rV, rG, todo, datos.idx, datos.toks,
                            op.raiz_taxon.lower() if op.raiz_taxon else None)
    if sp is None:
        log.warning("    [!] No common split to root both trees the same way; each tree is rooted on its most balanced split.")
        rV, rG = enraizar_balanceado(rV, todo, datos.idx), enraizar_balanceado(rG, todo, datos.idx)
    else:
        rV, rG = enraizar_en_split(rV, sp, todo, datos.idx), enraizar_en_split(rG, sp, todo, datos.idx)
    log.info(f"[+] Common display root: {motivo}")

    sV, sG = conj_splits(rV, todo), conj_splits(rG, todo)
    rf = len(sV ^ sG)
    nrf = rf / (len(sV) + len(sG)) if (len(sV) + len(sG)) else 0.0
    log.info(f"[+] {N} genomes | splits VAMPhyRE {len(sV)}, GTDB {len(sG)}, shared {len(sV & sG)} | RF = {rf} (normalised {nrf:.3f})")
    infoV = info_arbol(rV, rG, sG, todo)
    infoG = info_arbol(rG, rV, sV, todo)
    idsV, idsG = asignar_ids(rV, "V"), asignar_ids(rG, "G")
    ancV, ancG = ancestros_concordantes(rV, infoV), ancestros_concordantes(rG, infoG)

    pasos, loo = [], {}
    if op.max_rebeldes > 0:
        log.info(f"[-] Searching for rogue genomes (greedy removal, up to {op.max_rebeldes} steps)...")
        pasos, loo = rebeldes_voraz(list(sV), list(sG), todo, op.max_rebeldes)
        log.info(f"    [+] RF {rf} -> {pasos[-1][2] if pasos else rf} after removing {len(pasos)} genomes")

    if op.sin_desenredar:
        cr0 = cr1 = contar_cruces(rV, rG, N0)
    else:
        log.info("[-] Untangling (reordering branches to minimise crossings)...")
        cr0, cr1 = desenredar(rV, rG, N0)
        log.info(f"    [+] Line crossings: {cr0} -> {cr1}")

    return Analisis(rV=rV, rG=rG, conservar=conservar, todo=todo, N=N, motivo_raiz=motivo,
                    sV=sV, sG=sG, rf=rf, nrf=nrf, infoV=infoV, infoG=infoG, idsV=idsV, idsG=idsG,
                    ancV=ancV, ancG=ancG, pasos=pasos, loo=loo, cr0=cr0, cr1=cr1,
                    sin_desenredar=op.sin_desenredar)


def calcular_visibilidad(datos, conservar, resaltar=None, ocultar_lineas=None):
    """Which genomes are highlighted and which do not draw their connecting line."""
    N0 = datos.N0
    hl = np.ones(N0, dtype=bool)
    hay_res = bool(resaltar)
    if hay_res:
        sel = conj_taxa(datos.toks, resaltar)
        hl = np.zeros(N0, dtype=bool)
        hl[list(sel)] = True
        log.info(f"    [i] Highlighting {int(hl.sum())} genomes")
    vis_linea = np.ones(N0, dtype=bool)
    if ocultar_lineas:
        ocu = conj_taxa(datos.toks, ocultar_lineas)
        vis_linea[list(ocu)] = False
        log.info(f"    [i] Hidden lines: {len(ocu)} genomes")
    for i in range(N0):
        if i not in conservar:
            vis_linea[i] = False
    return Visibilidad(hl=hl, hay_res=hay_res, vis_linea=vis_linea)


def ejecutar(entradas, op=None, colapso=None, estilo=None, salida=None):
    """Full pipeline. Returns (datos, analisis)."""
    op = op or OpcionesAnalisis()
    colapso = colapso or OpcionesColapso()
    estilo = estilo or EstiloFigura()
    prefijo = salida.prefijo_efectivo
    if estilo.forma == "circular" and salida.figura == "tanglegram":
        raise ErrorDatos("The circular layout only draws single trees: use --figura vamphyre, gtdb or individuales.")
    os.makedirs(os.path.dirname(os.path.abspath(prefijo)) or ".", exist_ok=True)

    datos = preparar_datos(entradas, prefijo, op.color_rango, op.color_top)
    an = analizar(datos, op)
    colV, colG = planificar_colapsos(an.rV, an.rG, an.infoV, an.infoG, an.idsV, an.idsG,
                                     datos.tax, datos.toks, an.todo, colapso)
    vis = calcular_visibilidad(datos, an.conservar, op.resaltar, op.ocultar_lineas)

    contenido = ContenidoFigura(rV=an.rV, rG=an.rG, colV=colV, colG=colG, esquema=datos.esquema,
                                tax=datos.tax, nombres=datos.nombres, hl=vis.hl, hay_res=vis.hay_res,
                                vis_linea=vis.vis_linea, N=an.N, N0=datos.N0, rf=an.rf, nrf=an.nrf,
                                cr0=an.cr0, cr1=an.cr1, modo_colapso=colapso.modo,
                                conservar=an.conservar)
    for ruta, panel in salida.figuras():
        dibujante = DibujanteCircular if estilo.forma == "circular" else DibujanteTanglegram
        fig = dibujante(contenido, estilo, panel).dibujar()
        guardar_figura(fig, ruta, estilo.dpi, salida.formato_extra)

    escribir_tablas_y_reporte(datos, an, prefijo, op.max_rebeldes, op.rango_grupos, op.min_grupo,
                              colV, colG, colapso.modo)
    return datos, an
