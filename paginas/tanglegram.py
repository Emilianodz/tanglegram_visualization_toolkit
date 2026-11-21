"""App page "Tanglegram": compares the VAMPhyRE tree with the GTDB-Tk tree (or draws a single one).

Opened from `app.py` (`streamlit run app.py`). Interface only: all logic lives in `tanglegram_utils`
(this page calls `servicio.ejecutar_en_directorio`, just as the CLI calls `pipeline`). Files are
processed in a local temporary folder; nothing leaves your computer."""
import os
import shutil

import pandas as pd
import streamlit as st

from tanglegram_utils import __version__
from tanglegram_utils.collapse import OpcionesColapso
from tanglegram_utils.errors import ErrorDatos
from tanglegram_utils.models import Entradas, OpcionesAnalisis
from tanglegram_utils.plotting import EstiloFigura
from tanglegram_utils.servicio import ejecutar_en_directorio, empaquetar_zip
from tanglegram_utils.constants import MAX_DISTINGUIBLE
from tanglegram_utils.taxonomy import aviso_colores, rango_letra
from ui_comun import (NOMBRE_UI, RANGOS_UI, carpeta_trabajo, columnas_cacheadas, columnas_probables, compartir,
                      cobertura_ids, crear_selector_taxa, etiqueta_opcion, guardar_subida, lista, mostrar_imagen,
                      taxa_de_tabla)

FIGURAS = {"tanglegram": "Tanglegram (both trees facing each other)",
           "vamphyre": "VAMPhyRE tree only",
           "gtdb": "GTDB-Tk tree only",
           "individuales": "Both trees, each in its own image"}
MODOS = {"ninguno": "None", "taxonomia": "By taxonomy (pure clades)",
         "concordantes": "Concordant (identical in both trees)",
         "discordantes": "Discordant (conflicting clades)"}


# ------------------------------------------------------------------ sidebar
st.sidebar.title("Parameters")
st.sidebar.caption(f"tanglegram_utils v{__version__}")

with st.sidebar.expander("1 · Files", expanded=True):
    f_v = st.file_uploader("VAMPhyRE tree (Newick)", key="f_v")
    f_g = st.file_uploader("GTDB-Tk tree (Newick)", key="f_g")
    f_t = st.file_uploader("Taxonomy table (xlsx/csv/tsv)", type=["xlsx", "xls", "csv", "tsv", "txt"], key="f_t",
                           help="Required to colour and filter by taxonomy (e.g. gtdbtk.ar53.summary).")
    with st.expander("Harmonise names (optional)"):
        st.caption("Only needed if the table IDs do not match the tree leaf names.")
        f_o = st.file_uploader("Original GTDB tree (GCF_* names)", key="f_o",
                               help="Infers the ID map from the tree topology.")
        f_m = st.file_uploader("ID map (TSV: table_id ⇥ leaf_name)", key="f_m")
        fasta_dir = st.text_input("FASTA folder (path on this computer)", key="fasta_dir") or None

ruta_t = guardar_subida(f_t) if f_t else None
ruta_v = guardar_subida(f_v) if f_v else None
col_id, col_tax, taxa = "user_genome", "classification", {}
if ruta_t:
    try:
        cols = columnas_cacheadas(ruta_t, f_t.size)
        def_id, def_tax = columnas_probables(ruta_t, f_t.size)
        with st.sidebar.expander("Table columns"):
            st.caption("The GTDB-Tk summary (`gtdbtk.ar53.summary`) works as is: "
                       "**user_genome** = genome ID and **classification** = `d__;p__;…` taxonomy. "
                       "Change these only if your table uses other headers.")
            col_id = st.selectbox("Genome ID column", cols, index=cols.index(def_id))
            col_tax = st.selectbox("Taxonomy column", cols, index=cols.index(def_tax))
        taxa = taxa_de_tabla(ruta_t, f_t.size, col_id, col_tax)
        hay_ayuda = bool(f_o or f_m or fasta_dir)
        if ruta_v and not hay_ayuda:
            cob = cobertura_ids(ruta_v, f_v.size, ruta_t, f_t.size, col_id, col_tax)
            if cob < 0.9:
                st.sidebar.warning(
                    f"Only {cob:.0%} of the VAMPhyRE leaves match the table IDs. "
                    "Without extra help **there will be no colours or group names**. Open \"Harmonise names\" "
                    "and upload the **original GTDB tree** (the `.tree` file with GCF_* names), an ID map or the FASTA folder.")
        elif ruta_v:
            st.sidebar.success("IDs: the \"Harmonise names\" input will be used to link the table to the trees.")
    except Exception as e:  # unreadable table or wrong columns: report it and continue without taxonomy
        st.sidebar.warning(f"Could not read the taxonomy from the table: {e}")

compartir(vamphyre=ruta_v, gtdb=guardar_subida(f_g) if f_g else None, tabla=ruta_t, col_id=col_id, col_tax=col_tax)
selector_taxa, OPC, ETIQ = crear_selector_taxa(taxa)


with st.sidebar.expander("2 · Colour and selection", expanded=True):
    color_rango = st.selectbox("Taxonomic rank for colour", RANGOS_UI, index=4,
                               format_func=lambda r: NOMBRE_UI.get(r, r))
    color_top = int(st.number_input("Taxa with their own colour", 1, 100, 8,
                                    help="The rest are grey (\"other\"). Up to 8: colour-blind-safe palette; "
                                         "9-20: unvalidated extra colours; above 20: generated colours, hard to tell apart."))
    if aviso_colores(color_top):
        (st.error if color_top > MAX_DISTINGUIBLE else st.warning)(aviso_colores(color_top))
    resaltar = selector_taxa("Highlight taxa (the rest is dimmed)", "resaltar",
                             "A taxon outside the top N stays grey: raise \"Taxa with their own colour\".")
    ocultar = selector_taxa("Hide connecting lines of…", "ocultar")
    solo_taxa = selector_taxa("Zoom: only these taxa", "solo_taxa",
                              "Prunes both trees and recomputes everything on the subset (min. 4 genomes).")
    raiz_taxon = st.text_input("Root (display only) on taxon", help="Must be monophyletic in both trees.") or None

with st.sidebar.expander("3 · Clade collapsing"):
    modo = st.radio("Collapse", list(MODOS), format_func=MODOS.get)
    umbral, rango_col = None, None
    if modo in ("concordantes", "discordantes"):
        defecto = 1.0 if modo == "concordantes" else 0.5
        u = st.slider("Clade compatibility threshold", 0.0, 1.0, defecto, 0.05,
                      help="concordant: compatibility ≥ threshold · discordant: < threshold")
        umbral = None if u == defecto else u
    if modo == "taxonomia":
        r = st.selectbox("Rank to collapse", ["(same as colour)"] + RANGOS_UI,
                         format_func=lambda r: NOMBRE_UI.get(r, r))
        rango_col = None if r.startswith("(") else r
    min_col = st.number_input("Minimum clade size (leaves)", 2, 1000, 3)
    max_col = st.number_input("Maximum size (0 = no limit)", 0, 5000, 0) or None
    col_nodos = st.text_input("Always collapse these nodes", placeholder="V12,G40",
                              help="IDs from the node CSVs or shown with \"Show IDs of large nodes\".") or None
    exp_nodos = st.text_input("Keep these nodes expanded", placeholder="V7") or None
    col_taxa = selector_taxa("Always collapse these taxa", "col_taxa")
    exp_taxa = selector_taxa("Keep these taxa expanded", "exp_taxa")

with st.sidebar.expander("4 · Appearance"):
    figura = st.selectbox("Figure", list(FIGURAS), format_func=FIGURAS.get,
                          help="A single tree uses the same colour, collapse and label options. "
                               "Concordant/discordant collapsing is computed against the other tree even if it is not drawn.")
    forma = "rectangular"
    if figura != "tanglegram":
        forma = st.radio("Tree layout", ["rectangular", "circular"], horizontal=True,
                         help="The circular layout only applies when single trees are drawn.")
    arco = st.slider("Circle arc (°)", 180, 360, 360, 10,
                     help="The gap is at the top.") if forma == "circular" else 360.0
    soporte, soporte_min = "ninguno", 70.0
    if figura != "tanglegram":
        soporte = st.radio("Node support values", ["ninguno", "numero", "punto"], horizontal=True,
                           format_func=etiqueta_opcion,
                           help="Only if the tree has support values (e.g. GTDB-Tk). Support ≥ the minimum is shown.")
        if soporte != "ninguno":
            soporte_min = st.slider("Minimum support to show (0-100)", 0, 100, 70, 5)
    escala = st.radio("Branch lengths", ["cladograma", "longitud"], horizontal=True, format_func=etiqueta_opcion,
                      help="cladogram: equal branches · branch lengths: real lengths")
    etiquetas_hojas = st.checkbox("Write the name of each leaf", False,
                                  help="Unreadable with hundreds of genomes: use with \"Zoom: only these taxa\" (≲150 genomes).")
    texto_hoja = st.radio("Leaf text", ["especie", "nombre"], horizontal=True, format_func=etiqueta_opcion)
    etiquetas_grupos = st.checkbox("Label pure clades", True)
    min_etiqueta = st.number_input("Minimum size to label", 2, 500, 4)
    mostrar_ids = st.checkbox("Show IDs of large nodes", False)
    min_id = st.number_input("Minimum size to show an ID", 2, 500, 15)
    c1, c2 = st.columns(2)
    ancho_arbol = c1.number_input("Width of each tree (in, 0 = auto)", 0.0, 30.0, 0.0, 0.5,
                                  help="Automatic: 4.5 in the tanglegram, 8 when a single tree is drawn.") or None
    ancho_centro = c2.number_input("Connecting-line area width (in)", 1.0, 20.0, 5.0, 0.5)
    alto = st.number_input("Height (in, 0 = automatic)", 0.0, 200.0, 0.0, 1.0) or None
    dpi = st.slider("Resolution (dpi)", 50, 300, 110, 10, help="110 to explore; 300 for publication.")
    grosor_linea = st.slider("Connecting-line width", 0.3, 4.0, 1.3, 0.1)
    grosor_rama = st.slider("Branch width", 0.3, 3.0, 1.1, 0.1)
    esc_tri = st.slider("Triangle size", 0.5, 2.0, 1.0, 0.1)
    sin_desenredar = st.checkbox("Do not reorder branches (no untangling)", False)
    formatos = st.multiselect("Additional formats", ["svg", "pdf"])

with st.sidebar.expander("5 · Rogue analysis"):
    max_reb = st.number_input("Greedy removal steps (0 = skip)", 0, 200, 25)
    rangos_grupos = st.multiselect("Ranks for the group table", RANGOS_UI, default=["familia", "genero"],
                                   format_func=lambda r: NOMBRE_UI.get(r, r))
    min_grupo = st.number_input("Minimum group size", 2, 500, 3)

# ---------------------------------------------------------------------- main area
st.title("Tanglegram: VAMPhyRE vs GTDB-Tk")
st.caption("Topology comparison, rogue branches, and taxonomy-based collapsing and colouring. "
           "Everything is processed on your computer.")

listo = bool(f_v and f_g)
if not listo:
    st.info("**Start here:** drop the two Newick trees (and the taxonomy table if you want taxonomy colours) "
            "into the left panel, adjust the settings and press **Generate**.")

if st.button("Generate tanglegram", type="primary", disabled=not listo):
    entradas = Entradas(tree_vamphyre=ruta_v, tree_gtdb=guardar_subida(f_g), tabla=ruta_t,
                        col_id=col_id, col_tax=col_tax,
                        tree_gtdb_original=guardar_subida(f_o) if f_o else None,
                        fasta_dir=fasta_dir, mapa_ids=guardar_subida(f_m) if f_m else None)
    op = OpcionesAnalisis(color_rango=color_rango, color_top=color_top, resaltar=resaltar, ocultar_lineas=ocultar,
                          solo_taxa=solo_taxa, raiz_taxon=raiz_taxon, sin_desenredar=sin_desenredar,
                          max_rebeldes=int(max_reb), rango_grupos=",".join(rangos_grupos) or "familia,genero",
                          min_grupo=int(min_grupo))
    colapso = OpcionesColapso(modo=modo, umbral=umbral, min_colapso=int(min_col),
                              max_colapso=int(max_col) if max_col else None,
                              rango=rango_letra(rango_col or color_rango),
                              colapsar_nodos=col_nodos, expandir_nodos=exp_nodos,
                              colapsar_taxa=col_taxa, expandir_taxa=exp_taxa)
    estilo = EstiloFigura(escala=escala, etiquetas_hojas=etiquetas_hojas, texto_hoja=texto_hoja,
                          etiquetas_grupos=etiquetas_grupos, min_etiqueta=int(min_etiqueta),
                          mostrar_ids=mostrar_ids, min_id=int(min_id), ancho_arbol=ancho_arbol,
                          ancho_centro=ancho_centro, alto=alto, dpi=int(dpi), grosor_linea=grosor_linea,
                          grosor_rama=grosor_rama, escala_triangulo=esc_tri, forma=forma, arco=float(arco),
                          soporte=soporte, soporte_min=float(soporte_min))
    corrida = os.path.join(carpeta_trabajo(), "salida")
    shutil.rmtree(corrida, ignore_errors=True)
    st.session_state.pop("error", None)
    with st.spinner("Computing… (with ~500 genomes this takes from a few seconds to a minute)"):
        try:
            res_nuevo = ejecutar_en_directorio(entradas, op, colapso, estilo, corrida,
                                               formato_extra=lista(formatos), figura=figura)
            st.session_state["res"] = res_nuevo
            st.session_state["zip"] = empaquetar_zip(res_nuevo.archivos)
        except ErrorDatos as e:
            st.session_state.pop("res", None)
            st.session_state["error"] = (str(e), getattr(e, "log", ""))
        except Exception as e:  # unexpected error: show it without crashing the page
            st.session_state.pop("res", None)
            st.session_state["error"] = (f"{type(e).__name__}: {e}", "")

if "error" in st.session_state:
    msg, log = st.session_state["error"]
    st.error(msg)
    if log:
        with st.expander("Log"):
            st.code(log)

res = st.session_state.get("res")
if res:
    m = st.columns([1, 1.2, 1, 1.2, 2])
    m[0].metric("Genomes", res.n_genomas)
    m[1].metric("With taxonomy", f"{res.n_con_taxonomia}/{res.n_genomas}",
                help="Genomes linked to the table. If 0, there will be no colours or names: check \"Harmonise names\".")
    m[2].metric("RF distance", res.rf, help="Robinson-Foulds: splits not shared by the two trees.")
    m[3].metric("Normalised RF", f"{res.rf_norm:.3f}", help="0 = identical · 1 = no split in common")
    m[4].metric("Line crossings", f"{res.cruces_antes:,} → {res.cruces_despues:,}",
                help="Before (after ladderize) and after reordering branches.")
    avisos = [l.strip().replace("[!] ", "") for l in res.log.splitlines() if "[!]" in l]
    for a in avisos[:3]:
        st.warning(a)
    if len(avisos) > 3:
        with st.expander(f"{len(avisos) - 3} more warning(s)"):
            for a in avisos[3:]:
                st.write(a)

    t_fig, t_reb, t_rep, t_tab, t_des, t_log = st.tabs(
        ["Tanglegram", "Rogues", "Report", "Tables", "Downloads", "Log"])
    with t_fig:
        imgs = res.imagenes
        if len(imgs) > 1:
            ruta_fig = st.radio("Image", imgs, horizontal=True,
                                format_func=lambda r: os.path.basename(r)[len("tanglegram_"):-4])
        else:
            ruta_fig = imgs[0]
        mostrar_imagen(ruta_fig, "fig")
    with t_reb:
        reb = res.con_sufijo("_rebeldes.png")
        if reb:
            mostrar_imagen(reb[0], "reb")
        else:
            st.info("The rogue summary figure was not generated.")
    with t_rep:
        rep = res.con_sufijo("_reporte.md")
        if rep:
            with open(rep[0], encoding="utf-8") as f:
                st.markdown(f.read())
    with t_tab:
        csvs = res.con_sufijo(".csv")
        if csvs:
            elegido = st.selectbox("Table", csvs, format_func=lambda r: os.path.basename(r)[len("tanglegram_"):])
            st.dataframe(pd.read_csv(elegido), width="stretch", hide_index=True)
    with t_des:
        st.download_button("⬇ Everything as a ZIP", st.session_state["zip"], "tanglegram_results.zip",
                           "application/zip", type="primary", on_click="ignore")
        st.markdown("**Individual files**")
        for i, ruta in enumerate(res.archivos):
            c_a, c_b = st.columns([4, 1])
            c_a.write(os.path.basename(ruta))
            with open(ruta, "rb") as f:
                c_b.download_button("Download", f.read(), os.path.basename(ruta), key=f"dl_{i}", on_click="ignore")
        st.markdown("**Equivalent command** (to reproduce this figure from the terminal)")
        st.code(res.comando, language="bash")
        st.caption("Replace the paths with your own files; `--prefijo-salida` or `--out-plot` set where outputs are saved.")
    with t_log:
        st.code(res.log or "(no messages)")
