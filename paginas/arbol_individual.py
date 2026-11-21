"""App page "Single tree": draws ONE tree (VAMPhyRE, GTDB-Tk or a future consensus),
rectangular or circular, with taxonomy colours/collapsing, support values and iTOL export.

Interface only: the logic lives in `tanglegram_utils` (it calls `servicio.ejecutar_arbol_en_directorio`,
just as `generate_tree.py` calls `pipeline_arbol`)."""
import os
import shutil

import pandas as pd
import streamlit as st

from tanglegram_utils import __version__
from tanglegram_utils.collapse import OpcionesColapso
from tanglegram_utils.errors import ErrorDatos
from tanglegram_utils.models import EntradaArbol, OpcionesArbol
from tanglegram_utils.plotting import EstiloFigura
from tanglegram_utils.servicio import ejecutar_arbol_en_directorio, empaquetar_zip
from tanglegram_utils.constants import MAX_DISTINGUIBLE
from tanglegram_utils.taxonomy import aviso_colores, rango_letra
from ui_comun import (NOMBRE_UI, RANGOS_UI, carpeta_trabajo, columnas_cacheadas, columnas_probables,
                      compartidos, crear_selector_taxa, etiqueta_opcion, guardar_subida, info_arbol, lista,
                      mostrar_imagen, taxa_de_tabla)

MODOS = {"ninguno": "None", "taxonomia": "By taxonomy (pure clades)"}
RAICES = {"original": "As in the file", "punto-medio": "Midpoint of the tree",
          "taxon": "Outgroup (a taxon)"}
SUBIR = "Upload a file"

comp = compartidos()  # files already loaded on the "Tanglegram" page
res_tg = st.session_state.get("res")
mapa_tg = (res_tg.con_sufijo("_mapa_ids.tsv") or [None])[0] if res_tg is not None and hasattr(res_tg, "con_sufijo") else None

st.sidebar.title("Parameters")
st.sidebar.caption(f"tanglegram_utils v{__version__}")

# ------------------------------------------------------------------ 1 · files
ruta_a = ruta_t = None
with st.sidebar.expander("1 · Files", expanded=True):
    opciones_arbol = [SUBIR] + [{"vamphyre": "VAMPhyRE (loaded on Tanglegram page)",
                                 "gtdb": "GTDB-Tk (loaded on Tanglegram page)"}[k]
                                for k in ("vamphyre", "gtdb") if k in comp]
    origen = st.radio("Tree", opciones_arbol, key="ai_origen") if len(opciones_arbol) > 1 else opciones_arbol[0]
    if origen == SUBIR:
        f_a = st.file_uploader("Tree (Newick)", key="ai_arbol",
                               help="VAMPhyRE, GTDB-Tk or a consensus tree with bootstrap on internal nodes.")
        ruta_a = guardar_subida(f_a) if f_a else None
        tam_a = f_a.size if f_a else 0
    else:
        ruta_a = comp["vamphyre" if origen.startswith("VAMPhyRE") else "gtdb"]
        tam_a = os.path.getsize(ruta_a)
        st.caption(f"Using `{os.path.basename(ruta_a)}`.")

    if "tabla" in comp:
        usar_tabla = st.checkbox("Use the table loaded on the Tanglegram page", True, key="ai_usar_tabla")
    else:
        usar_tabla = False
    if usar_tabla:
        ruta_t, tam_t = comp["tabla"], os.path.getsize(comp["tabla"])
        st.caption(f"Table: `{os.path.basename(ruta_t)}`.")
    else:
        f_t = st.file_uploader("Taxonomy table (xlsx/csv/tsv)", type=["xlsx", "xls", "csv", "tsv", "txt"],
                               key="ai_tabla", help="For taxonomy colours and labels (e.g. gtdbtk.ar53.summary).")
        ruta_t = guardar_subida(f_t) if f_t else None
        tam_t = f_t.size if f_t else 0

    with st.expander("Harmonise names (optional)"):
        st.caption("Only needed if the table IDs do not match the tree leaf names.")
        usar_mapa_tg = bool(mapa_tg) and st.checkbox(
            "Use the ID map from the last tanglegram run", True, key="ai_usar_mapa_tg")
        f_m = st.file_uploader("ID map (TSV: table_id ⇥ leaf_name)", key="ai_mapa")
        f_o = st.file_uploader("Reference tree (same names as the table, same topology)", key="ai_ref",
                               help="Usually the original GTDB tree (GCF_* names).")
        fasta_dir = st.text_input("FASTA folder (path on this computer)", key="ai_fasta") or None

col_id, col_tax, taxa = comp.get("col_id", "user_genome"), comp.get("col_tax", "classification"), {}
if ruta_t:
    try:
        cols = columnas_cacheadas(ruta_t, tam_t)
        def_id, def_tax = columnas_probables(ruta_t, tam_t)
        with st.sidebar.expander("Table columns"):
            col_id = st.selectbox("Genome ID column", cols,
                                  index=cols.index(col_id if col_id in cols else def_id), key="ai_col_id")
            col_tax = st.selectbox("Taxonomy column", cols,
                                   index=cols.index(col_tax if col_tax in cols else def_tax), key="ai_col_tax")
        taxa = taxa_de_tabla(ruta_t, tam_t, col_id, col_tax)
    except Exception as e:  # unreadable table: report it and continue without taxonomy
        st.sidebar.warning(f"Could not read the taxonomy from the table: {e}")

selector_taxa, OPC, ETIQ = crear_selector_taxa(taxa)

# What the tree contains (leaves and support): shown when it is loaded
info = None
if ruta_a:
    try:
        info = info_arbol(ruta_a, os.path.getsize(ruta_a))
    except Exception as e:
        st.sidebar.error(f"Could not read the tree: {e}")

# ------------------------------------------------------------------ 2 · root and support
with st.sidebar.expander("2 · Root and support", expanded=True):
    raiz = st.radio("Root", list(RAICES), format_func=RAICES.get, key="ai_raiz",
                    help="Changes only the drawing and the number of fragments per group (not the unrooted topology).")
    raiz_taxon = None
    if raiz == "taxon":
        raiz_taxon = st.text_input("Outgroup (taxon)", placeholder="p__Thermoproteota",
                                   help="Must form a clade in the tree.") or None
    n_soporte = (info[1] or {}).get("con_soporte", 0) if info else 0
    contraer = None
    if info and not n_soporte:
        st.caption("This tree has no support values (names such as `Inner123` do not count).")
    elif info:
        s = info[1]
        st.caption(f"Support (0-100): median {s['mediana']:.0f}, minimum {s['minimo']:.0f} · "
                   f"{s['pct_ge_70']:.0f}% of nodes ≥ 70 · {s['pct_ge_90']:.0f}% ≥ 90.")
        if st.checkbox("Collapse low-support nodes (polytomies)", False, key="ai_contraer_on"):
            contraer = float(st.slider("Collapse nodes with support below…", 0, 100, 70, 5, key="ai_contraer"))

# ------------------------------------------------------------------ 3 · colour and selection
with st.sidebar.expander("3 · Colour and selection", expanded=True):
    color_rango = st.selectbox("Taxonomic rank for colour", RANGOS_UI, index=4, key="ai_color_rango",
                               format_func=lambda r: NOMBRE_UI.get(r, r))
    color_top = int(st.number_input("Taxa with their own colour", 1, 100, 8, key="ai_color_top",
                                    help="The rest are grey (\"other\"). Up to 8: colour-blind-safe palette; "
                                         "9-20: unvalidated extra colours; above 20: generated colours, hard to tell apart."))
    if aviso_colores(color_top):
        (st.error if color_top > MAX_DISTINGUIBLE else st.warning)(aviso_colores(color_top))
    resaltar = selector_taxa("Highlight taxa (the rest is dimmed)", "ai_resaltar")
    solo_taxa = selector_taxa("Zoom: only these taxa", "ai_solo_taxa",
                              "Prunes the tree and recomputes on the subset (min. 4 genomes).")

# ------------------------------------------------------------------ 4 · collapsing
with st.sidebar.expander("4 · Clade collapsing"):
    modo = st.radio("Collapse", list(MODOS), format_func=MODOS.get, key="ai_modo")
    rango_col = None
    if modo == "taxonomia":
        r = st.selectbox("Rank to collapse", ["(same as colour)"] + RANGOS_UI, key="ai_rango_col",
                         format_func=lambda r: NOMBRE_UI.get(r, r))
        rango_col = None if r.startswith("(") else r
    min_col = st.number_input("Minimum clade size (leaves)", 2, 1000, 3, key="ai_min_col")
    max_col = st.number_input("Maximum size (0 = no limit)", 0, 5000, 0, key="ai_max_col") or None
    col_nodos = st.text_input("Always collapse these nodes", placeholder="N12,N40", key="ai_col_nodos",
                              help="IDs from the node CSV or shown with \"Show IDs of large nodes\".") or None
    exp_nodos = st.text_input("Keep these nodes expanded", placeholder="N7", key="ai_exp_nodos") or None
    col_taxa = selector_taxa("Always collapse these taxa", "ai_col_taxa")
    exp_taxa = selector_taxa("Keep these taxa expanded", "ai_exp_taxa")

# ------------------------------------------------------------------ 5 · appearance
with st.sidebar.expander("5 · Appearance", expanded=True):
    forma = st.radio("Layout", ["circular", "rectangular"], horizontal=True, key="ai_forma")
    arco = st.slider("Circle arc (°)", 180, 360, 360, 10, key="ai_arco",
                     help="The gap is at the top.") if forma == "circular" else 360
    soporte = st.radio("Node support values", ["ninguno", "numero", "punto"], horizontal=True,
                       key="ai_soporte", format_func=etiqueta_opcion,
                       help="Only if the tree has support values. Support ≥ the minimum is shown.")
    soporte_min = st.slider("Minimum support to show (0-100)", 0, 100, 70, 5, key="ai_soporte_min") \
        if soporte != "ninguno" else 70
    escala = st.radio("Branch lengths", ["cladograma", "longitud"], horizontal=True, key="ai_escala",
                      format_func=etiqueta_opcion, help="cladogram: equal branches · branch lengths: real lengths")
    orden = st.radio("Branch order", ["ladderize", "original"], horizontal=True, key="ai_orden",
                     format_func=etiqueta_opcion)
    etiquetas_hojas = st.checkbox("Write the name of each leaf", False, key="ai_etq_hojas",
                                  help="Unreadable with hundreds of genomes: use with \"Zoom: only these taxa\".")
    texto_hoja = st.radio("Leaf text", ["especie", "nombre"], horizontal=True, key="ai_texto_hoja",
                          format_func=etiqueta_opcion)
    etiquetas_grupos = st.checkbox("Label pure clades", True, key="ai_etq_grupos")
    min_etiqueta = st.number_input("Minimum size to label", 2, 500, 4, key="ai_min_etq")
    mostrar_ids = st.checkbox("Show IDs of large nodes", False, key="ai_ids")
    min_id = st.number_input("Minimum size to show an ID", 2, 500, 15, key="ai_min_id")
    ancho = st.number_input("Tree radius / width (in, 0 = auto)", 0.0, 30.0, 0.0, 0.5, key="ai_ancho") or None
    alto = st.number_input("Height (in, 0 = automatic)", 0.0, 200.0, 0.0, 1.0, key="ai_alto") or None
    dpi = st.slider("Resolution (dpi)", 50, 300, 110, 10, key="ai_dpi", help="110 to explore; 300 for publication.")
    grosor_rama = st.slider("Branch width", 0.3, 3.0, 1.1, 0.1, key="ai_grosor")
    esc_tri = st.slider("Triangle size", 0.5, 2.0, 1.0, 0.1, key="ai_esc_tri")
    formatos = st.multiselect("Additional formats", ["svg", "pdf"], key="ai_formatos")
    itol = st.checkbox("Export for iTOL (tree + editable annotations)", True, key="ai_itol",
                       help="Writes _itol.tree, .txt annotation files and _itol.zip to upload to itol.embl.de.")

with st.sidebar.expander("6 · Group tables"):
    rangos_grupos = st.multiselect("Ranks for the monophyly table", RANGOS_UI, default=["familia", "genero"],
                                   format_func=lambda r: NOMBRE_UI.get(r, r), key="ai_rangos_grupos")
    min_grupo = st.number_input("Minimum group size", 2, 500, 3, key="ai_min_grupo")
titulo = st.sidebar.text_input("Figure title (optional)", key="ai_titulo") or None

# ---------------------------------------------------------------------- main area
st.title("Single tree")
st.caption("Draw one tree (VAMPhyRE, GTDB-Tk or a consensus) with colours, collapsing and support, "
           "in rectangular or circular layout. Everything is processed on your computer.")

if not ruta_a:
    st.info("**Start here:** upload a Newick tree in the left panel (or load both trees on the "
            "\"Tanglegram\" page and select one here). Add the taxonomy table for taxonomy colours.")
elif info:
    st.success(f"{info[0]} leaves" + (f" · support on {n_soporte} internal nodes" if n_soporte else " · no support values"))

if st.button("Generate tree", type="primary", disabled=not ruta_a, key="ai_generar"):
    ruta_mapa = guardar_subida(f_m) if f_m else (mapa_tg if usar_mapa_tg else None)
    entrada = EntradaArbol(arbol=ruta_a, tabla=ruta_t, col_id=col_id, col_tax=col_tax, mapa_ids=ruta_mapa,
                           fasta_dir=fasta_dir, arbol_referencia=guardar_subida(f_o) if f_o else None,
                           titulo=titulo)
    op = OpcionesArbol(color_rango=color_rango, color_top=color_top, resaltar=resaltar, solo_taxa=solo_taxa,
                       raiz=raiz, raiz_taxon=raiz_taxon, contraer_soporte=contraer, orden=orden,
                       rango_grupos=",".join(rangos_grupos) or "familia,genero", min_grupo=int(min_grupo),
                       itol=itol)
    colapso = OpcionesColapso(modo=modo, min_colapso=int(min_col), max_colapso=int(max_col) if max_col else None,
                              rango=rango_letra(rango_col or color_rango),
                              colapsar_nodos=col_nodos, expandir_nodos=exp_nodos,
                              colapsar_taxa=col_taxa, expandir_taxa=exp_taxa)
    estilo = EstiloFigura(escala=escala, etiquetas_hojas=etiquetas_hojas, texto_hoja=texto_hoja,
                          etiquetas_grupos=etiquetas_grupos, min_etiqueta=int(min_etiqueta),
                          mostrar_ids=mostrar_ids, min_id=int(min_id), ancho_arbol=ancho, alto=alto,
                          dpi=int(dpi), grosor_rama=grosor_rama, escala_triangulo=esc_tri, forma=forma,
                          arco=float(arco), soporte=soporte, soporte_min=float(soporte_min))
    corrida = os.path.join(carpeta_trabajo(), "salida_arbol")
    shutil.rmtree(corrida, ignore_errors=True)
    st.session_state.pop("ai_error", None)
    with st.spinner("Drawing…"):
        try:
            r = ejecutar_arbol_en_directorio(entrada, op, colapso, estilo, corrida, formato_extra=lista(formatos))
            st.session_state["ai_res"] = r
            st.session_state["ai_zip"] = empaquetar_zip(r.archivos)
        except ErrorDatos as e:
            st.session_state.pop("ai_res", None)
            st.session_state["ai_error"] = (str(e), getattr(e, "log", ""))
        except Exception as e:  # unexpected error: show it without crashing the page
            st.session_state.pop("ai_res", None)
            st.session_state["ai_error"] = (f"{type(e).__name__}: {e}", "")

if "ai_error" in st.session_state:
    msg, log = st.session_state["ai_error"]
    st.error(msg)
    if log:
        with st.expander("Log"):
            st.code(log)

r = st.session_state.get("ai_res")
if r:
    m = st.columns([1, 1.2, 1.4, 2])
    m[0].metric("Genomes", r.n_genomas)
    m[1].metric("With taxonomy", f"{r.n_con_taxonomia}/{r.n_genomas}",
                help="Genomes linked to the table. If 0, there will be no colours or names: check \"Harmonise names\".")
    m[2].metric("Collapsed nodes", r.n_contraidos, help="Turned into polytomies due to low support.")
    m[3].metric("Root", r.motivo_raiz or "—")
    avisos = [l.strip().replace("[!] ", "") for l in r.log.splitlines() if "[!]" in l]
    for a in avisos[:3]:
        st.warning(a)
    if len(avisos) > 3:
        with st.expander(f"{len(avisos) - 3} more warning(s)"):
            for a in avisos[3:]:
                st.write(a)

    t_fig, t_rep, t_tab, t_des, t_log = st.tabs(["Tree", "Report", "Tables", "Downloads", "Log"])
    with t_fig:
        mostrar_imagen(r.imagen, "ai_fig")
    with t_rep:
        rep = r.con_sufijo("_reporte.md")
        if rep:
            with open(rep[0], encoding="utf-8") as f:
                st.markdown(f.read())
        else:
            st.info("No report (it is generated with a taxonomy table).")
    with t_tab:
        csvs = r.con_sufijo(".csv")
        if csvs:
            elegido = st.selectbox("Table", csvs, format_func=lambda x: os.path.basename(x)[len("arbol_"):],
                                   key="ai_tabla_sel")
            st.dataframe(pd.read_csv(elegido), width="stretch", hide_index=True)
        else:
            st.info("No tables were generated.")
    with t_des:
        st.download_button("⬇ Everything as a ZIP", st.session_state["ai_zip"], "tree_results.zip",
                           "application/zip", type="primary", on_click="ignore", key="ai_dl_zip")
        itol_zip = r.con_sufijo("_itol.zip")
        if itol_zip:
            with open(itol_zip[0], "rb") as f:
                st.download_button("⬇ iTOL package", f.read(), os.path.basename(itol_zip[0]), "application/zip",
                                   on_click="ignore", key="ai_dl_itol")
            st.caption("On itol.embl.de: *Upload* → *Batch upload* with this ZIP, or upload `arbol_itol.tree` and "
                       "drag the `.txt` files onto the tree. Each layer stays editable in iTOL.")
        st.markdown("**Individual files**")
        for i, ruta in enumerate(r.archivos):
            c_a, c_b = st.columns([4, 1])
            c_a.write(os.path.basename(ruta))
            with open(ruta, "rb") as f:
                c_b.download_button("Download", f.read(), os.path.basename(ruta), key=f"ai_dl_{i}",
                                    on_click="ignore")
        st.markdown("**Equivalent command** (to reproduce this figure from the terminal)")
        st.code(r.comando, language="bash")
        st.caption("Replace the paths with your own files; `--out-plot` sets where outputs are saved.")
    with t_log:
        st.code(r.log or "(no messages)")
