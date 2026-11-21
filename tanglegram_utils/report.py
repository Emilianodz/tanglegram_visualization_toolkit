"""CSV tables, summary figure and Markdown report of the rogue-branch analysis."""
import collections
import logging
import os

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .bitutils import bits
from .constants import AZUL_SERIE, LETRA_A_NOMBRE, NARANJA_SERIE, TINTA
from .splits import canon, fragmentos
from .taxonomy import etiqueta_consenso, rango_letra, sin_prefijo
from .tree import preorden

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------- tables
def tabla_hojas(datos, an, max_rebeldes):
    """One row per genome: RF contribution, removal order and ancestor concordance."""
    tax, nombres = datos.tax, datos.nombres
    paso_de = {i: k + 1 for k, (i, _, _) in enumerate(an.pasos)}
    rf_tras = {i: r2 for i, _, r2 in an.pasos}
    filas = []
    for i in sorted(an.conservar):
        cV, tV = an.ancV.get(i, (0, 0))
        cG, tG = an.ancG.get(i, (0, 0))
        filas.append({
            "genoma": nombres[i], "especie": sin_prefijo(tax[i].get("s", "")),
            "genero": sin_prefijo(tax[i].get("g", "")), "familia": sin_prefijo(tax[i].get("f", "")),
            "orden": sin_prefijo(tax[i].get("o", "")), "filo": sin_prefijo(tax[i].get("p", "")),
            "delta_RF_al_quitarlo": an.loo.get(i, ""), "orden_rebelde_voraz": paso_de.get(i, ""),
            "RF_tras_quitarlo": rf_tras.get(i, ""),
            "ancestros_concordantes_V": f"{cV}/{tV}", "ancestros_concordantes_G": f"{cG}/{tG}",
            "frac_concordante_V": round(cV / tV, 3) if tV else "", "frac_concordante_G": round(cG / tG, 3) if tG else ""})
    df = pd.DataFrame(filas)
    if max_rebeldes > 0:
        df = df.assign(_o=pd.to_numeric(df["orden_rebelde_voraz"], errors="coerce").fillna(1e9),
                       _d=pd.to_numeric(df["delta_RF_al_quitarlo"], errors="coerce").fillna(-1)
                       ).sort_values(["_o", "_d"], ascending=[True, False]).drop(columns=["_o", "_d"])
    return df


def tabla_nodos(datos, an, raiz, info, ids_otro):
    """One row per internal branch of `raiz`: concordance and best match in the other tree."""
    ids_otra = {n: n.id for n in ids_otro.values()}
    hay_tax = datos.esquema.hay_tax
    filas = []
    for n in preorden(raiz):
        if n.es_hoja or n.padre is None:
            continue
        d = info[n]
        lab, frac = etiqueta_consenso(n.hojas, datos.tax) if hay_tax else ("", 0)
        sobran = faltan = ""
        if not d["conc"] and d["pareja"] is not None:
            sobran = "; ".join(datos.corto(i) for i in bits(n.mask & ~d["pareja"].mask & an.todo))
            faltan = "; ".join(datos.corto(i) for i in bits(d["pareja"].mask & ~n.mask & an.todo))
        filas.append({"id": n.id, "n_hojas": n.n, "concordante": bool(d["conc"]),
                      "compatibilidad_subarbol": round(d["comp"], 3), "jaccard_mejor_pareja": round(d["jmax"], 3),
                      "hojas_de_diferencia": d["dif"], "pareja_en_otro_arbol": ids_otra.get(d["pareja"], "") if d["pareja"] else "",
                      "taxon_consenso": sin_prefijo(lab), "pureza_consenso": round(frac, 2),
                      "sobran_vs_pareja": sobran[:300], "faltan_vs_pareja": faltan[:300]})
    return pd.DataFrame(filas).sort_values("n_hojas", ascending=False)


def tabla_grupos(datos, an, rangos, min_grupo):
    """Per taxonomic group metrics: monophyly, fragments, internal concordance."""
    tax = datos.tax
    top_rebeldes = {i for i, _, _ in an.pasos}
    filas = []
    for rl in rangos:
        miembros = collections.defaultdict(list)
        for i in an.conservar:
            c = tax[i].get(rl)
            if c:
                miembros[c].append(i)
        for c, lst in miembros.items():
            if len(lst) < min_grupo:
                continue
            S = 0
            for i in lst:
                S |= 1 << i
            fila = {"rango": LETRA_A_NOMBRE[rl], "grupo": sin_prefijo(c), "n": len(lst)}
            tot_c = tot_t = 0
            for key, raiz, info, spl in (("V", an.rV, an.infoV, an.sV), ("G", an.rG, an.infoG, an.sG)):
                cs = canon(S, an.todo)
                mono = True if cs is None else (cs in spl)
                mrca = min((n for n in preorden(raiz) if not n.es_hoja and (n.mask & S) == S), key=lambda n: n.n)
                internos = [n for n in preorden(raiz) if (not n.es_hoja) and n.padre is not None and (n.mask & ~S) == 0]
                cc = sum(1 for n in internos if info[n]["conc"])
                fila[f"monofiletico_{key}"] = mono
                fila[f"pureza_{key}"] = round(len(lst) / mrca.n, 3)
                fila[f"intrusos_{key}"] = mrca.n - len(lst)
                fila[f"fragmentos_{key}"] = fragmentos(raiz, S)
                fila[f"nodos_internos_{key}"] = len(internos)
                fila[f"concordantes_{key}"] = cc
                tot_c += cc
                tot_t += len(internos)
            fila["concordancia_interna"] = round(tot_c / tot_t, 3) if tot_t else np.nan
            fila["rebeldes_en_grupo"] = sum(1 for i in lst if i in top_rebeldes)
            filas.append(fila)
    return pd.DataFrame(filas)


# ---------------------------------------------------------------------- summary figure
def figura_resumen(datos, an, dfg, rangos, ruta):
    """Left panel: rogue genomes (ΔRF per step). Right panel: most fragmented groups."""
    pasos = an.pasos
    fig, axs = plt.subplots(1, 2, figsize=(16, max(5.5, 0.34 * max(len(pasos), 12) + 2)), facecolor="white")
    a = axs[0]
    if pasos:
        etq = [f"{k + 1}. {datos.corto(i)}" for k, (i, _, _) in enumerate(pasos)]
        vals = [r1 - r2 for _, r1, r2 in pasos]
        a.barh(range(len(pasos)), vals, color=AZUL_SERIE, height=0.62)
        a.set_yticks(range(len(pasos)))
        a.set_yticklabels(etq, fontsize=9, color=TINTA)
        a.invert_yaxis()
        for k, v in enumerate(vals):
            a.text(v + 0.1, k, str(v), va="center", fontsize=9, color=TINTA)
        a.set_xlabel("RF reduction when removing the genome (step by step)", color=TINTA)
        a.set_title(f"Rogue genomes · RF {an.rf} → {pasos[-1][2]} after {len(pasos)} genomes", loc="left", color=TINTA, fontsize=12)
    else:
        a.axis("off")
    a = axs[1]
    if not dfg.empty:
        ult = LETRA_A_NOMBRE[rangos[-1]]
        sub_ = dfg[(dfg["rango"] == ult)]
        peor = sub_.sort_values(["fragmentos_V", "n"], ascending=[False, False]).head(15)
        yy = np.arange(len(peor))
        a.barh(yy - 0.19, peor["fragmentos_V"], color=AZUL_SERIE, height=0.36, label="VAMPhyRE")
        a.barh(yy + 0.19, peor["fragmentos_G"], color=NARANJA_SERIE, height=0.36, label="GTDB-Tk")
        a.set_yticks(yy)
        a.set_yticklabels([f"{g} (n={n})" for g, n in zip(peor["grupo"], peor["n"])], fontsize=9, color=TINTA)
        a.invert_yaxis()
        for k, (fv, fg) in enumerate(zip(peor["fragmentos_V"], peor["fragmentos_G"])):
            a.text(fv + 0.15, k - 0.19, str(fv), va="center", fontsize=8.5, color=TINTA)
            a.text(fg + 0.15, k + 0.19, str(fg), va="center", fontsize=8.5, color=TINTA)
        a.set_xlabel("Fragments the group is split into (1 = monophyletic)", color=TINTA)
        a.legend(frameon=False, loc="lower right", fontsize=10)
        a.set_title(f"Most fragmented groups ({ult}) in VAMPhyRE", loc="left", color=TINTA, fontsize=12)
    else:
        a.axis("off")
    for a in axs:
        for sp_name in ("top", "right"):
            a.spines[sp_name].set_visible(False)
        a.tick_params(colors=TINTA)
    fig.tight_layout()
    fig.savefig(ruta, dpi=140, facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------------- report
def _bloques_identicos(datos, an):
    """VAMPhyRE subtrees identical in both trees (outermost only)."""
    ident = []
    for n in preorden(an.rV):
        if n.es_hoja or n.padre is None:
            continue
        pa = n.padre
        if an.infoV[n]["comp"] >= 1 - 1e-12 and n.n >= 3 and (pa.padre is None or an.infoV[pa]["comp"] < 1 - 1e-12):
            lab, fr = etiqueta_consenso(n.hojas, datos.tax) if datos.esquema.hay_tax else ("", 0)
            ident.append((n.n, n.id, sin_prefijo(lab), fr))
    ident.sort(reverse=True)
    return ident


def escribir_reporte(ruta, prefijo, datos, an, dfnV, dfnG, dfg, rangos, colV, colG,
                     modo_colapso, min_grupo):
    """Write the Markdown report (summary, rogues, branches, identical blocks, groups)."""
    tax, N, N0 = datos.tax, an.N, datos.N0
    rf, nrf, pasos, loo = an.rf, an.nrf, an.pasos, an.loo
    corto = datos.corto
    L = []
    L.append("# Rogue-branch report: VAMPhyRE vs GTDB-Tk\n")
    L.append("## 1. Summary\n")
    L.append(f"- Genomes compared: **{N}**" + (f" (of {N0} shared; `--solo-taxa` filter)" if N != N0 else ""))
    L.append(f"- Internal branches (splits): VAMPhyRE {len(an.sV)}, GTDB {len(an.sG)}, shared **{len(an.sV & an.sG)}**")
    L.append(f"- **RF = {rf}**, normalised = **{nrf:.3f}** (0 = identical trees, 1 = no branches in common)")
    L.append(f"- VAMPhyRE branches present in GTDB: {100 * len(an.sV & an.sG) / max(len(an.sV), 1):.1f}%")
    L.append(f"- Line crossings in the tanglegram: {an.cr0} → {an.cr1}" + ("" if not an.sin_desenredar else " (not untangled)"))
    L.append(f"- Common display root: {an.motivo_raiz} (affects the drawing only, not RF)")
    if colV or colG:
        L.append(f"- Collapsed clades (`--colapsar {modo_colapso}`): VAMPhyRE {len(colV)}, GTDB {len(colG)}")
    L.append("")
    L.append("## 2. Rogue genomes (greedy removal)\n")
    if pasos:
        L.append("Each row is the genome whose removal lowers RF the most at that step (removals are cumulative).\n")
        L.append("| # | Genome | Family | Genus | ΔRF | Remaining RF |")
        L.append("|---|---|---|---|---|---|")
        for k, (i, r1, r2) in enumerate(pasos):
            L.append(f"| {k + 1} | {corto(i)} | {sin_prefijo(tax[i].get('f', ''))} | {sin_prefijo(tax[i].get('g', ''))} | −{r1 - r2} | {r2} |")
        pct = 100 * (rf - pasos[-1][2]) / rf if rf else 0
        L.append(f"\nRemoving **{len(pasos)} genomes** ({100 * len(pasos) / N:.1f}% of the total) lowers RF from {rf} to {pasos[-1][2]} (**−{pct:.0f}%**).")
        n_cero = sum(1 for v in loo.values() if v > 0)
        L.append(f"{n_cero} genomes lower RF when removed individually; the rest do not change RF on their own.")
        if pct < 30:
            L.append(f"**Interpretation:** the disagreement is diffuse and not caused by a few genomes (each rogue accounts for only ~{(rf - pasos[-1][2]) / len(pasos):.0f} RF units). "
                     "The conflict is structural (arrangement of whole groups); sections 3 and 5 describe it better.")
        L.append("")
    else:
        L.append("(skipped or RF = 0)\n")
    L.append("## 3. Largest discordant branches\n")
    L.append("`hojas_de_diferencia` = genomes that would have to move for the branch to match its closest counterpart in the other tree.\n")
    for nombre_t, dfn in (("VAMPhyRE", dfnV), ("GTDB", dfnG)):
        d_ = dfn[~dfn["concordante"]].sort_values(["hojas_de_diferencia", "n_hojas"], ascending=[False, False]).head(10)
        L.append(f"**{nombre_t}** ({int((~dfn['concordante']).sum())} discordant branches of {len(dfn)})\n")
        L.append("| ID | Leaves | Consensus taxon | Jaccard | Diff. | Match | Extra / Missing |")
        L.append("|---|---|---|---|---|---|---|")
        for _, r in d_.iterrows():
            sf = (r["sobran_vs_pareja"][:60] + " / " + r["faltan_vs_pareja"][:60]).replace("|", "/")
            L.append(f"| {r['id']} | {r['n_hojas']} | {r['taxon_consenso']} | {r['jaccard_mejor_pareja']} | {r['hojas_de_diferencia']} | {r['pareja_en_otro_arbol']} | {sf} |")
        L.append("")
    L.append("## 4. Largest identical blocks (best resolved)\n")
    L.append("Branches whose whole subtree has the same topology in both trees (`compatibility = 1`), outermost only.\n")
    ident = _bloques_identicos(datos, an)
    L.append("| ID (VAMPhyRE) | Leaves | Consensus taxon (purity) |")
    L.append("|---|---|---|")
    for nn, nid, lab, fr in ident[:15]:
        L.append(f"| {nid} | {nn} | {lab} ({fr:.0%}) |")
    L.append(f"\nIdentical blocks with ≥3 leaves: {len(ident)}, covering {sum(x[0] for x in ident)} of {N} genomes.\n")
    L.append("## 5. Taxonomic groups\n")
    if dfg.empty:
        L.append("(no taxonomy table)\n")
    for rl in rangos:
        nr = LETRA_A_NOMBRE[rl]
        sub_ = dfg[dfg["rango"] == nr].dropna(subset=["concordancia_interna"]) if not dfg.empty else dfg
        if sub_.empty:
            continue
        L.append(f"### {nr.capitalize()} (groups with ≥{min_grupo} genomes: {len(sub_)})\n")
        ambos = int((sub_["monofiletico_V"] & sub_["monofiletico_G"]).sum())
        L.append(f"Monophyletic in both trees: **{ambos}/{len(sub_)}**. "
                 f"VAMPhyRE only: {int((sub_['monofiletico_V'] & ~sub_['monofiletico_G']).sum())}; "
                 f"GTDB only: {int((~sub_['monofiletico_V'] & sub_['monofiletico_G']).sum())}; "
                 f"neither: {int((~sub_['monofiletico_V'] & ~sub_['monofiletico_G']).sum())}.\n")
        mejor = sub_.sort_values(["concordancia_interna", "n"], ascending=[False, False]).head(8)
        peor = sub_.sort_values(["concordancia_interna", "n"], ascending=[True, False]).head(8)
        for titulo, t_ in (("Best resolved", mejor), ("Most conflicting", peor)):
            L.append(f"**{titulo}**\n")
            L.append("| Group | n | Mono V | Mono G | Fragments V/G | Internal concordance | Rogues |")
            L.append("|---|---|---|---|---|---|---|")
            for _, r in t_.iterrows():
                L.append(f"| {r['grupo']} | {r['n']} | {'yes' if r['monofiletico_V'] else 'no'} | {'yes' if r['monofiletico_G'] else 'no'} | "
                         f"{r['fragmentos_V']}/{r['fragmentos_G']} | {r['concordancia_interna']:.0%} | {r['rebeldes_en_grupo']} |")
            L.append("")
    L.append("## 6. Output files\n")
    for suf, desc in (("_hojas_rebeldes.csv", "per-genome ranking"), ("_nodos_vamphyre.csv", "VAMPhyRE branches (IDs V#)"),
                      ("_nodos_gtdb.csv", "GTDB branches (IDs G#)"), ("_grupos.csv", "per taxonomic group metrics"),
                      ("_mapa_ids.tsv", "ID map table"), ("_rebeldes.png", "summary figure")):
        if os.path.exists(prefijo + suf):
            L.append(f"- `{os.path.basename(prefijo + suf)}`: {desc}")
    L.append("\nTo collapse/expand specific branches, use the CSV IDs with `--colapsar-nodos V12,G40` / `--expandir-nodos`.")
    with open(ruta, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    log.info(f"[+] Report: {ruta}")


def escribir_tablas_y_reporte(datos, an, prefijo, max_rebeldes, rango_grupos, min_grupo,
                              colV, colG, modo_colapso):
    """Write CSVs, summary figure and report. Returns the main DataFrames."""
    tax_ok = datos.esquema.hay_tax
    tabla_hojas(datos, an, max_rebeldes).to_csv(prefijo + "_hojas_rebeldes.csv", index=False)
    dfnV = tabla_nodos(datos, an, an.rV, an.infoV, an.idsG)
    dfnG = tabla_nodos(datos, an, an.rG, an.infoG, an.idsV)
    dfnV.to_csv(prefijo + "_nodos_vamphyre.csv", index=False)
    dfnG.to_csv(prefijo + "_nodos_gtdb.csv", index=False)
    rangos = [rango_letra(r) for r in rango_grupos.split(",") if r.strip()] if tax_ok else []
    dfg = tabla_grupos(datos, an, rangos, min_grupo)
    if not dfg.empty:
        dfg.to_csv(prefijo + "_grupos.csv", index=False)
    if an.pasos or not dfg.empty:
        figura_resumen(datos, an, dfg, rangos, prefijo + "_rebeldes.png")
    escribir_reporte(prefijo + "_reporte.md", prefijo, datos, an, dfnV, dfnG, dfg, rangos,
                     colV, colG, modo_colapso, min_grupo)
    return dfnV, dfnG, dfg
