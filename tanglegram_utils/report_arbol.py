"""Tables and report for a single tree: nodes (support, size, taxon) and
monophyly of each taxonomic group in that tree."""
import pandas as pd

from .constants import LETRA_A_NOMBRE
from .splits import fragmentos
from .taxonomy import etiqueta_consenso, rango_letra, sin_prefijo
from .tree import preorden


def tabla_nodos_arbol(raiz, tax, hay_tax):
    """One row per internal node: ID, leaf count, support (0-100), branch length and, with a table, consensus taxon."""
    filas = []
    for n in preorden(raiz):
        if n.es_hoja or n.padre is None:
            continue
        fila = {"id": n.id, "hojas": n.n, "soporte": None if n.soporte is None else round(n.soporte, 2),
                "largo_rama": round(n.largo, 6)}
        if hay_tax:
            tok, frac = etiqueta_consenso(n.hojas, tax)
            fila["taxon_consenso"], fila["fraccion_taxon"] = sin_prefijo(tok), round(frac, 3)
        filas.append(fila)
    return pd.DataFrame(filas)


def tabla_grupos_arbol(raiz, tax, conservar, rangos, min_grupo):
    """Monophyly of each taxonomic group in this tree (with its current root).

    fragmentos = minimum number of clades that exactly cover the group (1 = monophyletic);
    pureza = fraction of the containing clade (MRCA) that belongs to the group; intrusos = the rest."""
    filas = []
    for rl in rangos:
        miembros = {}
        for i in conservar:
            g = tax[i].get(rl)
            if g:
                miembros.setdefault(g, []).append(i)
        for g, lst in miembros.items():
            if len(lst) < min_grupo:
                continue
            S = 0
            for i in lst:
                S |= 1 << i
            mrca = min((n for n in preorden(raiz) if (n.mask & S) == S), key=lambda n: n.n)
            frag = fragmentos(raiz, S)
            filas.append({"rango": LETRA_A_NOMBRE[rl], "grupo": sin_prefijo(g), "n": len(lst),
                          "monofiletico": frag == 1, "fragmentos": frag,
                          "pureza": round(len(lst) / mrca.n, 3), "intrusos": mrca.n - len(lst)})
    df = pd.DataFrame(filas)
    if len(df):
        df = df.sort_values(["rango", "fragmentos", "n"], ascending=[True, False, False]).reset_index(drop=True)
    return df


def escribir_reporte_arbol(ruta, titulo, N, motivo_raiz, resumen, n_contraidos, umbral_contraer, dfg, hay_tax):
    L = [f"# Tree report: {titulo}", "", f"- Genomes: **{N}**", f"- Root: {motivo_raiz}"]
    if resumen:
        L += [f"- Node support (0-100 scale): minimum {resumen['minimo']:.1f}, median {resumen['mediana']:.1f}; "
              f"{resumen['pct_ge_70']:.0f}% of nodes with support ≥ 70 and {resumen['pct_ge_90']:.0f}% ≥ 90 "
              f"({resumen['con_soporte']} of {resumen['nodos_internos']} internal nodes carry support)"]
    else:
        L += ["- Node support: the tree has no support values (bootstrap or consensus support "
              "will be read from internal node labels when available)."]
    if umbral_contraer is not None:
        L += [f"- Nodes with support < {umbral_contraer:g} collapsed into polytomies: **{n_contraidos}**"]
    if hay_tax and dfg is not None and len(dfg):
        L += ["", "## Taxonomic groups"]
        for rango, sub in dfg.groupby("rango", sort=False):
            mono = int(sub["monofiletico"].sum())
            L += ["", f"**{rango}**: {mono} of {len(sub)} groups are monophyletic.", "",
                  "| group | n | fragments | purity | intruders |", "|---|---|---|---|---|"]
            for _, r in sub[~sub["monofiletico"]].head(15).iterrows():
                L.append(f"| {r['grupo']} | {r['n']} | {r['fragmentos']} | {r['pureza']} | {r['intrusos']} |")
            if (~sub["monofiletico"]).sum() > 15:
                L.append(f"| … | | | | ({(~sub['monofiletico']).sum() - 15} more groups in the CSV) |")
    with open(ruta, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")


def rangos_de(texto):
    """'familia,genero' -> ['f', 'g']."""
    return [rango_letra(r) for r in texto.split(",") if r.strip()]
