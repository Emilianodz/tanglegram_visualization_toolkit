"""Identifier harmonisation: trees and taxonomy table in a single namespace
(the leaf names of the VAMPhyRE tree).

Mapping strategies, from most to least direct: exact match, alias TSV table
(--mapa-ids), FASTA headers (--fasta-dir) and topology isomorphism with the
original GTDB tree (--tree-gtdb-original)."""
import collections
import logging
import os

import pandas as pd

from .errors import ErrorDatos
from .names import acc_gcf, limpiar_etiqueta
from .taxonomy import genero_base, parsear_tax
from .tree import cargar_arbol, hojas_de, podar, postorden

log = logging.getLogger(__name__)


def detectar_separador(ruta):
    """Column separator of a text table, decided from its first (header) line.

    A tab always wins: GTDB-Tk summaries are tab-separated but some headers contain commas
    (e.g. `other_related_references(genome_id,species_name,...)`), which fools pandas'
    automatic sniffing into reading the file as a single CSV column. Then comma, then
    semicolon; None (let pandas guess) only if the header has none of them."""
    with open(ruta, encoding="utf-8-sig", errors="replace") as f:
        cabecera = f.readline()
    for sep in ("\t", ",", ";"):
        if sep in cabecera:
            return sep
    return None


def leer_tabla(ruta, nrows=None, **kw):
    """Read an xlsx/xls/csv/tsv table with a reliable separator (see `detectar_separador`)."""
    if str(ruta).lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(ruta, nrows=nrows, **kw)
    sep = detectar_separador(ruta)
    if sep is None:
        return pd.read_csv(ruta, sep=None, engine="python", nrows=nrows, **kw)
    return pd.read_csv(ruta, sep=sep, nrows=nrows, **kw)


def alias_desde_mapa(ruta, tips):
    df = leer_tabla(ruta, header=None, dtype=str)
    a = {}
    for x, y in zip(df.iloc[:, 0], df.iloc[:, 1]):
        x, y = limpiar_etiqueta(x), limpiar_etiqueta(y)
        if y in tips:
            a[x] = y
        elif x in tips:
            a[y] = x
    return a


def alias_desde_fasta(carpeta, tips):
    a = {}
    for archivo in os.listdir(carpeta):
        if archivo.endswith((".fna", ".fasta", ".fa")):
            with open(os.path.join(carpeta, archivo)) as f:
                linea = f.readline().strip()
            if linea.startswith(">"):
                enc = limpiar_etiqueta(linea)
                if enc in tips:
                    a[limpiar_etiqueta(os.path.splitext(archivo)[0])] = enc
    return a


def _firmas(raiz, dec=5):
    est, lon = {}, {}
    for n in postorden(raiz):
        if n.es_hoja:
            est[n] = "L"
            lon[n] = f"L{round(n.largo, dec):.{dec}f}"
        else:
            est[n] = "(" + ",".join(sorted(est[h] for h in n.hijos)) + ")"
            lon[n] = "(" + ",".join(sorted(lon[h] for h in n.hijos)) + ")"
    return est, lon


def alias_por_topologia(r_orig, r_canon):
    """Match the leaves of two trees with the same topology (e.g. the original GTDB
    tree with GCF_* names and its relabelled NZ_* version). Returns (alias, ambiguous)."""
    eA, lA = _firmas(r_orig)
    eB, lB = _firmas(r_canon)
    pares = []

    def emparejar(a, b, amb):
        if a.es_hoja:
            pares.append((a.nombre, b.nombre, amb))
            return
        ga, gb = collections.defaultdict(list), collections.defaultdict(list)
        for h in a.hijos:
            ga[eA[h]].append(h)
        for h in b.hijos:
            gb[eB[h]].append(h)
        if {k: len(v) for k, v in ga.items()} != {k: len(v) for k, v in gb.items()}:
            raise ValueError("different topologies")
        for k in ga:
            xa, xb = ga[k], gb[k]
            if len(xa) == 1:
                emparejar(xa[0], xb[0], amb)
                continue
            sa, sb = collections.defaultdict(list), collections.defaultdict(list)
            for h in xa:
                sa[lA[h]].append(h)
            for h in xb:
                sb[lB[h]].append(h)
            resto_a, resto_b = [], []
            for s in sa:
                u, v = sa[s], sb.get(s, [])
                m = min(len(u), len(v))
                for i in range(m):
                    emparejar(u[i], v[i], amb or len(u) > 1)
                resto_a += u[m:]
                resto_b += v[m:]
            for s in sb:
                if s not in sa:
                    resto_b += sb[s]
            for u, v in zip(resto_a, resto_b):
                emparejar(u, v, True)

    emparejar(r_orig, r_canon, False)
    alias = {x: y for x, y, _ in pares}
    ambiguos = {x for x, _, am in pares if am}
    return alias, ambiguos


def ligar_tabla(df, comunes, alias, ambiguos, prefijo):
    """Link the table rows [id, taxonomy] to the tree leaves.

    `comunes`: leaf names; `alias`: {alternative_name: leaf}; `ambiguos`: names with an
    uncertain match. Writes `<prefijo>_mapa_ids.tsv` and returns {leaf: taxonomy}.
    Shared by the two-tree comparison and the single-tree pipeline."""
    clasif, filas = {}, []
    por_acc = collections.defaultdict(set)
    for k, v in alias.items():
        a = acc_gcf(k)
        if a:
            por_acc[a].add(v)
    n_ok = 0
    for gid, cl in zip(df.iloc[:, 0], df.iloc[:, 1]):
        gl = limpiar_etiqueta(gid)
        if gl in comunes:
            canon_n = gl
        elif gl in alias and alias[gl] in comunes:
            canon_n = alias[gl]
        elif acc_gcf(gl) and len(por_acc.get(acc_gcf(gl), ())) == 1:
            canon_n = next(iter(por_acc[acc_gcf(gl)]))
        else:
            continue
        if canon_n in comunes:
            clasif[canon_n] = cl
            filas.append((gid, canon_n, "ambiguous" if gl in ambiguos else "ok"))
            n_ok += 1
    log.info(f"    [+] Taxonomy table: {n_ok} of {len(df)} rows linked to leaves ({len(comunes)} leaves; "
             f"{len(comunes) - len(clasif)} without taxonomy)")
    if n_ok < 0.5 * len(comunes):
        log.warning(f"    [!] Only {len(clasif)} of {len(comunes)} genomes have taxonomy: the table IDs do not "
                    "match the leaf names, so there will be no colours or group names. "
                    "Provide the original GTDB tree (--tree-gtdb-original), an ID map (--mapa-ids) "
                    "or the FASTA folder (--fasta-dir).")
    # sanity check: the genus from the table should appear in the leaf name
    ok = tot = 0
    for nm, cl in clasif.items():
        g = parsear_tax(cl).get("g")
        if g:
            tot += 1
            ok += genero_base(g).lower() in nm.lower()
    if tot:
        pct = 100.0 * ok / tot
        log.info(f"    [i] Mapping check: the table genus appears in the leaf name in {pct:.1f}% of cases")
        if pct < 80:
            log.warning("    [!] Low percentage: the ID mapping may be wrong. Check --mapa-ids / --fasta-dir.")
    n_amb = sum(1 for f in filas if f[2] == "ambiguous")
    if n_amb:
        log.info(f"    [i] {n_amb} leaves with an ambiguous match (identical sisters); RF is not affected.")
    pd.DataFrame(filas, columns=["id_tabla", "hoja_en_arboles", "emparejamiento"]).to_csv(
        prefijo + "_mapa_ids.tsv", sep="\t", index=False)
    log.info(f"    [+] ID map table: {prefijo}_mapa_ids.tsv")
    return clasif


def homogeneizar(rV, rG, df, prefijo, mapa_ids=None, fasta_dir=None, tree_gtdb_original=None):
    """Bring both trees and the table into the same namespace (the VAMPhyRE one)."""
    log.info("[-] Harmonising identifiers...")
    tipsV = {n.nombre for n in hojas_de(rV)}
    tipsG = {n.nombre for n in hojas_de(rG)}
    alias, ambiguos, fuentes = {}, set(), []

    def sumar(d, etiqueta):
        nuevos = {k: v for k, v in d.items() if k not in alias}
        alias.update(nuevos)
        fuentes.append(f"{etiqueta}: {len(nuevos)} alias")

    if mapa_ids:
        sumar(alias_desde_mapa(mapa_ids, tipsV), "--mapa-ids")
    if fasta_dir:
        sumar(alias_desde_fasta(fasta_dir, tipsV), "--fasta-dir")
    if tree_gtdb_original:
        r_orig = cargar_arbol(tree_gtdb_original)
        if {n.nombre for n in hojas_de(r_orig)} != tipsG and tipsG <= tipsV:
            try:
                al, amb = alias_por_topologia(r_orig, rG)
                sumar({k: v for k, v in al.items() if v in tipsV}, "topology (--tree-gtdb-original)")
                ambiguos |= amb
            except ValueError:
                log.warning("    [!] The original and relabelled GTDB trees do not share the same topology/root; "
                            "the map could not be inferred. Use --fasta-dir or --mapa-ids.")
    for f in fuentes:
        log.info("    [i] " + f)

    # GTDB tree -> canonical names
    sin_resolver = []
    for h in hojas_de(rG):
        if h.nombre in tipsV:
            continue
        nuevo = alias.get(h.nombre)
        if nuevo is None:
            sin_resolver.append(h.nombre)
        else:
            h.nombre = nuevo
    tipsG = {n.nombre for n in hojas_de(rG) if n.nombre not in sin_resolver}
    comunes = tipsV & tipsG
    log.info(f"    [+] Leaves VAMPhyRE: {len(tipsV)} | GTDB: {len(tipsG) + len(sin_resolver)} | shared: {len(comunes)}")
    if len(comunes) < 4:
        raise ErrorDatos("Fewer than 4 leaves shared between the trees: check the names (use --mapa-ids or --fasta-dir).")
    if len(comunes) < len(tipsV) or sin_resolver:
        log.warning(f"    [!] Pruning both trees to the {len(comunes)} shared leaves "
                    f"(unmatched VAMPhyRE: {len(tipsV - comunes)}, unmatched GTDB: {len(sin_resolver) + len(tipsG - comunes)}).")
    rV = podar(rV, comunes)
    rG = podar(rG, comunes)

    # Table -> canonical names
    clasif = ligar_tabla(df, comunes, alias, ambiguos, prefijo) if df is not None else {}
    return rV, rG, clasif, len(comunes)


def homogeneizar_arbol(raiz, df, prefijo, mapa_ids=None, fasta_dir=None, arbol_referencia=None):
    """Single tree: link the table [id, taxonomy] to the leaves of `raiz`. Returns {leaf: taxonomy}.

    Table names are resolved as with two trees: exact match, GCF/GCA accession,
    ID map, FASTA, or a reference tree with the same topology."""
    log.info("[-] Linking the taxonomy table to the leaves...")
    tips = {n.nombre for n in hojas_de(raiz)}
    alias, ambiguos, fuentes = {}, set(), []

    def sumar(d, etiqueta):
        nuevos = {k: v for k, v in d.items() if k not in alias}
        alias.update(nuevos)
        fuentes.append(f"{etiqueta}: {len(nuevos)} alias")

    if mapa_ids:
        sumar(alias_desde_mapa(mapa_ids, tips), "--mapa-ids")
    if fasta_dir:
        sumar(alias_desde_fasta(fasta_dir, tips), "--fasta-dir")
    if arbol_referencia:
        r_ref = cargar_arbol(arbol_referencia)
        try:
            al, amb = alias_por_topologia(r_ref, raiz)
            sumar({k: v for k, v in al.items() if v in tips}, "topology (--arbol-referencia)")
            ambiguos |= amb
        except ValueError:
            log.warning("    [!] The reference tree does not share the tree's topology: the map could not be "
                        "inferred. Use --mapa-ids (e.g. the *_mapa_ids.tsv from a tanglegram run).")
    for f in fuentes:
        log.info("    [i] " + f)
    return ligar_tabla(df, tips, alias, ambiguos, prefijo)


def cargar_tabla(ruta, col_id="user_genome", col_tax="classification"):
    """Read the taxonomy table (xlsx/xls/csv/tsv) and return a DataFrame [id, taxonomy]."""
    df = leer_tabla(ruta)
    df.columns = [str(c).strip() for c in df.columns]
    faltan = [(c, flag) for c, flag in ((col_id, "--col-id"), (col_tax, "--col-tax")) if c not in df.columns]
    if faltan:
        cols = list(df.columns)
        lineas = [f"The table has no column '{c}' (set it with {flag})." for c, flag in faltan]
        lineas.append(f"Available columns ({len(cols)}): {', '.join(map(str, cols[:15]))}"
                      + (", ..." if len(cols) > 15 else ""))
        if len(cols) == 1:
            lineas.append("Only one column was found: check that the file is really tab- or comma-separated.")
        sugerencias = []
        if any(f == "--col-id" for _, f in faltan) and cols:
            sugerencias.append(f"--col-id {cols[0]}")
        if any(f == "--col-tax" for _, f in faltan):
            con_tax = [c for c in cols if df[c].astype(str).str.contains("d__", regex=False).any()]
            if con_tax:
                sugerencias.append(f"--col-tax {con_tax[0]}")
        if sugerencias:
            lineas.append(f"Likely fix: {' '.join(sugerencias)} "
                          "(in the web app, choose them under \"Table columns\").")
        raise ErrorDatos(" ".join(lineas))
    return df[[col_id, col_tax]].dropna().astype(str)
