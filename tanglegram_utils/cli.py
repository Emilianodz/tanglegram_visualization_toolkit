"""Command-line interface: translates arguments into pipeline options."""
import argparse
import logging
import shlex
import sys

from .collapse import OpcionesColapso
from .errors import ErrorDatos
from .logs import configurar_logging
from .models import FIGURAS, Entradas, OpcionesAnalisis, Salida
from .pipeline import ejecutar
from .plotting import EstiloFigura
from .taxonomy import rango_letra


def construir_parser():
    p = argparse.ArgumentParser(
        prog="generate_tanglegram.py",
        description="VAMPhyRE vs GTDB-Tk tanglegram + rogue-branch analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    g = p.add_argument_group("Inputs")
    g.add_argument("--tree-vamphyre", required=True, help="VAMPhyRE Newick tree")
    g.add_argument("--tree-gtdb", required=True, help="GTDB-Tk Newick tree (preferably the relabelled one)")
    g.add_argument("--out-plot", required=True, help="Output PNG image")
    g.add_argument("--figura", default="tanglegram", choices=FIGURAS,
                   help="tanglegram: both trees facing each other · vamphyre / gtdb: a single tree · "
                        "individuales: one PNG per tree (<out>_vamphyre.png and <out>_gtdb.png)")
    g.add_argument("--tabla", help="Taxonomy table (xlsx/tsv/csv; e.g. gtdbtk.ar53.summary)")
    g.add_argument("--col-id", default="user_genome", help="Table column with the genome ID")
    g.add_argument("--col-tax", default="classification", help="Table column with the d__;p__;... taxonomy")
    g.add_argument("--tree-gtdb-original", help="Original GTDB tree (GCF_* names) to infer the ID map from topology")
    g.add_argument("--fasta-dir", help="Folder of .fna files: maps file name -> header (as in formate_tree_generation.py)")
    g.add_argument("--mapa-ids", help="2-column TSV (table_id <tab> leaf_name)")
    g.add_argument("--prefijo-salida", help="Prefix for report/CSVs (default: --out-plot without extension)")

    g = p.add_argument_group("Colour and taxonomic selection")
    g.add_argument("--color-rango", default="familia",
                   help="Rank used for colour: domain|phylum|class|order|family|genus|species "
                        "(Spanish names dominio..especie and letters d,p,c,o,f,g,s also accepted)")
    g.add_argument("--color-top", type=int, default=8, help="Number of taxa with their own colour (the rest: grey 'other'). ≤8: colour-blind safe; 9-20: unvalidated colours; >20: generated colours, hard to tell apart")
    g.add_argument("--resaltar", help="Comma-separated taxa to highlight (e.g. g__Halorubrum,Methanosarcinales); the rest is dimmed")
    g.add_argument("--ocultar-lineas", help="Taxa whose connecting lines are NOT drawn")
    g.add_argument("--solo-taxa", help="Prune both trees to these taxa (and recompute everything on the subset)")

    g = p.add_argument_group("Clade collapsing")
    g.add_argument("--colapsar", default="ninguno", choices=["ninguno", "concordantes", "discordantes", "taxonomia"],
                   help="concordantes: clades with identical/compatible topology in both trees; discordantes: conflicting ones; "
                        "taxonomia: clades made of a single taxon (see --rango-colapso)")
    g.add_argument("--umbral", type=float, help="Clade compatibility (0-1). concordantes: >= threshold (default 1.0); discordantes: < threshold (default 0.5)")
    g.add_argument("--min-colapso", type=int, default=3, help="Minimum size (leaves) to collapse a node")
    g.add_argument("--max-colapso", type=int, help="Maximum size (leaves) to collapse a node")
    g.add_argument("--rango-colapso", help="Rank for --colapsar taxonomia (default: --color-rango)")
    g.add_argument("--colapsar-nodos", help="Node IDs to always collapse (e.g. V12,G40; see the node CSVs or --mostrar-ids)")
    g.add_argument("--expandir-nodos", help="Node IDs that must stay expanded")
    g.add_argument("--colapsar-taxa", help="Taxa (pure clades) to always collapse")
    g.add_argument("--expandir-taxa", help="Taxa that must stay expanded")

    g = p.add_argument_group("Appearance")
    g.add_argument("--forma", default="rectangular", choices=["rectangular", "circular"],
                   help="Single-tree layout (--figura vamphyre|gtdb|individuales); circular does not apply to the tanglegram")
    g.add_argument("--arco", type=float, default=360.0, help="Circular tree: degrees spanned (the gap is at the top)")
    g.add_argument("--soporte", default="ninguno", choices=["ninguno", "numero", "punto"],
                   help="Node support values (bootstrap/SH): number or dot (single trees only)")
    g.add_argument("--soporte-min", type=float, default=70.0,
                   help="Only support >= this value is shown (0-100 scale; 0-1 is rescaled automatically)")
    g.add_argument("--escala", default="cladograma", choices=["cladograma", "longitud"], help="Equal branches (cladograma) or real branch lengths (longitud)")
    g.add_argument("--etiquetas-hojas", action="store_true", help="Write the name of each leaf")
    g.add_argument("--texto-hoja", default="especie", choices=["especie", "nombre"], help="Leaf label text: species (especie) or raw name (nombre)")
    g.add_argument("--sin-etiquetas-grupos", action="store_true", help="Do not label pure clades/triangles")
    g.add_argument("--min-etiqueta", type=int, default=4, help="Minimum size of a pure clade to label it")
    g.add_argument("--mostrar-ids", action="store_true", help="Write the ID of large nodes")
    g.add_argument("--min-id", type=int, default=15, help="Minimum size to show a node ID")
    g.add_argument("--ancho-arbol", type=float, default=None,
                   help="Width of each tree (inches); default 4.5 in a tanglegram and 8 for a single tree")
    g.add_argument("--ancho-centro", type=float, default=5.0, help="Width of the connecting-line area (inches)")
    g.add_argument("--alto", type=float, help="Figure height (inches); default depends on the number of rows")
    g.add_argument("--dpi", type=int, default=150)
    g.add_argument("--grosor-linea", type=float, default=1.3, help="Base width of the connecting lines")
    g.add_argument("--grosor-rama", type=float, default=1.1)
    g.add_argument("--escala-triangulo", type=float, default=1.0)
    g.add_argument("--sin-desenredar", action="store_true", help="Do not reorder branches to reduce crossings")
    g.add_argument("--raiz-taxon", help="Root both trees (display only) on this taxon if it is monophyletic in both")
    g.add_argument("--formato-extra", help="Additional figure formats: svg,pdf")

    g = p.add_argument_group("Analysis")
    g.add_argument("--max-rebeldes", type=int, default=25, help="Steps of the greedy rogue-genome removal (0 = skip)")
    g.add_argument("--rango-grupos", default="familia,genero", help="Ranks for the group table (comma-separated)")
    g.add_argument("--min-grupo", type=int, default=3, help="Minimum group size to evaluate")
    g.add_argument("--silencioso", action="store_true", help="Show warnings only")
    return p


def opciones_desde_args(args):
    """Convert the argparse Namespace into the pipeline dataclasses."""
    rc = rango_letra(args.color_rango)
    entradas = Entradas(tree_vamphyre=args.tree_vamphyre, tree_gtdb=args.tree_gtdb, tabla=args.tabla,
                        col_id=args.col_id, col_tax=args.col_tax,
                        tree_gtdb_original=args.tree_gtdb_original, fasta_dir=args.fasta_dir,
                        mapa_ids=args.mapa_ids)
    op = OpcionesAnalisis(color_rango=args.color_rango, color_top=args.color_top, resaltar=args.resaltar,
                          ocultar_lineas=args.ocultar_lineas, solo_taxa=args.solo_taxa,
                          raiz_taxon=args.raiz_taxon, sin_desenredar=args.sin_desenredar,
                          max_rebeldes=args.max_rebeldes, rango_grupos=args.rango_grupos,
                          min_grupo=args.min_grupo)
    colapso = OpcionesColapso(modo=args.colapsar, umbral=args.umbral, min_colapso=args.min_colapso,
                              max_colapso=args.max_colapso,
                              rango=rango_letra(args.rango_colapso) if args.rango_colapso else rc,
                              colapsar_nodos=args.colapsar_nodos, expandir_nodos=args.expandir_nodos,
                              colapsar_taxa=args.colapsar_taxa, expandir_taxa=args.expandir_taxa)
    estilo = EstiloFigura(escala=args.escala, etiquetas_hojas=args.etiquetas_hojas, texto_hoja=args.texto_hoja,
                          etiquetas_grupos=not args.sin_etiquetas_grupos, min_etiqueta=args.min_etiqueta,
                          mostrar_ids=args.mostrar_ids, min_id=args.min_id, ancho_arbol=args.ancho_arbol,
                          ancho_centro=args.ancho_centro, alto=args.alto, dpi=args.dpi,
                          grosor_linea=args.grosor_linea, grosor_rama=args.grosor_rama,
                          escala_triangulo=args.escala_triangulo, forma=args.forma, arco=args.arco,
                          soporte=args.soporte, soporte_min=args.soporte_min)
    salida = Salida(out_plot=args.out_plot, prefijo=args.prefijo_salida, formato_extra=args.formato_extra,
                    figura=args.figura)
    return entradas, op, colapso, estilo, salida


def args_desde_opciones(entradas, op, colapso, estilo, salida):
    """Inverse of `opciones_desde_args`: returns {argparse_dest: value}."""
    return {
        "tree_vamphyre": entradas.tree_vamphyre, "tree_gtdb": entradas.tree_gtdb,
        "out_plot": salida.out_plot, "tabla": entradas.tabla, "col_id": entradas.col_id,
        "col_tax": entradas.col_tax, "tree_gtdb_original": entradas.tree_gtdb_original,
        "fasta_dir": entradas.fasta_dir, "mapa_ids": entradas.mapa_ids,
        "prefijo_salida": salida.prefijo, "color_rango": op.color_rango, "color_top": op.color_top,
        "resaltar": op.resaltar, "ocultar_lineas": op.ocultar_lineas, "solo_taxa": op.solo_taxa,
        "colapsar": colapso.modo, "umbral": colapso.umbral, "min_colapso": colapso.min_colapso,
        "max_colapso": colapso.max_colapso,
        "rango_colapso": colapso.rango if colapso.rango != rango_letra(op.color_rango) else None,
        "colapsar_nodos": colapso.colapsar_nodos, "expandir_nodos": colapso.expandir_nodos,
        "colapsar_taxa": colapso.colapsar_taxa, "expandir_taxa": colapso.expandir_taxa,
        "escala": estilo.escala, "etiquetas_hojas": estilo.etiquetas_hojas,
        "texto_hoja": estilo.texto_hoja, "sin_etiquetas_grupos": not estilo.etiquetas_grupos,
        "min_etiqueta": estilo.min_etiqueta, "mostrar_ids": estilo.mostrar_ids, "min_id": estilo.min_id,
        "ancho_arbol": estilo.ancho_arbol, "ancho_centro": estilo.ancho_centro, "alto": estilo.alto,
        "dpi": estilo.dpi, "grosor_linea": estilo.grosor_linea, "grosor_rama": estilo.grosor_rama,
        "escala_triangulo": estilo.escala_triangulo, "sin_desenredar": op.sin_desenredar,
        "forma": estilo.forma, "arco": estilo.arco, "soporte": estilo.soporte, "soporte_min": estilo.soporte_min,
        "raiz_taxon": op.raiz_taxon, "formato_extra": salida.formato_extra, "figura": salida.figura,
        "max_rebeldes": op.max_rebeldes, "rango_grupos": op.rango_grupos, "min_grupo": op.min_grupo,
    }


def formatear_comando(valores, defaults, obligatorios, programa):
    """Build the command line with flags whose value differs from the default (required ones always)."""
    partes = []
    for dest, val in valores.items():
        flag = "--" + dest.replace("_", "-")
        if dest in obligatorios:
            partes.append(f"{flag} {shlex.quote(str(val))}")
        elif val is None or val == defaults[dest]:
            continue
        elif val is True:
            partes.append(flag)
        else:
            partes.append(f"{flag} {shlex.quote(str(val))}")
    return programa + " \\\n  " + " \\\n  ".join(partes)


def comando_equivalente(entradas, op, colapso, estilo, salida, programa="python generate_tanglegram.py"):
    """Command line that reproduces these options (only flags that differ from the
    defaults), recording how a figure was generated."""
    obligatorios = ("tree_vamphyre", "tree_gtdb", "out_plot")
    base = vars(construir_parser().parse_args(["--tree-vamphyre", "x", "--tree-gtdb", "x", "--out-plot", "x"]))
    return formatear_comando(args_desde_opciones(entradas, op, colapso, estilo, salida), base, obligatorios, programa)


def main(argv=None):
    """Entry point. Returns 0 on success, 1 on a data error."""
    args = construir_parser().parse_args(argv)
    configurar_logging(silencioso=args.silencioso)
    try:
        ejecutar(*opciones_desde_args(args))
    except ErrorDatos as e:
        logging.getLogger("tanglegram_utils").error(f"[!] {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
