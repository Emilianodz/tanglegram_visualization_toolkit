"""Command-line interface for the single tree (generate_tree.py)."""
import argparse
import logging
import sys

from .cli import formatear_comando
from .collapse import OpcionesColapso
from .errors import ErrorDatos
from .logs import configurar_logging
from .models import EntradaArbol, OpcionesArbol, Salida
from .pipeline_arbol import ejecutar_arbol
from .plotting import EstiloFigura
from .taxonomy import rango_letra


def construir_parser():
    p = argparse.ArgumentParser(
        prog="generate_tree.py",
        description="Single tree (rectangular or circular) with taxonomy colours, clade collapsing, support values and iTOL export",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    g = p.add_argument_group("Inputs")
    g.add_argument("--arbol", required=True, help="Newick tree (VAMPhyRE, GTDB-Tk, consensus...)")
    g.add_argument("--out-plot", required=True, help="Output image (png; also svg/pdf with --formato-extra)")
    g.add_argument("--tabla", help="Taxonomy table (xlsx/tsv/csv; e.g. gtdbtk.ar53.summary) for colours and groups")
    g.add_argument("--col-id", default="user_genome", help="Table column with the genome ID")
    g.add_argument("--col-tax", default="classification", help="Table column with the d__;p__;... taxonomy")
    g.add_argument("--mapa-ids", help="TSV table_id <tab> leaf_name (the *_mapa_ids.tsv written by the tanglegram works)")
    g.add_argument("--fasta-dir", help="Folder of .fna files: maps file name -> header")
    g.add_argument("--arbol-referencia", help="Tree with the table names and the SAME topology (e.g. the original GTDB tree)")
    g.add_argument("--titulo", help="Figure title (default: the file name)")
    g.add_argument("--prefijo-salida", help="Prefix for tables and report (default: --out-plot without extension)")

    g = p.add_argument_group("Tree")
    g.add_argument("--raiz", default="original", choices=["original", "punto-medio", "taxon"],
                   help="original: as in the file · punto-medio: midpoint of the longest path · taxon: outgroup")
    g.add_argument("--raiz-taxon", help="Outgroup when --raiz taxon (must form a clade)")
    g.add_argument("--contraer-soporte", type=float,
                   help="Nodes with support < this value (0-100 scale) become polytomies, as in a consensus tree")
    g.add_argument("--orden", default="ladderize", choices=["ladderize", "original"],
                   help="ladderize: branches sorted by size · original: as in the file")

    g = p.add_argument_group("Colour and taxonomic selection")
    g.add_argument("--color-rango", default="familia", help="Rank used for colour: domain|phylum|class|order|family|genus|species (Spanish names and letters d,p,c,o,f,g,s also accepted)")
    g.add_argument("--color-top", type=int, default=8, help="Number of taxa with their own colour (the rest: grey 'other'). ≤8: colour-blind safe; 9-20: unvalidated colours; >20: generated colours, hard to tell apart")
    g.add_argument("--resaltar", help="Comma-separated taxa to highlight; the rest is dimmed")
    g.add_argument("--solo-taxa", help="Draw only these taxa (prunes the tree)")

    g = p.add_argument_group("Clade collapsing")
    g.add_argument("--colapsar", default="ninguno", choices=["ninguno", "taxonomia"],
                   help="taxonomia: clades made of a single taxon (see --rango-colapso)")
    g.add_argument("--min-colapso", type=int, default=3, help="Minimum size (leaves) to collapse a node")
    g.add_argument("--max-colapso", type=int, help="Maximum size (leaves) to collapse a node")
    g.add_argument("--rango-colapso", help="Rank for --colapsar taxonomia (default: --color-rango)")
    g.add_argument("--colapsar-nodos", help="Node IDs to always collapse (e.g. N12,N40; see the node CSV)")
    g.add_argument("--expandir-nodos", help="Node IDs that must stay expanded")
    g.add_argument("--colapsar-taxa", help="Taxa (pure clades) to always collapse")
    g.add_argument("--expandir-taxa", help="Taxa that must stay expanded")

    g = p.add_argument_group("Appearance")
    g.add_argument("--forma", default="circular", choices=["rectangular", "circular"])
    g.add_argument("--arco", type=float, default=360.0, help="Circular: degrees spanned by the tree (the gap is at the top)")
    g.add_argument("--soporte", default="ninguno", choices=["ninguno", "numero", "punto"],
                   help="Node support values (bootstrap/SH): number (numero) or dot (punto)")
    g.add_argument("--soporte-min", type=float, default=70.0,
                   help="Only support >= this value is shown (0-100 scale; 0-1 is rescaled automatically)")
    g.add_argument("--escala", default="cladograma", choices=["cladograma", "longitud"], help="Equal branches (cladograma) or real branch lengths (longitud)")
    g.add_argument("--etiquetas-hojas", action="store_true", help="Write the name of each leaf")
    g.add_argument("--texto-hoja", default="especie", choices=["especie", "nombre"], help="Leaf label text: species (especie) or raw name (nombre)")
    g.add_argument("--sin-etiquetas-grupos", action="store_true", help="Do not label pure clades/triangles")
    g.add_argument("--min-etiqueta", type=int, default=4, help="Minimum size of a pure clade to label it")
    g.add_argument("--mostrar-ids", action="store_true", help="Write the ID of large nodes")
    g.add_argument("--min-id", type=int, default=15, help="Minimum size to show a node ID")
    g.add_argument("--ancho-arbol", type=float, help="Radius (circular) or width (rectangular) in inches; default: automatic")
    g.add_argument("--alto", type=float, help="Figure height (inches); default: automatic")
    g.add_argument("--dpi", type=int, default=150)
    g.add_argument("--grosor-rama", type=float, default=1.1)
    g.add_argument("--escala-triangulo", type=float, default=1.0)
    g.add_argument("--formato-extra", help="Additional figure formats: svg,pdf")

    g = p.add_argument_group("Tables")
    g.add_argument("--rango-grupos", default="familia,genero", help="Ranks for the monophyly table (comma-separated)")
    g.add_argument("--min-grupo", type=int, default=3, help="Minimum group size to evaluate")
    g.add_argument("--sin-itol", action="store_true",
                   help="Do not export iTOL files (_itol.tree, .txt annotations and _itol.zip)")
    g.add_argument("--silencioso", action="store_true", help="Show warnings only")
    return p


def opciones_desde_args(args):
    """argparse Namespace -> (EntradaArbol, OpcionesArbol, OpcionesColapso, EstiloFigura, Salida)."""
    entrada = EntradaArbol(arbol=args.arbol, tabla=args.tabla, col_id=args.col_id, col_tax=args.col_tax,
                           mapa_ids=args.mapa_ids, fasta_dir=args.fasta_dir,
                           arbol_referencia=args.arbol_referencia, titulo=args.titulo)
    op = OpcionesArbol(color_rango=args.color_rango, color_top=args.color_top, resaltar=args.resaltar,
                       solo_taxa=args.solo_taxa, raiz=args.raiz, raiz_taxon=args.raiz_taxon,
                       contraer_soporte=args.contraer_soporte, orden=args.orden,
                       rango_grupos=args.rango_grupos, min_grupo=args.min_grupo, itol=not args.sin_itol)
    colapso = OpcionesColapso(modo=args.colapsar, min_colapso=args.min_colapso, max_colapso=args.max_colapso,
                              rango=rango_letra(args.rango_colapso or args.color_rango),
                              colapsar_nodos=args.colapsar_nodos, expandir_nodos=args.expandir_nodos,
                              colapsar_taxa=args.colapsar_taxa, expandir_taxa=args.expandir_taxa)
    estilo = EstiloFigura(escala=args.escala, etiquetas_hojas=args.etiquetas_hojas, texto_hoja=args.texto_hoja,
                          etiquetas_grupos=not args.sin_etiquetas_grupos, min_etiqueta=args.min_etiqueta,
                          mostrar_ids=args.mostrar_ids, min_id=args.min_id, ancho_arbol=args.ancho_arbol,
                          alto=args.alto, dpi=args.dpi, grosor_rama=args.grosor_rama,
                          escala_triangulo=args.escala_triangulo, forma=args.forma, arco=args.arco,
                          soporte=args.soporte, soporte_min=args.soporte_min)
    salida = Salida(out_plot=args.out_plot, prefijo=args.prefijo_salida, formato_extra=args.formato_extra)
    return entrada, op, colapso, estilo, salida


def args_desde_opciones(entrada, op, colapso, estilo, salida):
    """Inverse of `opciones_desde_args`: {argparse_dest: value}."""
    return {
        "arbol": entrada.arbol, "out_plot": salida.out_plot, "tabla": entrada.tabla, "col_id": entrada.col_id,
        "col_tax": entrada.col_tax, "mapa_ids": entrada.mapa_ids, "fasta_dir": entrada.fasta_dir,
        "arbol_referencia": entrada.arbol_referencia, "titulo": entrada.titulo,
        "prefijo_salida": salida.prefijo, "raiz": op.raiz, "raiz_taxon": op.raiz_taxon,
        "contraer_soporte": op.contraer_soporte, "orden": op.orden,
        "color_rango": op.color_rango, "color_top": op.color_top, "resaltar": op.resaltar, "solo_taxa": op.solo_taxa,
        "colapsar": colapso.modo, "min_colapso": colapso.min_colapso, "max_colapso": colapso.max_colapso,
        "rango_colapso": colapso.rango if colapso.rango != rango_letra(op.color_rango) else None,
        "colapsar_nodos": colapso.colapsar_nodos, "expandir_nodos": colapso.expandir_nodos,
        "colapsar_taxa": colapso.colapsar_taxa, "expandir_taxa": colapso.expandir_taxa,
        "forma": estilo.forma, "arco": estilo.arco, "soporte": estilo.soporte, "soporte_min": estilo.soporte_min,
        "escala": estilo.escala, "etiquetas_hojas": estilo.etiquetas_hojas, "texto_hoja": estilo.texto_hoja,
        "sin_etiquetas_grupos": not estilo.etiquetas_grupos, "min_etiqueta": estilo.min_etiqueta,
        "mostrar_ids": estilo.mostrar_ids, "min_id": estilo.min_id, "ancho_arbol": estilo.ancho_arbol,
        "alto": estilo.alto, "dpi": estilo.dpi, "grosor_rama": estilo.grosor_rama,
        "escala_triangulo": estilo.escala_triangulo, "formato_extra": salida.formato_extra,
        "rango_grupos": op.rango_grupos, "min_grupo": op.min_grupo, "sin_itol": not op.itol,
    }


def comando_equivalente(entrada, op, colapso, estilo, salida, programa="python generate_tree.py"):
    """Command line that reproduces these options (only flags that differ from the defaults)."""
    base = vars(construir_parser().parse_args(["--arbol", "x", "--out-plot", "x"]))
    return formatear_comando(args_desde_opciones(entrada, op, colapso, estilo, salida), base,
                             ("arbol", "out_plot"), programa)


def main(argv=None):
    """Entry point. Returns 0 on success, 1 on a data error."""
    args = construir_parser().parse_args(argv)
    configurar_logging(silencioso=args.silencioso)
    try:
        ejecutar_arbol(*opciones_desde_args(args))
    except ErrorDatos as e:
        logging.getLogger("tanglegram_utils").error(f"[!] {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
