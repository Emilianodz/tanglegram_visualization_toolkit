"""Shared constants: taxonomic ranks, palettes and neutral colours."""

# GTDB ranks (prefix letter) and their names.
# Rank names are accepted in English or Spanish on input; they are displayed in English.
RANGOS = ["d", "p", "c", "o", "f", "g", "s"]
NOMBRE_RANGO = {"dominio": "d", "filo": "p", "clase": "c", "orden": "o",
                "familia": "f", "genero": "g", "género": "g", "especie": "s",
                "domain": "d", "phylum": "p", "class": "c", "order": "o",
                "family": "f", "genus": "g", "species": "s"}
LETRA_A_NOMBRE = {"d": "domain", "p": "phylum", "c": "class", "o": "order",
                  "f": "family", "g": "genus", "s": "species"}

# Label for genomes with no taxon at the chosen rank
SIN = "unclassified"

# Reference categorical palette (8 validated colours, fixed order, never cycled).
# The extra colours (9+) are NOT validated for colour-vision deficiency: use them with --resaltar.
PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
          "#6250d6", "#e34948"]
PALETA_EXTRA = ["#0e7490", "#a16207", "#be185d", "#4d7c0f", "#7c3aed", "#b91c1c",
                "#0369a1", "#9a3412", "#475569", "#15803d", "#c026d3", "#ca8a04"]

# Warning thresholds for --color-top:
#  > MAX_SEGURO_DALTONISMO : colours are no longer validated for colour-vision deficiency.
#  > MAX_DISTINGUIBLE      : beyond ~20 categories colours are hard to tell apart even with
#                            normal vision; prefer highlighting, collapsing or a higher rank.
MAX_SEGURO_DALTONISMO = len(PALETA)                       # 8
MAX_DISTINGUIBLE = len(PALETA) + len(PALETA_EXTRA)        # 20

# Neutral colours
GRIS_OTROS = "#9ca3af"
GRIS_SIN = "#d1d5db"
GRIS_ATENUADO = "#c9cdd3"
GRIS_RAMA = "#4b5563"
GRIS_MIXTO = "#6b7280"
TINTA = "#1f2937"
TINTA_SUAVE = "#52514e"
AZUL_SERIE = "#2a78d6"
NARANJA_SERIE = "#eb6834"
