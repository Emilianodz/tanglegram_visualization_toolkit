"""Normalisation of genome/leaf identifiers."""
import re


def limpiar_etiqueta(nombre):
    """Standardise labels (same rules as formate_tree_generation.py)."""
    nombre = str(nombre).strip()
    nombre = re.sub(r"^\[\s*\d+\s*\]\s*", "", nombre)
    nombre = nombre.lstrip(">").replace("#", "").replace("=", "_").replace(" ", "_")
    return re.sub(r"_{2,}", "_", nombre).strip("_")


def acc_gcf(s):
    m = re.search(r"(GC[AF]_\d+\.\d+)", str(s))
    return m.group(1) if m else None
