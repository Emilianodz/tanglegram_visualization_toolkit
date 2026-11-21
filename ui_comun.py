"""Interface components (Streamlit) shared by the app pages.

Interface only: the logic lives in `tanglegram_utils`."""
import io
import os
import tempfile

import streamlit as st
from PIL import Image

from tanglegram_utils.constants import RANGOS
from tanglegram_utils.servicio import (adivinar_columnas, columnas_tabla, diagnosticar_ids,
                                       inspeccionar_arbol, listar_taxa)

# Internal rank values (passed to the pipeline) and their English display names
RANGOS_UI = ["dominio", "filo", "clase", "orden", "familia", "genero", "especie"]
NOMBRE_UI = {"dominio": "domain", "filo": "phylum", "clase": "class", "orden": "order",
             "familia": "family", "genero": "genus", "especie": "species"}
# English labels for option values that the CLI keeps in Spanish
OPCIONES_UI = {"ninguno": "none", "numero": "number", "punto": "dot", "cladograma": "cladogram",
               "longitud": "branch lengths", "especie": "species", "nombre": "raw name",
               "ladderize": "ladderize", "original": "original"}
# Selectable ranks are phylum→genus: domain (everything) and species (1 genome each) are useless as filters
RANGOS_SELECCION = [r for r in RANGOS if r not in ("d", "s")]


def etiqueta_opcion(valor):
    """English label for an option value (values themselves stay as the CLI expects)."""
    return OPCIONES_UI.get(valor, valor)


def carpeta_trabajo():
    if "workdir" not in st.session_state:
        st.session_state["workdir"] = tempfile.mkdtemp(prefix="tanglegram_app_")
    return st.session_state["workdir"]


def guardar_subida(archivo):
    """Write the uploaded file into the temporary folder and return its path."""
    destino = os.path.join(carpeta_trabajo(), "entradas")
    os.makedirs(destino, exist_ok=True)
    ruta = os.path.join(destino, os.path.basename(archivo.name))
    if not os.path.exists(ruta) or os.path.getsize(ruta) != archivo.size:
        with open(ruta, "wb") as f:
            f.write(archivo.getbuffer())
    return ruta


@st.cache_data(show_spinner=False)
def columnas_cacheadas(ruta, tam):
    return columnas_tabla(ruta)


@st.cache_data(show_spinner=False)
def columnas_probables(ruta, tam):
    return adivinar_columnas(ruta)


@st.cache_data(show_spinner=False)
def cobertura_ids(ruta_arbol, tam_arbol, ruta_tabla, tam_tabla, col_id, col_tax):
    return diagnosticar_ids(ruta_arbol, ruta_tabla, col_id, col_tax)


@st.cache_data(show_spinner=False)
def taxa_de_tabla(ruta, tam, col_id, col_tax):
    return listar_taxa(ruta, col_id, col_tax)


@st.cache_data(show_spinner=False)
def info_arbol(ruta, tam):
    return inspeccionar_arbol(ruta)


@st.cache_data(show_spinner=False)
def vista_previa(ruta, mtime, ancho_max=2400):
    """Downscaled PNG for display in the browser (the original is downloaded in full)."""
    im = Image.open(ruta)
    im.thumbnail((ancho_max, im.height))
    buf = io.BytesIO()
    im.save(buf, "PNG")
    return buf.getvalue()


def lista(valores):
    """['a','b'] -> 'a,b'; empty -> None."""
    return ",".join(valores) if valores else None


def mostrar_imagen(ruta, clave):
    ajustar = st.toggle("Fit to page width", True, key=f"aj_{clave}")
    ancho = "stretch" if ajustar else st.slider("Width (px)", 600, 2400, 1400, 100, key=f"an_{clave}")
    st.image(vista_previa(ruta, os.path.getmtime(ruta)), width=ancho)


def crear_selector_taxa(taxa):
    """Return (selector_taxa, options, labels) for the taxa of the loaded table."""
    opciones = [t for r in RANGOS_SELECCION for t, _ in taxa.get(r, [])]
    etiquetas = {t: f"{t}  ·  {n}" for r in RANGOS_SELECCION for t, n in taxa.get(r, [])}

    def selector_taxa(etiqueta, clave, ayuda=None):
        if not opciones:
            st.caption(f"{etiqueta}: upload the table to enable this.")
            return None
        return lista(st.multiselect(etiqueta, opciones, key=clave, format_func=etiquetas.get, help=ayuda))

    return selector_taxa, opciones, etiquetas


def compartir(**archivos):
    """Make the files loaded on this page available to the other pages (only those that exist)."""
    st.session_state["compartido"] = {k: v for k, v in archivos.items() if v}


def compartidos():
    return st.session_state.get("compartido", {})
