"""Web app (Streamlit): VAMPhyRE vs GTDB-Tk tanglegram and single trees.

    streamlit run app.py

This file only defines navigation; each page lives in `paginas/` and all logic in
`tanglegram_utils`. Files are processed in a local temporary folder; nothing leaves your computer."""
import matplotlib

matplotlib.use("Agg")  # before pyplot is imported (by the package): no GUI window

import streamlit as st  # noqa: E402

st.set_page_config(page_title="VAMPhyRE / GTDB-Tk trees", page_icon="🌿", layout="wide")

pagina = st.navigation([
    st.Page("paginas/tanglegram.py", title="Tanglegram (two trees)", icon="🔀", default=True),
    st.Page("paginas/arbol_individual.py", title="Single tree", icon="🌳"),
])
pagina.run()
