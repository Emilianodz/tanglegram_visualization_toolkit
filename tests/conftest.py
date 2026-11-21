import sys
from pathlib import Path

import pytest

# allows running `pytest` from the project root without installing the package
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tanglegram_utils.tree import desde_cadena, hojas_de, indexar  # noqa: E402

# 12 leaves, trifurcating root (like the real trees)
NWK_V = ("(((a1:1,a2:1):1,(a3:1,a4:1):1):1,((b1:1,b2:1):1,(b3:1,b4:1):1):1,"
         "((c1:1,c2:1):1,(c3:1,c4:1):1):1);")
# same as V but with a3 <-> b1 swapped
NWK_G = ("(((a1:1,a2:1):1,(b1:1,a4:1):1):1,((a3:1,b2:1):1,(b3:1,b4:1):1):1,"
         "((c1:1,c2:1):1,(c3:1,c4:1):1):1);")


def preparar(nwk_a, nwk_b):
    """Return (rA, rB, idx, todo) with both trees indexed."""
    rA, rB = desde_cadena(nwk_a), desde_cadena(nwk_b)
    nombres = sorted(h.nombre for h in hojas_de(rA))
    idx = {n: i for i, n in enumerate(nombres)}
    indexar(rA, idx)
    indexar(rB, idx)
    return rA, rB, idx, (1 << len(nombres)) - 1


@pytest.fixture
def arboles():
    return preparar(NWK_V, NWK_G)


@pytest.fixture
def datos_sinteticos(tmp_path):
    """Minimal input files for the full pipeline."""
    (tmp_path / "v.newick").write_text(NWK_V + "\n")
    (tmp_path / "g.newick").write_text(NWK_G + "\n")
    filas = ["user_genome,classification"]
    for p, fam in (("a", "Fa"), ("b", "Fb"), ("c", "Fc")):
        for k in range(1, 5):
            filas.append(f"{p}{k},d__Archaea;p__P1;c__C1;o__O1;f__{fam};g__G{p};s__G{p} sp{k}")
    (tmp_path / "tabla.csv").write_text("\n".join(filas) + "\n")
    return tmp_path
