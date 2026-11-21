"""Tree geometry: x/y positions and triangles for collapsed clades."""
import math

from .tree import preorden


def layout_arbol(raiz, colapsados, modo, lado, x_raiz, w, esc_tri):
    """Compute tree positions (with triangles for collapsed nodes)."""
    y, y0, y1, xs = {}, {}, {}, {}
    unidades, cursor = [], [0.0]
    max_alt = max(raiz.altura, 1)
    max_prof = max(raiz.profmax, 1e-12)

    def fx(frac):
        return x_raiz + frac * w if lado == "izq" else x_raiz - frac * w

    def frac_nodo(n):
        return 1.0 - n.altura / max_alt if modo == "cladograma" else n.prof / max_prof

    def frac_base(n):
        return 1.0 if modo == "cladograma" else n.profmax / max_prof

    def alto_tri(n):
        return esc_tri * min(5.0, 1.0 + 0.4 * math.log2(max(n.n, 2)))

    def rec(n):
        if n.es_hoja or n in colapsados:
            h = 1.0 if n.es_hoja else alto_tri(n)
            yc = cursor[0] + h / 2
            cursor[0] += h
            y[n], y0[n], y1[n] = yc, yc - h / 2, yc + h / 2
            xs[n] = fx(frac_nodo(n))
            unidades.append({"nodo": n, "tri": not n.es_hoja, "y": yc, "h": h,
                             "x_apice": xs[n], "x_base": fx(frac_base(n)), "hojas": n.hojas})
            return
        for h in n.hijos:
            rec(h)
        y[n] = (y[n.hijos[0]] + y[n.hijos[-1]]) / 2
        y0[n], y1[n] = y0[n.hijos[0]], y1[n.hijos[-1]]
        xs[n] = fx(frac_nodo(n))

    rec(raiz)
    segmentos = []
    for n in preorden(raiz):
        if n.es_hoja or n in colapsados or n not in xs:
            continue
        if n.padre is not None and n.padre in colapsados:
            continue
        for h in n.hijos:
            segmentos.append((xs[n], y[n], xs[n], y[h], h))
            segmentos.append((xs[n], y[h], xs[h], y[h], h))
    return {"unidades": unidades, "segmentos": segmentos, "y": y, "y0": y0, "y1": y1,
            "x": xs, "filas": cursor[0]}


def escalar_y(lay, factor):
    for u in lay["unidades"]:
        u["y"] *= factor
        u["h"] *= factor
    lay["segmentos"] = [(a, b * factor, c, d * factor, n) for a, b, c, d, n in lay["segmentos"]]
    for k in ("y", "y0", "y1"):
        lay[k] = {n: v * factor for n, v in lay[k].items()}
