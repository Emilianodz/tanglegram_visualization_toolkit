"""Tools to compare phylogenetic trees (VAMPhyRE vs GTDB-Tk) and draw tanglegrams
with rogue-branch analysis, plus single-tree figures and iTOL export.

Modules reusable outside the tanglegram:
  names     genome identifier cleaning
  taxonomy  GTDB taxonomy parsing, taxon selection, colour scheme
  idmap     auxiliary table loading and ID harmonisation across files
  tree      Newick I/O and tree operations
  splits    Robinson-Foulds distance and per-branch discordance
  itol      iTOL tree and annotation export
"""
__version__ = "1.0.0"
