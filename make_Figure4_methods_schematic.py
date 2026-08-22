import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- output/input paths, resolved relative to this file so the script can be
# --- run from anywhere in the repository -------------------------------------
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / 'data'      # input files live in maelstrom/data
FIGURE = HERE.parent / 'figure'  # output figures go to maelstrom/figure

SAVE_PNG = False   # set True to also write a .png alongside each .pdf


def save(fig, stem):
    """Write `fig` to figure/<stem>.pdf, matching the filenames used in the
    dissertation. Set SAVE_PNG above to also emit a preview .png."""
    FIGURE.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE / f'{stem}.pdf', bbox_inches='tight')
    if SAVE_PNG:
        fig.savefig(FIGURE / f'{stem}.png', bbox_inches='tight', dpi=200)
    print(f'wrote {FIGURE / (stem + ".pdf")}')

from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

fig, ax = plt.subplots(figsize=(9.5, 11.5), dpi=200)
ax.set_xlim(0, 10)
ax.set_ylim(0, 15.5)
ax.axis('off')

def box(x, y, w, h, title, body, facecolor, edgecolor, fontsize_title=10.5, fontsize_body=8.7):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.02,rounding_size=0.12',
                                  facecolor=facecolor, edgecolor=edgecolor, lw=1.4, zorder=2))
    ax.text(x + w/2, y + h - 0.28, title, ha='center', va='top', fontsize=fontsize_title,
            fontweight='bold', zorder=3)
    ax.text(x + w/2, y + h - 0.55, body, ha='center', va='top', fontsize=fontsize_body,
            zorder=3, linespacing=1.45)

def arrow(x0, y0, x1, y1, color='#444444', lw=1.6):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=14,
                          color=color, lw=lw, zorder=1)
    ax.add_patch(a)

# ---- Stage 1: dataset assembly ----
box(1.0, 13.3, 8.0, 1.7, '1. Sequence retrieval and dataset assembly',
    'PANTHER \u00b7 OrthoDB \u00b7 UniProt \u00b7 NCBI BLAST (Metazoa)   +   EukProt (Amoebozoa)\n'
    'domain-architecture check (UniProt/InterPro) \u00b7 balancing tree (MEGA12)\n'
    'scope check: EukProt-wide + PhylomeDB (choanoflagellates, other opisthokonts)',
    '#EAF3FF', '#4a90c2')

arrow(5.0, 13.3, 5.0, 12.75)

# ---- Stage 2: similarity search ----
box(1.0, 11.55, 8.0, 1.2, '2. Similarity searches',
    'PSI-BLAST (NCBI nr + EBI UniProtKB), seeded with H. sapiens Q96JY0\n'
    'filtered by the same domain-architecture check as above',
    '#EAF3FF', '#4a90c2')

arrow(5.0, 11.55, 5.0, 11.05)

# ---- 81-taxon dataset marker ----
box(3.0, 9.95, 4.0, 0.95, '', '81-taxon dataset\n(72 Metazoa + 9 Amoebozoa)', '#FFF3B0', '#8a7a00',
    fontsize_body=9.3)

arrow(4.2, 9.95, 2.6, 9.55)
arrow(5.8, 9.95, 7.4, 9.55)

# ---- Stage 3 (left branch): MSA + motifs ----
box(0.6, 7.55, 4.0, 2.0, '3. Alignment and motif ID',
    'MAFFT L-INS-i\n'
    '\u2193\n'
    'full-length alignment (1288 col)\n'
    'domain-only alignment (336 col)\n'
    'ECHC/DEDD verified vs H. sapiens',
    '#E9F7EA', '#2E7D32', fontsize_body=8.3)

arrow(2.6, 7.55, 2.6, 7.05)

box(0.6, 5.5, 4.0, 1.55, '5. Phylogenetic analysis',
    'IQ-TREE3 (ModelFinder, UFBoot2,\nSH-aLRT) on both alignments,\nrerooted on Amoebozoa',
    '#E9F7EA', '#2E7D32', fontsize_body=8.3)

arrow(2.6, 5.5, 2.6, 5.0)

box(0.6, 3.85, 4.0, 1.15, '', 'full-length tree  +  domain-only tree', '#FFF3B0', '#8a7a00',
    fontsize_body=9.0)

# ---- Stage 4 (right branch): structural comparison ----
box(5.4, 7.55, 4.0, 2.0, '4. Structural comparison',
    'AlphaFold DB structures\n'
    '(H. sapiens, D. melanogaster,\nA. castellanii)\n'
    'PyMOL superposition \u00b7 Foldseek search',
    '#FBEAEA', '#B22222', fontsize_body=8.3)

arrow(7.4, 7.55, 7.4, 7.05)

box(5.4, 5.85, 4.0, 1.2, '', 'RMSD comparisons\n(structural homology)', '#FFF3B0', '#8a7a00',
    fontsize_body=9.0)

# ---- convergence into Results ----
arrow(2.6, 3.85, 4.6, 3.15)
arrow(7.4, 5.85, 5.4, 3.15)

box(2.5, 1.7, 5.0, 1.3, 'Results',
    'which tree is primary; motif conservation\nacross taxa; structural homology; HMG-box',
    '#F3F3F3', '#555555', fontsize_body=8.8)

fig.tight_layout()
save(fig, 'Figure4_methods_schematic')
