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


# Data taken directly from the verified draft text / Table 1 (fungal Foldseek
# candidates) and the AlphaFold superposition results (Amorphea comparisons).
amorphea = [
    ('A. castellanii vs H. sapiens\n(AlphaFold superposition)', 0.821),
    ('D. melanogaster vs H. sapiens\n(AlphaFold superposition)', 1.240),
]
fungal = [
    ('Paramarasmius palmivorus', 3.329),
    ('Mycena alexandri', 4.352),
    ('Cladophialophora carrionii', 5.563),
    ('Piromyces finnis', 8.545),
    ('Sporothrix schenckii', 10.897),
    ('Schizosaccharomyces pombe', 11.432),
    ('Saccharomyces cerevisiae', 12.866),
    ('Cytospora leucostoma', 17.095),
]
fungal_sorted = sorted(fungal, key=lambda x: x[1])

rows = amorphea + [('', None)] + fungal_sorted
labels = [r[0] for r in rows]
values = [r[1] for r in rows]

fig, ax = plt.subplots(figsize=(8.5, 5.2), dpi=200)
y = list(range(len(rows)))[::-1]

colors = ['#2E7D32'] * len(amorphea) + ['none'] + ['#B22222'] * len(fungal_sorted)

for yi, val, col in zip(y, values, colors):
    if val is None:
        continue
    ax.plot([0, val], [yi, yi], color=col, lw=2.2, zorder=1, alpha=0.55)
    ax.scatter([val], [yi], color=col, s=90, zorder=3, edgecolor='black', linewidth=0.6)
    ax.text(val + 0.35, yi, f'{val:.3f} \u00c5', va='center', fontsize=9.5)

ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=9.5)
ax.set_xlabel('RMSD (\u00c5) to Bombyx mori / query MAEL domain', fontsize=10.5)
ax.set_xlim(0, max(v for v in values if v is not None) * 1.22)
ax.spines[['top', 'right']].set_visible(False)
ax.axvspan(0, 1.5, color='#2E7D32', alpha=0.06, zorder=0)
ax.text(0.05, len(rows) - 0.3, 'Amorphea (structural homology confirmed)',
         fontsize=9, color='#2E7D32', fontweight='bold', va='bottom')
ax.text(0.05, len(fungal_sorted) - 0.4, 'Fungal Foldseek candidates (no confirmed homology)',
         fontsize=9, color='#B22222', fontweight='bold', va='bottom')

fig.tight_layout()
save(fig, 'Figure6_RMSD_comparison')
plt.close(fig)
