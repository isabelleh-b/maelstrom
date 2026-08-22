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

from matplotlib.patches import Rectangle

seqs = {}
name = None
with open(DATA / 'domain_only_aln.fasta') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split('/')[0]
            seqs[name] = ''
        else:
            seqs[name] += line

# Verified column positions (Homo sapiens / H. sapiens numbering), domain-only
# alignment coordinates (= full-length column - 424), cross-checked directly
# against known residues in Methods above.
POS = {'D119': 16, 'E121': 18, 'Glu138': 36, 'Cys283': 233, 'His286': 236,
       'Cys294': 248, 'Asp298': 252}
DISPLAY_LABEL = {'D119': 'D119', 'E121': 'E121', 'Glu138': 'E138 (ECHC)',
                  'Cys283': 'C283 (ECHC)', 'His286': 'H286 (ECHC)',
                  'Cys294': 'C294 (ECHC)', 'Asp298': 'D298'}

reps = ['Homo_sapiens', 'Mus_musculus', 'Gallus_gallus', 'Taeniopygia_guttata',
        'Corvus_moneduloides', 'Danaus_plexippus', 'Strongylocentrotus_purpuratus',
        'Ciona_savignyi']
amz = ['Acanthamoeba_castellani', 'Rhizamoeba_saxonica', 'Flabellula_baltica',
       'Dictyostelium_discoideum', 'Soliformovum_irregulare', 'Ceratiomyxa_fruticulosa',
       'Filamoeba_P010183', 'Filamoeba_P014044', 'Dermamoeba_algensis']
taxa = reps + amz

# Build three windows with 4-column padding either side of the tracked positions
def make_window(cols, pad=4):
    lo = min(cols) - pad
    hi = max(cols) + pad
    return lo, hi

winA_lo, winA_hi = make_window([POS['D119'], POS['E121']])
winB_lo, winB_hi = make_window([POS['Glu138']])
winC_lo, winC_hi = make_window([POS['Cys283'], POS['His286'], POS['Cys294'], POS['Asp298']])
windows = [('Window A', winA_lo, winA_hi), ('Window B', winB_lo, winB_hi),
           ('Window C', winC_lo, winC_hi)]

tracked_cols = set(POS.values())
gap_between_windows = 2.2

# x-position lookup: assign each alignment column (within a window) an x slot
col_x = {}
x_cursor = 0.0
window_spans = []
for wname, lo, hi in windows:
    start_x = x_cursor
    for c in range(lo, hi + 1):
        col_x[c] = x_cursor
        x_cursor += 1.0
    window_spans.append((wname, start_x, x_cursor - 1.0, lo, hi))
    x_cursor += gap_between_windows

n_taxa = len(taxa)
row_gap_after = 8  # extra visual gap between metazoan block and amoebozoa block

def row_y(i):
    return i if i < 8 else i + 1  # insert one blank row's worth of gap after index 7

fig_w = x_cursor + 4.5
fig_h = 0.5 * (n_taxa + 1) + 2.0
fig = plt.figure(figsize=(fig_w, fig_h), dpi=200)
ax = fig.add_axes([0.24, 0.06, 0.74, 0.86])

cell_w, cell_h = 0.92, 0.82

for i, tx in enumerate(taxa):
    y = row_y(i)
    seq = seqs[tx]
    label = tx.replace('_', ' ')
    ax.text(-0.6, y, label, fontsize=11, style='italic', ha='right', va='center')
    for wname, sx, ex, lo, hi in window_spans:
        for c in range(lo, hi + 1):
            x = col_x[c]
            resid = seq[c - 1] if 0 <= c - 1 < len(seq) else '-'
            is_tracked = c in tracked_cols
            if resid == '-':
                face = '#FFFFFF'
                hatch = '////' if is_tracked else None
                edge = '#9E9E9E'
            elif is_tracked:
                face = '#FFF2A8'  # highlight tracked column
                hatch = None
                edge = '#8a7a00'
            else:
                face = '#FFFFFF'
                hatch = None
                edge = '#DDDDDD'
            rect = Rectangle((x - cell_w / 2, y - cell_h / 2), cell_w, cell_h,
                              facecolor=face, edgecolor=edge, lw=1.1 if is_tracked else 0.6,
                              hatch=hatch, zorder=2)
            ax.add_patch(rect)
            ax.text(x, y, resid, fontsize=10.5, ha='center', va='center',
                     family='monospace', fontweight='bold' if is_tracked else 'normal',
                     zorder=3)

# column headers: tracked-position labels above their column, per window
for label, c in POS.items():
    x = col_x[c]
    ax.text(x, -1.5, DISPLAY_LABEL[label], fontsize=10, ha='left', va='bottom',
             fontweight='bold', rotation=40)

# window group labels
for wname, sx, ex, lo, hi in window_spans:
    ax.text((sx + ex) / 2, -3.6, wname, fontsize=11, ha='center', va='top', fontweight='bold')
    ax.plot([sx - 0.5, ex + 0.5], [-3.2, -3.2], color='#888888', lw=1.0)

# divider + group labels between Metazoa and Amoebozoa blocks
div_y = row_y(7) + 1.0
ax.plot([col_x[windows[0][1]] - 0.6, x_cursor - gap_between_windows], [div_y, div_y],
         color='#888888', lw=0.8, linestyle=(0, (4, 3)))
ax.text(-3.4, 3.5, 'Metazoa', fontsize=11.5, fontweight='bold', ha='center', va='center',
         rotation=90)
ax.text(-3.4, row_y(8) + 4.0, 'Amoebozoa /\nOutgroup', fontsize=11.5, fontweight='bold',
         ha='center', va='center', rotation=90)

ax.set_xlim(-6.5, x_cursor)
ax.set_ylim(row_y(n_taxa - 1) + 1.2, -4.2)
ax.axis('off')

# legend
lx0 = x_cursor - gap_between_windows + 0.6
ax.add_patch(Rectangle((lx0, row_y(0) - 0.4), 0.9, 0.8, facecolor='#FFF2A8', edgecolor='#8a7a00', lw=1.1))
ax.text(lx0 + 1.1, row_y(0), 'Tracked ECHC / DEDD column', fontsize=10, va='center')
ax.add_patch(Rectangle((lx0, row_y(1) - 0.4), 0.9, 0.8, facecolor='#FFFFFF', edgecolor='#DDDDDD', lw=0.6))
ax.text(lx0 + 1.1, row_y(1), 'Flanking alignment column (context)', fontsize=10, va='center')
ax.add_patch(Rectangle((lx0, row_y(2) - 0.4), 0.9, 0.8, facecolor='#FFFFFF', edgecolor='#9E9E9E', lw=0.6, hatch='////'))
ax.text(lx0 + 1.1, row_y(2), 'Alignment gap', fontsize=10, va='center')

save(fig, 'Figure9_active_site_alignment')
plt.close(fig)
