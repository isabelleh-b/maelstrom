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

from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.5, 5.8), dpi=200,
                                 gridspec_kw={'height_ratios': [1.3, 1.6]})

# ---------------- Panel 1: Maelstrom domain map (Homo sapiens, Q96JY0, 434 aa) ----------------
total_len = 434
hmg_start, hmg_end = 5, 65          # matches the HMG_box_human selection in the PyMOL session
msd_start, msd_end = 130, 325       # InterPro-supported MSD boundary
amb_start, amb_end = 106, 130       # AlphaFold/PyMOL model treats this as MSD; InterPro does not

ax1.plot([1, total_len], [0, 0], color='#888888', lw=2, zorder=1)
ax1.add_patch(Rectangle((hmg_start, -0.18), hmg_end - hmg_start, 0.36,
               facecolor='#A9D6F5', edgecolor='#4a90c2', lw=1.3, zorder=2))
ax1.text((hmg_start + hmg_end) / 2, 0, 'HMG', ha='center', va='center',
          fontsize=10.5, fontweight='bold', zorder=3)
ax1.add_patch(Rectangle((amb_start, -0.18), amb_end - amb_start, 0.36,
               facecolor='#3F9B3F', edgecolor='#2E7D32', lw=1.3, zorder=2))
ax1.add_patch(Rectangle((msd_start, -0.18), msd_end - msd_start, 0.36,
               facecolor='#8CE28C', edgecolor='#2E7D32', lw=1.3, zorder=2))
ax1.text((msd_start + msd_end) / 2, 0, 'MAEL-specific domain (MSD)', ha='center', va='center',
          fontsize=10.5, fontweight='bold', zorder=3)

# DEDD positions are labelled by ANCESTRAL residue identity, the convention used
# throughout the dissertation. H. sapiens carries substitutions at all three
# (Asn119, Phe121, Val298); this is stated explicitly in the figure caption.
dedd_pos = {'Asp119': 119, 'Glu121': 121, 'Asp298': 298}
echc_pos = {'Glu138': 138, 'Cys283': 283, 'His286': 286, 'Cys294': 294}

# stagger levels chosen by hand so labels in the crowded 283-294 cluster don't collide
dedd_level = {'Asp119': 0, 'Glu121': 1, 'Asp298': 0}
echc_level = {'Glu138': 0, 'Cys283': 0, 'His286': 1, 'Cys294': 2}

for label, pos in dedd_pos.items():
    lvl = dedd_level[label]
    top = 0.62 + 0.30 * lvl
    ax1.plot([pos, pos], [0.18, top], color='#C62828', lw=1.2, zorder=4)
    ax1.scatter([pos], [top], color='#C62828', s=40, zorder=5, edgecolor='black', lw=0.5)
    xoff = {'Asp119': -7, 'Glu121': 0, 'Asp298': 0}.get(label, 0)
    ha = 'right' if xoff < 0 else 'center'
    ax1.text(pos + xoff, top + 0.08, label, ha=ha, va='bottom', fontsize=8.3, color='#C62828')

for label, pos in echc_pos.items():
    lvl = echc_level[label]
    bot = -0.18 - 0.34 * (lvl + 1)
    ax1.plot([pos, pos], [-0.18, bot], color='#8a7a00', lw=1.2, zorder=4)
    ax1.scatter([pos], [bot], color='#FFD54F', s=40, zorder=5, edgecolor='#8a7a00', lw=0.7)
    ax1.text(pos, bot - 0.08, label, ha='center', va='top', fontsize=8.3, color='#8a7a00')

for tick in [1, 100, 200, 300, 400, total_len]:
    ax1.plot([tick, tick], [-0.02, 0.02], color='#888888', lw=1)
    ax1.text(tick, -1.75, str(tick), ha='center', fontsize=8.3, color='#555555')

ax1.text(1, 1.35, 'DEDD (absent catalytic tetrad; substituted in H. sapiens, see legend)',
          color='#C62828', fontsize=9.3, va='bottom')
ax1.text(1, -1.9, 'ECHC (Zn-coordinating structural motif; retained across taxa)', color='#8a7a00',
          fontsize=9.3, va='top')

# key for the two MSD shades
key_y = -2.45
ax1.add_patch(Rectangle((6, key_y - 0.10), 16, 0.20, facecolor='#3F9B3F',
               edgecolor='#2E7D32', lw=0.9, zorder=3))
ax1.text(28, key_y, 'residues 106\u2013130: modelled as part of the MSD by AlphaFold/PyMOL, '
          'but not assigned to the domain by InterPro',
          fontsize=8.0, va='center', color='#333333')

ax1.set_xlim(-10, total_len + 10)
ax1.set_ylim(-2.85, 1.75)
ax1.axis('off')
ax1.set_title('Maelstrom domain architecture (Homo sapiens, Q96JY0, 434 aa)\n'
               'HMG-box residues 5\u201365; MSD residues 106/130\u2013325 (boundary used for structural superposition)',
               fontsize=10, pad=2)

# ---------------- Panel 2: active-site comparison schematic (not to scale / not co-linear) ----------------
ax2.axis('off')
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 4)

ax2.add_patch(Rectangle((0.3, 2.5), 4.4, 1.0, facecolor='#F5C6C6', edgecolor='#8B0000', lw=1.3))
ax2.text(2.5, 3.0, 'Active DnaQ-H nuclease\n(e.g. E. coli \u03b5-subunit, DnaQ)', ha='center',
          va='center', fontsize=9.6, fontweight='bold')
dnaq_resids = ['Asp12', 'Glu14', 'Asp103', 'Asp167']
for i, r in enumerate(dnaq_resids):
    x = 0.75 + i * 1.0
    ax2.scatter([x], [2.3], color='#C62828', s=60, edgecolor='black', lw=0.5, zorder=5)
    ax2.text(x, 2.02, r, ha='center', fontsize=8, color='#C62828')
ax2.text(2.5, 1.65, 'catalytic DEDD tetrad, coordinates 2 \u00d7 Mn\u00b2\u207a, hydrolyses RNA/DNA',
          ha='center', fontsize=8.3, style='italic', color='#555555')

ax2.add_patch(Rectangle((5.3, 2.5), 4.4, 1.0, facecolor='#FFF3B0', edgecolor='#8a7a00', lw=1.3))
ax2.text(7.5, 3.0, 'Metazoan Maelstrom\n(pseudonuclease)', ha='center', va='center',
          fontsize=9.6, fontweight='bold')
echc_resids = ['Glu138', 'Cys283', 'His286', 'Cys294']
for i, r in enumerate(echc_resids):
    x = 5.75 + i * 1.0
    ax2.scatter([x], [2.3], color='#FFD54F', s=60, edgecolor='#8a7a00', lw=0.7, zorder=5)
    ax2.text(x, 2.02, r, ha='center', fontsize=8, color='#8a7a00')
ax2.text(7.5, 1.65, 'ECHC motif, coordinates 1 \u00d7 Zn\u00b2\u207a, structural rather than catalytic',
          ha='center', fontsize=8.3, style='italic', color='#555555')

arrow = FancyArrowPatch((4.75, 3.0), (5.25, 3.0), arrowstyle='-|>', mutation_scale=16,
                          color='#444444', lw=1.6)
ax2.add_patch(arrow)
ax2.text(5.0, 3.28, 'active-site switch', ha='center', fontsize=8.6, color='#444444',
          style='italic')

ax2.text(5.0, 0.9, 'Same DnaQ-H structural scaffold; independent taxa retain different subsets of the\n'
                    'ancestral DEDD tetrad (Discussion), but the ECHC motif is conserved at 96\u201399% across all 81 taxa sampled here.',
          ha='center', fontsize=8.6, color='#333333')

fig.tight_layout()
save(fig, 'Figure1_domain_map_activesite')
plt.close(fig)
