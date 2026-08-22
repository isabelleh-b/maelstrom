import csv
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

from matplotlib.patches import Circle, Rectangle
from matplotlib.lines import Line2D

from newick import parse_newick, compute_layout, ladderize, prune_leaf, rescale_y
import tree_plot as tp

# ---------- load conservation data ----------
with open(DATA / 'conservation_full.csv') as f:
    rows = list(csv.reader(f))
header = rows[0]
data = rows[1:]

def key(name):
    return name.replace(' ', '_')

# Look up column positions by header name rather than hardcoding indices,
# to avoid off-by-one errors if the sheet layout shifts.
idx_group = header.index('Taxonomic Group')
idx_clade = header.index('Clade/Order')
idx_d119 = header.index('DEDD Asp119(anc.)')
idx_d119m = idx_d119 + 1  # 'match' column immediately follows each residue column
idx_e121 = header.index('DEDD Glu121(anc.)')
idx_e121m = idx_e121 + 1
idx_d298 = header.index('DEDD Asp298(anc.)')
idx_d298m = idx_d298 + 1
assert header[idx_d119m] == 'match' and header[idx_e121m] == 'match' and header[idx_d298m] == 'match'

cons = {}
for r in data:
    name = key(r[0])
    cons[name] = {
        'group': r[idx_group],
        'clade': r[idx_clade],
        'DEDD_Asp119': r[idx_d119], 'DEDD_Asp119_match': r[idx_d119m],
        'DEDD_Glu121': r[idx_e121], 'DEDD_Glu121_match': r[idx_e121m],
        'DEDD_Asp298': r[idx_d298], 'DEDD_Asp298_match': r[idx_d298m],
    }

groups = sorted(set(r[1] for r in data))
cmap = plt.get_cmap('tab20')
group_color = {g: cmap(i / max(1, len(groups) - 1) * 0.95) for i, g in enumerate(groups)}
# make Amoebozoa a strong, unambiguous dark red rather than a tab20 sample
group_color['Amoebozoa/Outgroup'] = '#B22222'

def dedd_state(name, pos):
    c = cons[name]
    match = c[f'DEDD_{pos}_match']
    resid = c[f'DEDD_{pos}']
    if resid == '-':
        return 'gap'
    if match == 'Y':
        return 'retained'
    return 'substituted'

STATE_COLOR = {'retained': '#2E7D32', 'substituted': '#BDBDBD', 'gap': '#FFFFFF'}
STATE_EDGE = {'retained': '#1B5E20', 'substituted': '#757575', 'gap': '#9E9E9E'}
DEDD_POS_LABELS = {'Asp119': 'D119', 'Glu121': 'E121', 'Asp298': 'D298'}
DEDD_POS_ORDER = ['Asp119', 'Glu121', 'Asp298']


def draw_dedd_columns(ax, leaves, x0, col_width, header_y_offset=1.2, marker_r=0.32,
                       label_fontsize=8, header_fontsize=9):
    for j, pos in enumerate(DEDD_POS_ORDER):
        cx = x0 + j * col_width
        ax.text(cx, leaves[0].y - header_y_offset, DEDD_POS_LABELS[pos],
                 fontsize=header_fontsize, ha='center', va='bottom', rotation=0, fontweight='bold')
        for leaf in leaves:
            state = dedd_state(leaf.name, pos)
            hatch = '////' if state == 'gap' else None
            circ = Circle((cx, leaf.y), marker_r, facecolor=STATE_COLOR[state],
                           edgecolor=STATE_EDGE[state], lw=0.9, zorder=4, hatch=hatch)
            ax.add_patch(circ)
    return x0 + (len(DEDD_POS_ORDER) - 1) * col_width


# ========== FIGURE 7: domain-only tree, all taxa, with DEDD status ==========
with open(DATA / 'domain_only_tree_REROOTED.treefile') as f:
    root2 = parse_newick(f.read())
root2 = prune_leaf(root2, 'Macaca_mulatta')
ladderize(root2)
leaves2 = compute_layout(root2)
n2 = len(leaves2)
print('Figure 7 (domain-only tree), taxa:', n2)

max_x2 = max(l.x for l in leaves2)

fig_h = 0.235 * n2 + 1.6
fig = plt.figure(figsize=(11.5, fig_h), dpi=200)
ax = fig.add_axes([0.02, 0.03, 0.60, 0.93])
ax.set_facecolor('white')

def branch_color_2(node):
    tips = node.leaves()
    grp = set(cons[t.name]['group'] for t in tips)
    if len(grp) == 1:
        return group_color[grp.pop()]
    return '#8a8a8a'

tp.draw_tree(ax, root2, leaves2, tip_color=lambda l: group_color[cons[l.name]['group']],
             branch_color=branch_color_2, label_fontsize=8.6, linewidth=1.1,
             label_xpad=0.04 * max_x2, node_support_fontsize=None)

# DEDD columns
label_width_guess = 0.62 * max_x2  # generous room for longest binomial name
x0 = max_x2 + label_width_guess
col_w = 0.14 * max_x2
last_x = draw_dedd_columns(ax, leaves2, x0, col_w, header_y_offset=2.2,
                            marker_r=0.30, header_fontsize=9)

# scale bar
import math
raw = max_x2 * 0.25
mag = 10 ** math.floor(math.log10(raw))
nice_opts = [1, 2, 5, 10]
scale_len = min(nice_opts, key=lambda m: abs(m * mag - raw)) * mag
unit_word2 = 'substitution' if scale_len == 1 else 'substitutions'
tp.scale_bar(ax, 0, n2 + 3.0, scale_len, f'{scale_len:g} {unit_word2}/site', fontsize=9)

ax.set_xlim(-0.02 * max_x2, last_x + col_w * 1.3)
ax.set_ylim(n2 + 5.0, -3.0)
ax.axis('off')

# legend: taxonomic groups, placed in dedicated panel to the right, not overlapping tree
lax = fig.add_axes([0.66, 0.03, 0.33, 0.93])
lax.axis('off')
lax.set_xlim(0, 1)
lax.set_ylim(0, 1)
lax.text(0.0, 0.995, 'Taxonomic group', fontsize=11, fontweight='bold', va='top')
n_groups = len(groups)
taxo_top, taxo_bottom = 0.96, 0.46  # reserve the lower ~46% of the panel for the DEDD legend
row_h = (taxo_top - taxo_bottom) / n_groups
for i, g in enumerate(groups):
    y = taxo_top - i * row_h
    lax.add_patch(Rectangle((0.0, y - 0.012), 0.045, 0.024, facecolor=group_color[g],
                             edgecolor='none'))
    lax.text(0.065, y, g, fontsize=9.3, va='center')

y0 = taxo_bottom - 0.05
lax.text(0.0, y0, 'DEDD motif (ancestral identity, H. sapiens numbering)', fontsize=10.3,
          fontweight='bold', va='top')
states = ['retained', 'substituted', 'gap']
state_labels = {'retained': 'Ancestral residue retained', 'substituted': 'Substituted',
                'gap': 'Alignment gap'}
for i, st in enumerate(states):
    y = y0 - 0.07 - i * 0.075
    circ = Circle((0.022, y), 0.018, transform=lax.transData, facecolor=STATE_COLOR[st],
                   edgecolor=STATE_EDGE[st], lw=1.0, hatch='////' if st == 'gap' else None)
    lax.add_patch(circ)
    lax.text(0.065, y, state_labels[st], fontsize=9.3, va='center')

save(fig, 'Figure7_domain_tree_DEDD')
plt.close(fig)

# ================= FIGURE 8: Amoebozoa clade zoom =================
with open(DATA / 'domain_only_tree_REROOTED.treefile') as f:
    root_full = parse_newick(f.read())
amz_names = [n for n, c in cons.items() if c['group'] == 'Amoebozoa/Outgroup']
clade = root_full.get_common_ancestor(amz_names)
clade.parent = None  # detach: this subtree is its own root for layout purposes
ladderize(clade)
leaves4 = compute_layout(clade)
rescale_y(clade, 1.6)
n4 = len(leaves4)
print('Figure 8 (Amoebozoa clade), taxa:', n4)
max_x4 = max(l.x for l in leaves4)
max_y4 = max(l.y for l in leaves4)

fig2 = plt.figure(figsize=(11.0, 0.62 * max_y4 + 2.0), dpi=200)
ax2 = fig2.add_axes([0.03, 0.06, 0.60, 0.86])

tp.draw_tree(ax2, clade, leaves4, tip_color=lambda l: '#B22222',
             branch_color=lambda n: '#B22222', label_fontsize=12, linewidth=1.6,
             label_xpad=0.05 * max_x4, node_support_fontsize=None)

# Manually place internal-node support labels with simple collision avoidance:
# nodes whose y-values fall close together get staggered onto different
# horizontal offsets so the text doesn't overlap.
support_nodes = [n for n in clade.iter_nodes() if (not n.is_leaf()) and n.support]
support_nodes.append(clade)  # include the clade root's own split support
support_nodes = list({id(n): n for n in support_nodes}.values())
support_nodes.sort(key=lambda n: n.y)

placed = []  # list of (y, offset_level)
y_thresh = 0.95
base_offset = 0.055 * max_x4
for n in support_nodes:
    level = 0
    for (py, plevel) in placed:
        if abs(py - n.y) < y_thresh:
            level = max(level, plevel + 1)
    placed.append((n.y, level))
    xoff = base_offset * (1 + level)
    yjit = 0.0 if level == 0 else (0.32 if level % 2 else -0.32)
    ax2.text(n.x - xoff, n.y - 0.30 + yjit, n.support, fontsize=10, color='#555555',
              ha='right', va='center',
              bbox=dict(boxstyle='round,pad=0.06', fc='white', ec='none', alpha=0.9), zorder=5)

label_width_guess4 = 1.25 * max_x4
x0b = max_x4 + label_width_guess4
col_w4 = 0.22 * max_x4
last_x4 = draw_dedd_columns(ax2, leaves4, x0b, col_w4, header_y_offset=0.9,
                             marker_r=0.28, header_fontsize=11, label_fontsize=10)

raw4 = max_x4 * 0.3
mag4 = 10 ** math.floor(math.log10(raw4))
scale_len4 = min(nice_opts, key=lambda m: abs(m * mag4 - raw4)) * mag4
unit_word = 'substitution' if scale_len4 == 1 else 'substitutions'
tp.scale_bar(ax2, 0, max_y4 + 1.6, scale_len4, f'{scale_len4:g} {unit_word}/site', fontsize=10.5)

ax2.set_xlim(-0.02 * max_x4, last_x4 + col_w4 * 1.3)
ax2.set_ylim(max_y4 + 2.6, -2.2)
ax2.axis('off')
ax2.text((0 + max_x4) / 2, -1.4, 'SH-aLRT / UFBoot support shown at each node',
          fontsize=10, ha='center', style='italic', color='#444444')

lax2 = fig2.add_axes([0.67, 0.06, 0.31, 0.86])
lax2.axis('off')
lax2.set_xlim(0, 1); lax2.set_ylim(0, 1)
lax2.text(0.0, 0.98, 'DEDD motif\n(ancestral identity,\nH. sapiens numbering)', fontsize=11,
           fontweight='bold', va='top')
for i, st in enumerate(states):
    y = 0.72 - i * 0.11
    circ = Circle((0.03, y), 0.028, transform=lax2.transData, facecolor=STATE_COLOR[st],
                   edgecolor=STATE_EDGE[st], lw=1.2, hatch='////' if st == 'gap' else None)
    lax2.add_patch(circ)
    lax2.text(0.1, y, state_labels[st], fontsize=10.5, va='center')

save(fig2, 'Figure8_amoebozoa_clade')
plt.close(fig2)
