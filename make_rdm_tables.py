"""
Generate the LaTeX accession tables from the taxon list.

Run prep_data.py first. Produces two files:

    RDM_table_amoebozoa.tex -- Table 1.1, the nine-taxon Amoebozoa outgroup,
                               included inline in the Research Data Management section
    RDM_table_full.tex      -- Table A.1, all 81 taxa, included in Appendix A

Both are longtable environments, so \\usepackage{longtable} is required in the preamble.
"""

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / 'data'
TABLES = HERE.parent / 'tables'   # LaTeX tables for the dissertation

IN_CSV = DATA / 'sequence_list_full.csv'


def esc(s):
    """Escape the LaTeX special characters which occur in this dataset."""
    if s is None:
        return ''
    return (str(s).replace('&', r'\&').replace('%', r'\%')
                  .replace('_', r'\_').replace('#', r'\#'))


def short_source(notes):
    """Collapse the long source filenames into a short, readable label."""
    n = (notes or '').strip()
    if n in ('reviewed', 'unreviewed'):
        return f'UniProt ({n})'
    if 'amoebic_maelstrom_full_sequences' in n:
        return 'EukProt (custom)'
    return esc(n) if n else '--'


def uniprot_cell(code):
    # Amoebic sequences without a UniProt entry carry an em-dash placeholder.
    if code and str(code).startswith('\u2014'):
        return 'none'
    return esc(code) if code else '--'


def load_rows():
    with open(IN_CSV) as f:
        rows = list(csv.reader(f))
    # Drop the header and any row without a species name (the source sheet
    # contains one stray blank Amoebozoa/Discosea row).
    return [r for r in rows[1:] if r[2] and r[2].strip()]


def longtable(colspec, caption, label, header_cells, body_lines):
    out = [f'\\begin{{longtable}}{{{colspec}}}']
    out.append(f'\\caption*{{{caption}}} \\label{{{label}}} \\\\')
    hdr = ' & '.join(f'\\textbf{{{h}}}' for h in header_cells) + r' \\'
    for block in ('\\endfirsthead', '\\endhead'):
        out += [r'\toprule', hdr, r'\midrule', block]
        if block == '\\endfirsthead':
            out.insert(len(out) - 1, '')  # keep the two header blocks visually separate
    out = [ln for ln in out if ln != '']
    # rebuild cleanly: firsthead then head
    out = [f'\\begin{{longtable}}{{{colspec}}}',
           f'\\caption*{{{caption}}} \\label{{{label}}} \\\\',
           r'\toprule', hdr, r'\midrule', r'\endfirsthead',
           r'\toprule', hdr, r'\midrule', r'\endhead',
           r'\bottomrule', r'\endfoot']
    out += body_lines
    out.append(r'\end{longtable}')
    return '\n'.join(out)


def main():
    data = load_rows()
    print(f'{len(data)} taxa loaded')

    # ---- Table 1.1: Amoebozoa outgroup only, no group column needed ----
    amz = [r for r in data if r[0] == 'Amoebozoa/Outgroup']
    body = [
        f'\\textit{{{esc(r[2])}}} & {uniprot_cell(r[3])} & '
        f'{esc(r[4]) if r[4] else "--"} & {short_source(r[5])} \\\\'
        for r in amz
    ]
    tex = longtable(
        '@{} p{4.3cm} p{2.6cm} p{1.4cm} p{2.3cm} @{}',
        f'Table 1.1: \\textit{{Amoebozoa/Outgroup}} accessions ({len(amz)} taxa).',
        'tab:geneids_amoebozoa',
        ['Species', 'UniProt code', 'Length (aa)', 'Source / notes'],
        body)
    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / 'RDM_table_amoebozoa.tex', 'w') as f:
        f.write(tex + '\n')
    print(f'wrote RDM_table_amoebozoa.tex ({len(amz)} taxa)')

    # ---- Table A.1: all taxa, with the taxonomic group column ----
    body = [
        f'{esc(r[0])} & \\textit{{{esc(r[2])}}} & {uniprot_cell(r[3])} & '
        f'{esc(r[4]) if r[4] else "--"} & {short_source(r[5])} \\\\'
        for r in data
    ]
    caption = ('Table A.1: Full sequence accession table. Taxa, UniProt accessions, and '
               'sequence lengths for the 81-taxon dataset used throughout this dissertation '
               '(see Table 1.1, Research Data Management, for the nine-taxon Amoebozoa '
               'outgroup specifically).')
    tex = longtable(
        '@{} p{3.3cm} p{3.9cm} p{2.5cm} p{1.4cm} p{2.1cm} @{}',
        caption, 'tab:geneids_full',
        ['Taxonomic group', 'Species', 'UniProt code', 'Length (aa)', 'Source / notes'],
        body)
    with open(TABLES / 'RDM_table_full.tex', 'w') as f:
        f.write(tex + '\n')
    print(f'wrote RDM_table_full.tex ({len(data)} taxa)')


if __name__ == '__main__':
    main()
