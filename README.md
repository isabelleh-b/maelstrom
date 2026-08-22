# maelstrom

Evolutionary analysis of the piRNA-pathway pseudonuclease Maelstrom, whose MAEL-specific domain retains a DnaQ-H nuclease fold while having substituted the catalytic DEDD tetrad for a structural zinc-binding ECHC motif.

This repository forms the supplementary data and analysis for the dissertation:

*Lost Catalysis: Evolution of an Active-Site Switch in the Pseudoenzyme Maelstrom*; student exam number B301434, supervised by A.G. Cook and E.W.J. Wallace. MSc Biotechnology, School of Biological Sciences, University of Edinburgh, 2025–26.

The analysis covers 81 taxa: 72 Metazoa and a nine-taxon Amoebozoa outgroup. It tests whether the DEDD-to-ECHC active-site switch occurred once at the base of Metazoa or repeatedly across lineages, combining multiple sequence alignment, maximum-likelihood phylogenetics, per-residue motif conservation, and structural superposition. All motif positions are reported in *Homo sapiens* full-length numbering (UniProt Q96JY0): ECHC at Glu138, Cys283, His286 and Cys294, and the ancestral DEDD identities at Asp119, Glu121 and Asp298.

# Contents

## data

Data files, both input (retrieved sequences) and output (aligned sequences, trees, motif conservation summaries).

### data/fasta

* `protein-sequencesbeforealignment.fasta` — candidate sequences retrieved from UniProt before domain-architecture filtering and taxon selection
* `full_length_aln.fasta` — full-length alignment, MAFFT L-INS-i, 81 sequences, 1288 columns
* `domain_only_aln.fasta` — domain-only alignment, columns 425–760 of the full-length alignment sliced to the MAEL-specific domain boundary, 81 sequences, 336 columns

### data/treefile

Maximum-likelihood trees from IQ-TREE3, rerooted on the monophyletic Amoebozoa outgroup. Internal node labels give support as `SH-aLRT/UFBoot`.

* `domain_only_tree_REROOTED.treefile` — inferred under LG+I+G4, used as the primary phylogeny throughout, since it recovers the Amoebozoa clade with strong support (99.9/99)
* `full_length_tree_REROOTED.treefile` — inferred under JTT+R5, used only for the preliminary HMG-box annotation, since it places the outgroup weakly (14.2/70)

### data/csv

* `Full_data_81_taxa.csv` — per-taxon residue calls at all seven tracked positions, with taxonomic group and clade or order, and a match column against the expected identity
* `ECHC_DEDD_conservation_final.csv` — summary of the same data by position: alignment column, expected residue, taxa matching, and conservation percentage

## analysis

Python scripts which generate figures and tables from the files in `data`.

* `make_Figure1_domain_map_activesite.py` — domain map and active-site switch schematic
* `make_Figure4_methods_schematic.py` — analysis pipeline overview
* `make_Figure6_RMSD_comparison.py` — RMSD of candidate structural matches to the *Bombyx mori* query
* `make_Figure8_trees.py` — domain-only phylogeny with DEDD status, and the Amoebozoa clade at higher resolution with support values
* `make_Figure9_active_site_alignment.py` — alignment columns at the tracked ECHC and DEDD positions
* `make_rdm_tables.py` — LaTeX accession tables for the Research Data Management section and Appendix A

Alignment and phylogenetic inference were run directly at the command line:

```
mafft --localpair --maxiterate 1000 full_length_sequences.fasta > full_length_aln.fasta
iqtree3 -s [alignment] -m MFP --alrt 1000 -B 1000 -T AUTO
```

Both treefiles were then rerooted on the Amoebozoa outgroup with a custom Python Newick script.

## figure

Figures used in the dissertation, numbered to match. Figures 1, 4, 6, 7, 8 and 9 are output by the scripts in `analysis`. Figures 2 and 5 are output by PyMOL, and Figures 10 and 11 by iTOL.

* `Figure1_domain_map_activesite.pdf`
* `Figure2_Homo_MAEL_domain_colouring.png`
* `Figure3_hypothesis_sketch.pdf`
* `Figure4_methods_schematic.pdf`
* `Figure5_alphafold_superposition.png`
* `Figure6_RMSD_comparison.pdf`
* `Figure7_domain_tree_DEDD.pdf`
* `Figure8_amoebozoa_clade.pdf`
* `Figure9_active_site_alignment.pdf`
* `Figure10_Full_sequence_HMG_predictions.pdf`
* `Figure11_Domain_Only_DEDD_ECHC_motif.pdf`

### figure/pymol

* `Figure2HomoSapien_DomainColourCode.pse` — PyMOL session underlying Figure 2
* `align/Alignment code human fly acanthamoeba` — commands for the three-way superposition in Figure 5, with *H. sapiens* 130–326, *D. melanogaster* 120–335 and *A. castellanii* 256–337 superposed onto the human domain
* `align/fungi align or super` — template commands used to superpose each filtered Foldseek fungal candidate onto the human reference, with accession and residue range substituted per candidate

# License

Original code and data analysis in this repository is covered by an Apache 2.0 licence, see LICENCE.

Data copied from other sources are covered by their own licences. Sequence data, accessions and taxonomic information were retrieved from UniProt, OrthoDB, PANTHER, NCBI and EukProt. Structural models were retrieved from AlphaFold DB.
