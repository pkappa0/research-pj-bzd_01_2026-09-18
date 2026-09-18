# Strict 3-drug PoC QC report

This report is an audit trail. No values are averaged, no missing interactions are imputed, and an absent interaction is distinct from absent interaction data.

- ChEMBL comparable Ki pair rows: **122** (pair-level; `pair_priority=1` means same document, type and unit).
- Docking pose rows: **54**.
- PLIP interaction rows: **238**.
- Residue mapping rows requiring QC: **0**.

## Structural QC

6HUP is an exact human α1β3γ2L receptor with a diazepam copy assigned to alpha1 chain D at the alpha1/gamma2 site; its coordinates define the alpha1 box. 7QNE is an independent exact α1β3γ2 structure with Ro15-4513 and is retained as an alternative candidate.

No exact human α2β3γ2 full pentamer with a BZD-site co-crystal ligand was found in the investigated RCSB set. The alpha2 docking receptor is therefore a provisional local-site construct consisting of chains C/D/E (β3/α2/γ2) extracted from mixed native 9CTJ. It is not equivalent to an exact α2β3γ2 pentamer. The DZP-derived box is transferred by a 207-residue gamma2 C→E Cα superposition (RMSD recorded in `docking_boxes.csv`); this transfer is reproducible but biologically uncertain.

The strict mechanistic interpretation must therefore treat alpha2 docking, residue mapping, and any IFP comparison as exploratory/provisional. They do not establish subtype-specific causality.

## Docking QC

All six jobs use the same Vina seed, exhaustiveness, number of modes, energy range, and ligand preparation protocol. PDBQT conversion uses a transparent local writer with RDKit Gasteiger charges and element/aromatic atom types; receptor partial charges are zero because Open Babel/Meeko were not available. This is a methodological limitation.

## Interaction QC

PLIP is run from the vendored package when available. Interactions are retained pose-by-pose. Features that cannot be mapped to a common position are excluded from the fingerprint tables and remain visible as `qc_required` in `residue_mapping.csv`; they are not auto-resolved.
