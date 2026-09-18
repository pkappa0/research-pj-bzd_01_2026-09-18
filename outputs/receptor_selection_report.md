# Receptor selection report

## Scope

The strict PoC compares diazepam, triazolam and zolpidem against ChEMBL targets CHEMBL2094121 (human α1β3γ2) and CHEMBL2094130 (human α2β3γ2). Existing ChEMBL raw activity is kept unchanged under `data/raw/pharmacology/chembl/`.

## Candidate structures

```
pdb_id                                                                         title       composition co_crystal_ligand  resolution_angstrom                                                bzd_site_assessment bzd_ligand_present                      poc_decision
  6HUP    human full-length alpha1beta3gamma2L GABA(A)R with diazepam, GABA and Mb38 alpha1beta3gamma2    DZP (diazepam)                 3.58 exact alpha1beta3gamma2 receptor; BZD ligand in alpha1/gamma2 site                yes selected alpha1; docking scaffold
  7QNE human full-length synaptic alpha1beta3gamma2 GABA(A)R with Ro15-4513 and Mb38 alpha1beta3gamma2   EIE (Ro15-4513)                 2.70                  exact alpha1beta3gamma2 receptor; BZD-site ligand                yes      alternative alpha1 candidate
```

```
pdb_id                                                                 title                      composition co_crystal_ligand  resolution_angstrom                                                          bzd_site_assessment bzd_ligand_present                                             poc_decision
  9CTJ native human GABAA receptor beta2-alpha1-beta3-alpha2-gamma2 assembly beta2-alpha1-beta3-alpha2-gamma2              none                 3.74 mixed alpha1/alpha2 and beta2/beta3; contains alpha2-beta3-gamma2 local site                 no selected only as provisional alpha2 local-site construct
  9CX7 native human GABAA receptor beta3-alpha1-gamma2-beta3-alpha2 assembly beta3-alpha1-gamma2-beta3-alpha2              none                 3.30                 mixed alpha1/alpha2; not an exact alpha2beta3gamma2 receptor                 no                                             not selected
  9CXC native human GABAA receptor beta3-alpha1-gamma2-beta2-alpha2 assembly beta3-alpha1-gamma2-beta2-alpha2              none                 3.30                                          mixed alpha1/alpha2 and beta2/beta3                 no                                             not selected
  9CSB native human GABAA receptor beta3-alpha1-beta2-alpha2-gamma2 assembly beta3-alpha1-beta2-alpha2-gamma2              none                 3.34                                          mixed alpha1/alpha2 and beta2/beta3                 no                                             not selected
```

## Selection

6HUP was selected for α1β3γ2 because it is an exact receptor composition and provides a bound diazepam molecule at the α1/γ2 BZD site. The docking box is its DZP coordinate extent plus a fixed 6 Å margin (minimum 20 Å per axis), recorded in `docking_boxes.csv`.

An exact α2β3γ2 structure was not identified. 9CTJ was selected only as a provisional local-site source because it contains human α2, β3 and γ2 chains and preserves a native α2/γ2 interface, but it also contains α1 and β2. Chains C/D/E are extracted and the box is transferred by gamma2 alignment. This is a structural approximation and is explicitly not called a clean α2β3γ2 full receptor.

All structural limitations, box transfer details and receptor preparation limitations are repeated in `QC_REPORT.md` and `FINAL_REPORT.md`.