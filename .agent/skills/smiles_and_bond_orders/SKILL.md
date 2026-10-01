---
name: smiles_and_bond_orders
description: >
  Why the SMILES strings in breakdown logs only describe connectivity, how atoms and bonds are decided from
  coordinates, and how to turn a log SMILES into a proper RDKit molecule without RDKit adding hydrogens to radicals.
  Read before touching bond_adder.py, create_smiles_from_dataframe.py, molecule_utils.py, drawing or xyz generation.
---

# SMILES strings and bond orders

## Facts about the simulations
- Molecules are (radical) cations, MD at high energy, so radicals and open shell fragments are normal.
- Every hydrogen is an explicit atom in the trajectory. RDKit must never add hydrogens.
- Only the geometry is known, bond orders are not. A stretched bond looks like a longer bond, not a lower order.

## Bond finding (`bond_adder.py`)
- Pair distances for all steps at once (chunked). A pair bonds below the formation distance (1.63 Å, C-S 1.9 Å) and stays
  bonded until above the break distance (2.5 Å): `bound[t] = form[t] or (bound[t-1] and keep[t])`, the only loop is over time.
- Hydrogen: keeps only its closest currently bonded non-hydrogen atom, if there is none the closest bonded hydrogen
  with a lower index (so H2 is one bond). Works for any atom order.
- Break distance is one global value (`--break_distance`), there are no element specific options on purpose.

## SMILES generation (`create_smiles_from_dataframe.makeSmiles`)
- Bond orders written: 2 only for a terminal O or S bonded to a heavy atom (`C=O`, `C=S`), 1 for everything else.
  Carbon-carbon bonds are all single, so a SMILES is really a labelled graph, and two states are "the same isomer" when
  their graphs are isomorphic (same canonical string).
- Hydrogen atoms get `SetNoImplicit(True)`, otherwise a lone H is written `[HH]`. Other atoms must *not* get it, that changes the
  canonical strings compared with all existing logs (which is how logs made by different versions stay mergeable).
- RDKit writes `[SH]`-like brackets for atoms with odd valence (S with three neighbours). That H is not real. `countAtoms`
  ignores it, `moleculeFromSmiles` removes it.
- One RDKit call per distinct bond pattern (bytes of the packed bond row are the cache key), states that recur are free.
- Event labels: more `.` than the step before = Fragmentation, fewer = Recombination, else Isomerization.
  `--min_lifetime N` replaces states that last fewer than N steps by the state before them (noise filter).

## Getting a real molecule from a SMILES: `molecule_utils.moleculeFromSmiles`
`Chem.MolFromSmiles(smiles)` sanitizes and gives every C/O/S with too few bonds extra implicit hydrogens, which is wrong here.
`moleculeFromSmiles` instead: sets no implicit Hs, makes all bonds of O/S with 2+ neighbours single (this also cleans legacy
strings with `=O=`), raises the order between neighbouring atoms that both have free valence using a maximum matching
(networkx), marks what remains as radical electrons, then sanitizes (aromaticity gets perceived). Use it for drawings and 3D
coordinates (`coordinates_file_maker`, which embeds in 3D and optimises with MMFF, or UFF when MMFF has no parameters).
Kekulé structure choice is arbitrary but deterministic, do not use it to decide whether two states are equal, compare SMILES strings.

## Known limits
- Bond orders and aromaticity are not tracked, so states differing only in bond order are one state.
- Valences known: H 1, C 4, N 3/5, O 2, S 2/4/6 (`VALENCES`). Add other elements there and in `NOMINAL_MASSES` if needed.
- The functional group detection in `isomer_analyser.py` is connectivity based and only knows O and S groups.
