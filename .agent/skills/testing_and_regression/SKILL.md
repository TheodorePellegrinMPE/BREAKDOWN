---
name: testing_and_regression
description: >
  How to run and test breakdown: the conda environment, the example trajectories and logs
  that act as regression references, the compare tool, and the pitfalls (huge files, overwriting
  example logs, which checkout is imported). Use before changing bond finding, SMILES generation,
  statistics or plots.
---

# Testing breakdown

## Environment

```bash
source ~/anaconda3/etc/profile.d/conda.sh && conda activate aac
export PYTHONPATH=/home/tpelleg/Documents/autoastrochemslopified/autoastrochem/src MPLBACKEND=Agg
```
`aac` has the old `autoastrochem` package installed editable from another checkout (no `breakdown` in it), so without `PYTHONPATH` you test the wrong code or get an import error.
Run the CLI as `python -m breakdown.main <command>`. Command names are lowercase.

## Unit tests

`pytest` from the project folder (the one with `pyproject.toml`) (config in `pyproject.toml`). They are fast (about 1 s) and generate their own data:
`tests/test_bonds_and_smiles.py` (bond hysteresis, hydrogens, labels, `moleculeFromSmiles`),
`tests/test_pipeline.py` (mol -> log -> full_log on a toy molecule),
`tests/test_summary_and_survival.py` (Kaplan-Meier, `run_summary.json` vs streaming `full_log.txt`, stale summary),
`tests/test_diverse_molecules.py` (RDKit-built molecules: alcohol, ketone, acid, ether, epoxide, thio analogues, disulfide, pyridine, amide, F/Cl/Br/I compounds;
scripted H shifts and fragmentation; shuffled atom order).
Add a case there whenever a bug is fixed. `writeMol` in `tests/conftest.py` writes the deMon `.mol` layout.

## Real reference data (`../example_out_and_log_files/`)

Per molecule: `<name>_50ps.inp` (input), `<name>_traj7.mol` (trajectory, 10001 frames, C12H8 has 50001),
`<name>_traj7.out` (energies, not coordinates), `<name>_traj7_log.txt` (reference output of the original code),
and `anthracene_15ev_full_log.txt`, `phenol_5ev_full_log.txt` (100 runs each).

**Never write into that folder** (logs would be overwritten). Copy the `.mol` files into a scratch directory
(the Claude scratchpad) first. The `.mol` files are 16–118 MB, and `.out` and `full_log.txt` files are huge: do not `cat`/Read them,
use `head`, `grep -c`, or scripts that stream.

### Regression recipe

```bash
mkdir /tmp/x && cp ../example_out_and_log_files/*.mol /tmp/x && cd /tmp/x
python -m breakdown.main analysedemonfolder . -j 4          # 9 files in about 3 s
python tools/compare_logs.py /tmp/x/phenol_5ev_traj7_log.txt ../example_out_and_log_files/phenol_5ev_traj7_log.txt
```
`tools/compare_logs.py` compares event labels and the bonded graph, not the exact string. With default settings
all 9 example runs give 0 mismatches. The only expected differences are steps where the old code said `Isomerization`
and the new code says `Recombination` (fragments joined), and final formulas that now include all elements
(the old logs were made with the default `-e C,H,O`, so S/N were dropped).
One more known difference: the fluorene reference has an H with two bonds at step 2015 (`[H][H]C1...`), the new code never lets a hydrogen have two bonds, so it reports that fragmentation one step earlier (3 event mismatches for fluorene instead of 1).
If more differ, bond finding or SMILES generation changed behaviour, find out why before going on.
Exact text equality only holds for the SMILES strings when no S-H/O-H states are involved, see `smiles_and_bond_orders`.

### Statistics regression

Split a full log into per-run logs (`re.split('(?=^In file: )', text_before_'Automatically generated molecule names.')`),
write each as `<name>_log.txt` in a scratch folder, run `makechartsandstats . NAME MOL -s <carbons>` and compare the
part of `full_log.txt` after `Automatically generated molecule names.` between old and new code. `git stash`/`git worktree`
the old code and run it with its own `PYTHONPATH`. File order matters for tie-breaking of isomer names, the code sorts the file names.
`-f <group>` also exercises the functional group statistics.

## Performance expectations

One 10k-step run: about 0.3–0.6 s. 100 runs of `makechartsandstats`: a few seconds. If something takes tens of seconds
a per-step Python loop or a repeated whole-file scan has crept in.

## Checking plots

PDFs: `pdftoppm -png -r 50 file.pdf out` then look at the PNG. Flow charts need graphviz `dot` on the PATH.
