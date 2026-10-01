---
name: full_log_interpreter
description: >
  Teaches an AI how to extract structured information from breakdown's
  full_log.txt output files produced by MD simulation analysis. Covers the
  complete file format, every line type, every statistics dictionary, and
  provides guidance on writing Python scripts to query the file safely without
  loading it into memory.
---

# full_log.txt Interpreter Skill

## ⚠️ Critical Rule: Never Read full_log.txt Directly

`full_log.txt` is **extremely large** — often millions of lines — because it
concatenates the per-step SMILES records for every simulation run and then
appends many JSON statistics dictionaries.  
**Never use `open()` or any tool to read the whole file into memory.**  
Instead, **always write a new Python script** that processes the file line-by-
line (using a streaming approach) or uses `str.split` / regex on targeted
chunks, and run that script to get the answer.

---

## File Generation Pipeline

The pipeline is:

```
.mol files (MD simulation output)
  -> AutoDetectIsomerization.autoAnalyseFile()    [per-run log.txt]
  -> GenerateChartsStats.autoAnalyseInFolder()    [full_log.txt]
```

`full_log.txt` is the concatenation of all individual `*_log.txt` files
followed by a set of appended statistics blocks.

---

## full_log.txt Format — Complete Specification

### 1. Run Block (repeated once per simulation run)

Each run starts immediately after the previous run ends. There is no explicit
separator *between* runs other than the `In file:` header that starts each one.

```
In file: <filename>.mol                      ← run header
At step 0:               <smiles>            ← first step, no event
At step 1:               <smiles>
At step 2: Isomerization <smiles>            ← molecule changed (same # of fragments)
<isomer analysis line>                       ← immediately after Isomerization or Fragmentation
At step 3: Isomerization <smiles>
<isomer analysis line>
At step 4: Fragmentation <smiles>            ← new fragment appeared (more '.' in smiles)
At step 9: Recombination <smiles>            ← fragments joined (fewer '.' in smiles)
<isomer analysis line>
...
At step N:               <smiles>            ← last recorded step
Final molecules present: <formula>, <formula>, ...   ← chemical-formula list of end fragments
```

Blank lines may appear between run blocks.

#### `In file:` line
- Format: `In file: <filename>.mol `  (trailing space present)
- Marks the beginning of a new simulation run.

#### `At step` lines
- Format: `At step <N>: <event> <smiles>`
- `<N>` — integer step index, 0-based
- `<event>` — one of:
  - 13 spaces (`'             '`) — no change since last step
  - `'Isomerization'` — SMILES changed, same number of `.`-separated fragments
  - `'Fragmentation'` — SMILES changed, *more* `.`-separated fragments than before
  - `'Recombination'` — SMILES changed, *fewer* `.`-separated fragments than before
    (added in newer versions, older logs never have it; treat it like the other
    two when counting events)
- `<smiles>` — the **full SMILES string** for this step. It encodes which atoms are
  bonded, bond orders are only a guess, and hydrogens are separate `[H]` atoms. An
  `H` inside other brackets, like `[SH]`, is added by RDKit for unusual valences and is
  not a real hydrogen. To draw it or get 3D coordinates use
  `breakdown.molecule_utils.moleculeFromSmiles`, never `Chem.MolFromSmiles`, which adds
  hydrogens to radicals. If the molecule has
  fragmented, the SMILES contains dot-separated sub-SMILES, e.g.
  `CC(=O)O.CH4`. The SMILES is always the **last whitespace-separated token**
  on the line.

> **Parsing tip:** `smiles = line.split()[-1]`

#### Isomer analysis line (optional, immediately after Isomerization / Fragmentation / Recombination)
- Appears on the line *immediately following* an `At step ... Isomerization`,
  `At step ... Fragmentation` or `At step ... Recombination` line.
- Single-line concatenation of detected structural features, e.g.:
  ```
  Number of 3-carbon rings: 1. Hydroxyl group count: 2. Ketone count: 1. Ether count: 1.
  ```
- Possible fields (all optional, only shown when count > 0):
  - `Number of <N>-carbon rings: <count>.` — ring sizes observed (3, 4, 5, 6, …)
  - `Hydroxyl group count: <count>.`
  - `Methylene group count: <count>.`
  - `Tertiary carbon count: <count>.`
  - `Ketone count: <count>.`
  - `Ketene count: <count>.`
  - `Ether count: <count>.`
  - `Aldehyde count: <count>.`
  - `Carboxyl count: <count>.`
  - `Epoxide count: <count>.`
  - Sulfur analogues: `Thiol group count`, `Thioketone count`, `Thioaldehyde count`,
    `Thioether count`, `Thiirane count`
- If no detectable features, an empty string is stored (the line may be blank
  or absent).

#### `Final molecules present:` line
- Format: `Final molecules present: <formula1>, <formula2>, ...`
- Chemical formulas (e.g. `C8H8O`, `C2H2`) for all fragments at the *last
  recorded step*, sorted largest-first by atom count.
- One per run, at the very end of the run block.

---

### 2. Statistics Blocks (appended after all runs)

After all run blocks, `full_log.txt` has one header line followed by a JSON
block for each statistics dictionary. The pattern is always:

```
<Human-readable header sentence>\n
<JSON object, indented 4 spaces>\n
```

The statistics blocks appear in this **fixed order** (using the exact header
strings as delimiters):

| Order | Header line | JSON content |
|-------|-------------|--------------|
| 1 | `Automatically generated molecule names.` | `{ "<smiles>": "<name>", ... }` |
| 2 | `Total number of steps at which an isomer is present.` | `{ "<name>": <int>, ... }` |
| 3 | `Total number of runs in which an isomer appeared.` | `{ "<name>": <int>, ... }` |
| 4 | `Number of runs at which each chemical fragment was present at the end.` | `{ "<formula>": <int>, ... }` |
| 5 | `Total number of transitions to an isomer from an isomer.` | `{ "<name>": { "<name>": <int>, ... }, ... }` |
| 6 | `Net number of transitions to an isomer from an isomer.` | `{ "<name>": { "<name>": <int>, ... }, ... }` |
| 7 | `Number of transitions only including fragmentation.` | `{ "<name>": { "<name>": <int>, ... }, ... }` |
| 8 | `Chemical names of fragments of this isomerization state.` | `{ "<name>": ["<formula>", ...], ... }` |
| 9 | `Molecular masses of different observed isomers.` | `{ <mass_int>: ["<formula>", ...], ... }` |
| 10 (optional) | `Runs containing original functional group <group> at step.` | `{ "<step_int>": <int>, ... }` |
| 11 (optional) | `Runs with each <heteroatom> in each final group and/or fragment.` | `{ "<desc>": <int>, ... }` |
| 12 | `First step of each event type in each run.` | `{ "<file>.mol": {"isomerization": <step or null>, "fragmentation": ..., "recombination": ..., "anyChange": ..., "lastStep": <int>} }` (newer versions only) |

> Blocks 10 and 11 only exist when `functionalGroup` was set during analysis.

#### Isomer Naming Convention
Isomers are automatically named with the pattern:
`<moleculeName>_<letter><number>`

- `<moleculeName>` — the molecule name passed to the pipeline
- `<letter>` — `A` = no carbon loss from the major fragment, `B` = 1 carbon
  lost, `C` = 2 carbons lost, etc. (alphabetically by carbon loss in minor
  fragments)
- `<number>` — integer, 1-indexed; identifies which isomer within the same
  carbon-loss group, ordered by prevalence (most common = 1)

Example: `propanol_A1` is the most common intact propanol isomer;
`propanol_B1` is the most common isomer that has lost one carbon.

---

## Key Parsing Patterns

### Streaming over `At step` lines
```python
with open('full_log.txt', 'r') as f:
    for line in f:
        if line.startswith('At step'):
            smiles = line.split()[-1].rstrip()
            step = int(line.split()[2].rstrip(':'))
            event = 'Fragmentation' if 'Fragmentation' in line else \
                    'Isomerization' if 'Isomerization' in line else None
            # process ...
```

### Extracting a statistics block by header
```python
import json, re

def extract_stats_block(path, header_line):
    """Stream file, return the JSON dict following the given header line."""
    with open(path, 'r') as f:
        capture = False
        json_lines = []
        for line in f:
            if capture:
                json_lines.append(line)
                # JSON objects end when we return to zero brace depth
                text = ''.join(json_lines)
                try:
                    return json.loads(text)
                except json.JSONDecodeError:
                    continue  # keep accumulating
            if line.strip() == header_line:
                capture = True
    return None
```

### Iterating over run blocks
```python
def iter_runs(path):
    """Yield (filename, lines_in_run) for each simulation run."""
    with open(path, 'r') as f:
        current_file = None
        current_lines = []
        for line in f:
            if line.startswith('In file:'):
                if current_file is not None:
                    yield current_file, current_lines
                current_file = line.strip().split('In file: ')[1]
                current_lines = []
            elif current_file is not None:
                current_lines.append(line)
        if current_file is not None:
            yield current_file, current_lines
```

### Checking if a line is an isomer analysis line
```python
def is_analysis_line(line):
    """True if the line is the functional-group analysis following an event."""
    markers = [
        'count:', 'rings:', 'Hydroxyl', 'Ketone', 'Ketene',
        'Ether', 'Aldehyde', 'Carboxyl', 'Epoxide', 'Methylene', 'Tertiary'
    ]
    return any(m in line for m in markers)
```

---

## Common Query Recipes

When the user asks for something, **write a Python script** that does one of
the following patterns and run it. Never load the whole file into a variable.

### "How many runs are there?"
Stream and count `In file:` lines.

### "What SMILES appeared most often overall?"
Stream `At step` lines, accumulate a `collections.Counter` on the last token.

### "What was the molecule at step N in run X?"
Stream with `iter_runs`, find the matching file, look for `At step N:`.

### "What are the top-K most common isomers by step count?"
Extract the `Total number of steps at which an isomer is present.` JSON block.

### "What fraction of runs ended with fragmentation?"
Count `Final molecules present:` lines where the value contains a comma (i.e.
more than one fragment).

### "How many isomerization events happened in total?"
Stream and count lines where `'Isomerization' in line` (excluding analysis
lines by checking `line.startswith('At step')`).

### "What functional groups were present at transitions?"
Stream pairs: when you see `At step ... Isomerization`, read the *next* line
as the analysis line and parse `<group> count: <N>` with regex.

### "Get the full transition matrix (isomerizationStats)"
Extract the `Total number of transitions to an isomer from an isomer.` JSON
block — already fully computed and stored in the file.

---

## Complete Statistics Dictionary Reference

### `isomerStats` (block 2)
```json
{
  "propanol_A1": 18432,
  "propanol_A2": 4210,
  ...
}
```
Keys: isomer names. Values: total `At step` lines in the file where that SMILES
appeared (summed across all runs and all steps).

### `isomerRuns` (block 3)
```json
{
  "propanol_A1": 198,
  ...
}
```
Keys: isomer names. Values: number of *distinct runs* in which the isomer
appeared at least once.

### `endFragmentStats` (block 4)
```json
{
  "C3H8O": 120,
  "C2H6": 30,
  ...
}
```
Keys: chemical formula strings. Values: number of runs in which that fragment
was present in the `Final molecules present:` line.

### `isomerizationStats` (block 5)
```json
{
  "propanol_A2": {
    "propanol_A1": 45
  },
  ...
}
```
Nested dict: `isomerizationStats[destination][source]` = total number of
transitions *from* source *to* destination (including both Isomerization and
Fragmentation events).

### `netIsomerizationStats` (block 6)
Same structure as block 5 but only entries where
`transitions(A→B) > transitions(B→A)`, storing the *net* surplus.

### `fragmentationStats` (block 7)
Same structure, but only entries where the transition crossed a fragmentation
boundary (different number of `.` in source vs destination SMILES).

### `isomerChemicalNames` (block 8)
```json
{
  "propanol_A1": ["C3H8O"],
  "propanol_B1": ["C2H6O", "CH2"],
  ...
}
```
Maps each isomer name to the list of constituent fragment chemical formulas.

### `isomerMasses` (block 9)
```json
{
  "60": ["C3H8O"],
  "46": ["C2H6O"],
  ...
}
```
Maps integer molecular masses to lists of chemical formulas observed at that
mass. Keys are integer strings; values are lists. Sorted descending by mass.

### `groupSteps` (block 10 — optional)
```json
{
  "0": 200,
  "1": 198,
  "50": 145,
  ...
}
```
Keys: step numbers (as strings). Values: number of runs in which the tracked
functional group was *still present* at that step.

### `endingGroups` (block 11 — optional)
```json
{
  "hydroxyl": 80,
  "ketone / CO": 35,
  "propanol_A1": 14,
  ...
}
```
Keys describe what the oxygen ended up in at the end of each run (either a
named functional group, a fragment formula, or the original group name if no
change). Values: run counts. Groups are taken from the analysis of the state after the
run's *last* event only. A run with events but no tracked group or heteroatom fragment at
the end is counted as `"other"`. Older versions used the last analysis line that had any
tracked group, which could report a group the molecule had already lost.

---

## Script Template

Use this as a starting point for any new extraction script:

```python
#!/usr/bin/env python3
"""
Extract <description> from full_log.txt.
Usage: python extract_<name>.py /path/to/full_log.txt
"""
import sys
import re
import json
from collections import defaultdict

LOG_PATH = sys.argv[1] if len(sys.argv) > 1 else 'full_log.txt'

# ---- streaming pass ----
results = defaultdict(int)

with open(LOG_PATH, 'r') as f:
    prev_line = ''
    for line in f:
        line = line.rstrip('\n')
        
        if line.startswith('At step'):
            smiles = line.split()[-1]
            step   = int(line.split()[2].rstrip(':'))
            is_iso = 'Isomerization' in line
            is_fra = 'Fragmentation' in line
            # <your logic here>

        elif line.startswith('In file:'):
            run_name = line.strip().split('In file: ')[1]
            # <reset per-run accumulators here>

        elif line.startswith('Final molecules present:'):
            fragments = line.split(': ')[1].split(', ')
            # <end-of-run logic here>

        prev_line = line

# ---- print / save results ----
print(json.dumps(dict(results), indent=2))
```

---

## Notes and Gotchas

1. **Trailing whitespace**: The `In file:` line has a trailing space before
   `\n`. Strip with `.strip()` before parsing.
2. **Blank lines**: Appear between run blocks. Guard with `if not line.strip(): continue`.
3. **Step numbering is per-run**: Step 0 always means the very first frame of
   *that* run; it restarts for every `In file:` block.
4. **Statistics keys are already isomer names, not raw SMILES**: In all
   statistics blocks, raw SMILES strings have been replaced with friendly names
   like `propanol_A1`. Only the isomer-names dict (block 1) maps raw SMILES →
   name.
5. **Molecular mass calculation** uses integer masses: H=1, C=12, N=14, O=16,
   S=32.
6. **Early termination**: A run is cut short (fewer `At step` lines) if the
   molecule accumulates more than 3 carbon-containing fragments for more than
   20 consecutive steps. The `Final molecules present:` line still appears at
   the end of such truncated runs.
7. **The analysis line** (after Isomerization/Fragmentation) is a *single line*
   with no leading `At step` prefix. Never confuse it with the next step line.
8. **Dot-separated SMILES**: Fragments within one step are joined by `.`.
   Counting dots tells you how many fragment boundaries there are
   (`num_fragments = smiles.count('.') + 1`).

---

## run_summary.json (shortcut)

Newer versions write `run_summary.json` next to `full_log.txt`. It holds every statistics block above
under `sections` (key = heading text, value = the parsed dictionary), plus `firstEvents`, `timePerStepFs`, `numRuns`
and `elementalIntact`. **Prefer it over `full_log.txt`** when you only need statistics: `json.load` it, or use
`breakdown.run_summary.loadSection(folder, heading)`, which uses the summary when it is not older than `full_log.txt`
and otherwise streams `full_log.txt` to the heading. Per-step SMILES are only in the run logs.
