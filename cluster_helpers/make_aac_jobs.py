"""Make the Slurm script of the analysis of each molecule and energy, in LOCAL_BASE/<molecule>/.

One job does everything: it analyses all finished trajectories in parallel (`breakdown analysedemonfolder`, about a second per
file) and then makes the statistics, plots, drawings of ALL structures and the coordinates (`breakdown makechartsandstats`).
It uses the environment made by setup_cluster.py and installs nothing itself.
"""
import os

import cluster_config as cfg
from cluster_common import energyText, localDir, makeParser, nodeLine, runDirName, selectedJobs, timePerStepFs


def analysisScript(mol, abbr, energy):
    ev = energyText(energy)
    runDir = f'{cfg.REMOTE_BASE}/{mol}/{runDirName(energy)}'
    group = cfg.FUNCTIONAL_GROUPS.get(abbr)
    groupOption = f' -f {group}' if group else ''
    bond = cfg.BOND_OPTIONS
    modules = '\n'.join(f'module load {m}' for m in cfg.MODULES)
    return f"""#!/bin/bash
#SBATCH --job-name=bd_{abbr}_{ev}ev
{nodeLine()}#SBATCH --nodes 1
#SBATCH --cpus-per-task={cfg.ANALYSIS_CPUS}
#SBATCH --time={cfg.ANALYSIS_TIME}
#SBATCH --mem={cfg.ANALYSIS_MEMORY}
#SBATCH --output={cfg.REMOTE_BASE}/error_files/analysis_%j.out
#SBATCH --error={cfg.REMOTE_BASE}/error_files/analysis_%j.err

set -euo pipefail
module purge
{modules}
source "{cfg.REMOTE_BASE}/env/bin/activate"

RUN_DIR="{runDir}"
breakdown --version

# How many trajectories are there? This job may start after some of the MD runs failed.
N_MOL=$(ls "$RUN_DIR"/deMon_[0-9]*.mol 2>/dev/null | wc -l)
echo "$N_MOL of {cfg.N_RUNS} trajectories found in $RUN_DIR"
if [ "$N_MOL" -lt {cfg.N_RUNS} ]; then
    echo "WARNING: only $N_MOL of {cfg.N_RUNS} runs finished, the statistics use those" >&2
fi
# Runs that start in the same second can get the same random velocities: look for identical first frames
DUPLICATES=$(for f in "$RUN_DIR"/deMon_[0-9]*.mol; do head -n 20 "$f" | md5sum; done | sort | uniq -d | wc -l)
if [ "$DUPLICATES" -gt 0 ]; then
    echo "WARNING: $DUPLICATES groups of trajectories have an identical first frame, they may be copies of each other" >&2
fi

# Bond distances and options are recorded in every log and in full_log.txt
breakdown analysedemonfolder "$RUN_DIR" -p 'deMon_[0-9]*.mol' -j "$SLURM_CPUS_PER_TASK" \\
    -fd {bond['formation_distance']} -bd {bond['break_distance']} -ml {bond['min_lifetime']} \\
    || echo "WARNING: some trajectories could not be analysed, see above" >&2

# All structures are drawn on purpose, rare ones can be important transitions in the flow charts
breakdown makechartsandstats "$RUN_DIR" {mol}_{ev}ev {abbr}{groupOption} \\
    -t {timePerStepFs()} --draw --coordinates
"""


def main():
    args = makeParser('Make the analysis Slurm scripts.').parse_args()
    for mol, abbr, energy in selectedJobs(args):
        if not os.path.exists(localDir(mol)):
            print(f'{mol}: folder {localDir(mol)} not found, skipped')
            continue
        name = f'analysis_{abbr}_{energyText(energy)}ev.slurm'
        if args.dry_run:
            print(f'  [dry run] would write {os.path.join(localDir(mol), name)}')
            continue
        with open(os.path.join(localDir(mol), name), 'w') as file:
            file.write(analysisScript(mol, abbr, energy))
        print(f'{mol}, {energyText(energy)} eV: wrote {name}')


if __name__ == '__main__':
    main()
