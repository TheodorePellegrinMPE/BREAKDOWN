"""Helpers shared by the cluster scripts: settings, running commands, and the record of submitted jobs."""
import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from datetime import datetime, timezone

try:
    import cluster_config as cfg
except ImportError:
    sys.exit('cluster_config.py is missing. Copy cluster_config.example.py to cluster_config.py and edit it.')

STATE = {'dry_run': False}


def makeParser(description, energies=True, molecules=True):
    """Argument parser with the options every script has: --dry-run, and which energies and molecules to use."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--dry-run', action='store_true', help='print what would be done, change nothing')
    if energies:
        parser.add_argument('--ev', type=float, nargs='+', default=None,
                            help=f'kinetic energies in eV to use, default is ENERGIES_EV = {cfg.ENERGIES_EV}')
    if molecules:
        parser.add_argument('--molecule', nargs='+', default=None,
                            help=f'molecule folders to use, default is all of MOLECULES = {list(cfg.MOLECULES)}')
    return parser


def selectedJobs(args):
    """The (molecule, short name, energy) combinations chosen by the options and the settings."""
    STATE['dry_run'] = args.dry_run
    energies = args.ev if getattr(args, 'ev', None) else cfg.ENERGIES_EV
    molecules = getattr(args, 'molecule', None) or list(cfg.MOLECULES)
    unknown = [m for m in molecules if m not in cfg.MOLECULES]
    if unknown:
        sys.exit(f'Unknown molecules {unknown}, add them to MOLECULES in cluster_config.py')
    return [(mol, cfg.MOLECULES[mol], energy) for mol in molecules for energy in energies]


def energyText(energy):
    """25.0 -> '25', 12.5 -> '12.5', for file names."""
    return f'{energy:g}'


def timePerStepFs():
    """Simulated time between the frames of deMon.mol, the time per step of the analysis."""
    return cfg.MD_TIMESTEP_FS * cfg.MD_OUTPUT_INTERVAL


def totalPs():
    return cfg.MD_STEPS * cfg.MD_TIMESTEP_FS / 1000


def inpName(mol, energy):
    """File name of the MD input, with the simulated time in it."""
    return f'{mol}_{energyText(energy)}ev_{totalPs():g}ps.inp'


def runDirName(energy):
    return f'run0_{energyText(energy)}ev'


def localDir(mol):
    return os.path.join(cfg.LOCAL_BASE, mol)


def remoteDir(mol):
    return f'{cfg.REMOTE_BASE}/{mol}'


def run(command, check=True, capture=False):
    """Run a local command (a list). Stops with a clear message if it fails, unless check is False."""
    printable = ' '.join(shlex.quote(part) for part in command)
    if STATE['dry_run']:
        print(f'  [dry run] {printable}')
        return subprocess.CompletedProcess(command, 0, stdout='DRYRUN', stderr='')
    result = subprocess.run(command, capture_output=capture, text=True)
    if check and result.returncode != 0:
        sys.exit(f'Command failed ({result.returncode}): {printable}\n{result.stderr or ""}')
    return result


def ssh(command, check=True):
    """Run a shell command on the cluster, and return its output."""
    return run(['ssh', cfg.REMOTE_HOST, command], check=check, capture=True).stdout.strip()


def q(path):
    """Quote a path or name for the remote shell, names can contain spaces, commas and brackets."""
    return shlex.quote(str(path))


def upload(localPaths, remoteFolder):
    """Copy files to a folder on the cluster."""
    run(['rsync', '-az', '--partial', *map(str, localPaths), f'{cfg.REMOTE_HOST}:{remoteFolder}/'])


def download(remotePath, localFolder, recursive=False):
    """Copy a file or folder from the cluster into a local folder."""
    if not STATE['dry_run']:
        os.makedirs(localFolder, exist_ok=True)
    source = f'{cfg.REMOTE_HOST}:{remotePath}' + ('/' if recursive else '')
    target = os.path.join(localFolder, os.path.basename(remotePath)) if recursive else localFolder
    if recursive and not STATE['dry_run']:
        os.makedirs(target, exist_ok=True)
    run(['rsync', '-az', '--partial', source, target + ('/' if recursive else '')])


def nodeLine():
    """The #SBATCH line that pins a job to a node, or an empty string."""
    return f'#SBATCH -w {cfg.SLURM_NODE}\n' if cfg.SLURM_NODE else ''


def parseJobId(text):
    """Job id from `sbatch --parsable` output (it can be `123` or `123;cluster`)."""
    if text == 'DRYRUN':
        return 'DRYRUN'
    match = re.match(r'(\d+)', text.strip())
    if not match:
        sys.exit(f'Could not read a job id from sbatch output: {text!r}')
    return match.group(1)


def jobRecordPath(mol, energy):
    return os.path.join(localDir(mol), f'jobs_{energyText(energy)}ev.json')


def loadJobRecord(mol, energy):
    """The job ids saved when the jobs of a molecule and energy were submitted."""
    path = jobRecordPath(mol, energy)
    if os.path.exists(path):
        with open(path) as file:
            return json.load(file)
    return {}


def saveJobRecord(mol, energy, **entries):
    """Add entries (like md_job_id) to the job record of a molecule and energy."""
    if STATE['dry_run']:
        print(f'  [dry run] would save {entries} to {jobRecordPath(mol, energy)}')
        return
    record = loadJobRecord(mol, energy)
    record.update(entries)
    record['updated'] = datetime.now(timezone.utc).isoformat(timespec='seconds')
    with open(jobRecordPath(mol, energy), 'w') as file:
        json.dump(record, file, indent=2)
