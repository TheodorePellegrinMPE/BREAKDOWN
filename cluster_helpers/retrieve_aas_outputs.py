"""Download the results of the analysis jobs. By default only the small files are fetched: run_summary.json (all that the
compare and flow chart commands need) and the plots. Choose more with --what."""
import os

import cluster_config as cfg
from cluster_common import (download, energyText, loadJobRecord, localDir, makeParser, remoteDir, runDirName, selectedJobs,
                            ssh)

# What can be fetched: name -> (path in the run folder, is it a folder)
ITEMS = {
    'summary': ('run_summary.json', False),
    'plots': ('plots', True),
    'drawings': ('drawings', True),
    'coordinates': ('coordinates', True),
    'full_log': ('full_log.txt', False),     # 100 MB or more, the compare commands do not need it
    'logs': (None, False),                   # every per-run *_log.txt
}


def jobState(jobId):
    """State of a job according to Slurm (COMPLETED, RUNNING, FAILED, ...), or unknown."""
    output = ssh(f'sacct -j {jobId} --format=State -n -X -P', check=False)
    states = sorted({line.strip().split()[0] for line in output.splitlines() if line.strip()})
    return ','.join(states) or 'unknown'


def main():
    parser = makeParser('Download the analysis results.')
    parser.add_argument('--what', nargs='+', choices=list(ITEMS), default=['summary', 'plots'],
                        help='what to download, default: summary plots')
    parser.add_argument('--force', action='store_true', help='download even if the analysis job has not completed')
    args = parser.parse_args()
    print(f'Retrieving analysis results from {cfg.REMOTE_HOST}...')

    for mol, abbr, energy in selectedJobs(args):
        record = loadJobRecord(mol, energy)
        jobId = record.get('analysis_job_id')
        if jobId and not args.force and not args.dry_run:
            state = jobState(jobId)
            print(f'\n--- {mol}, {energyText(energy)} eV: analysis job {jobId} is {state} ---')
            if state != 'COMPLETED':
                print('  not completed, skipped (use --force to download anyway)')
                continue
        else:
            print(f'\n--- {mol}, {energyText(energy)} eV ---')

        remoteRun = f'{remoteDir(mol)}/{runDirName(energy)}'
        localRun = os.path.join(localDir(mol), runDirName(energy))
        for item in args.what:
            path, isFolder = ITEMS[item]
            print(f'  Downloading {item}...')
            if item == 'logs':
                download(f'{remoteRun}/*_log.txt', localRun)
            else:
                download(f'{remoteRun}/{path}', localRun, recursive=isFolder)

    print('\nDone. The compare commands can now be run on the downloaded run_summary.json files.')


if __name__ == '__main__':
    main()
