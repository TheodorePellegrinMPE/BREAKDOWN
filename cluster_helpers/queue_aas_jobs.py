"""Send the analysis Slurm scripts to the cluster and submit them. The analysis waits for the MD job of the same molecule and
energy (its id is read from jobs_<energy>ev.json, saved by queue_md_jobs.py)."""
import os

import cluster_config as cfg
from cluster_common import (energyText, loadJobRecord, localDir, makeParser, parseJobId, q, remoteDir, saveJobRecord,
                            selectedJobs, ssh, upload)


def main():
    parser = makeParser('Upload and submit the analysis jobs.')
    parser.add_argument('--md-job-id', help='wait for this job instead of the one in the job record')
    parser.add_argument('--no-dependency', action='store_true', help='start the analysis without waiting for MD jobs')
    args = parser.parse_args()
    print(f'Deploying analysis jobs to {cfg.REMOTE_HOST}...')

    for mol, abbr, energy in selectedJobs(args):
        scriptName = f'analysis_{abbr}_{energyText(energy)}ev.slurm'
        localScript = os.path.join(localDir(mol), scriptName)
        if not os.path.exists(localScript) and not args.dry_run:
            print(f'Skipping {mol}, {energyText(energy)} eV: {scriptName} not found, run make_aac_jobs.py first')
            continue

        mdJobId = None if args.no_dependency else (args.md_job_id or loadJobRecord(mol, energy).get('md_job_id'))
        if not args.no_dependency and not mdJobId:
            print(f'Skipping {mol}, {energyText(energy)} eV: no MD job id known. Run queue_md_jobs.py first, give '
                  '--md-job-id, or use --no-dependency if the runs are already finished.')
            continue

        print(f'\n--- Analysis of {mol}, {energyText(energy)} eV ---')
        upload([localScript], remoteDir(mol))
        # afterany: the analysis also runs if some of the MD runs failed, it reports how many trajectories it found
        dependency = f'--dependency=afterany:{mdJobId} ' if mdJobId else ''
        output = ssh(f'cd {q(remoteDir(mol))} && sbatch --parsable {dependency}{q(scriptName)}')
        jobId = parseJobId(output)
        print(f'  Analysis job submitted: {jobId}' + (f' (after MD job {mdJobId})' if mdJobId else ''))
        saveJobRecord(mol, energy, analysis_job_id=jobId)

    print('\nDone. retrieve_aas_outputs.py downloads the results.')


if __name__ == '__main__':
    main()
