"""Send the MD inputs to the cluster and submit the MD array jobs. The job ids are saved in jobs_<energy>ev.json so that
queue_aas_jobs.py can make the analysis wait for them."""
import os

import cluster_config as cfg
from cluster_common import (energyText, inpName, localDir, makeParser, parseJobId, q, remoteDir, runDirName, saveJobRecord,
                            selectedJobs, ssh, upload)


def main():
    args = makeParser('Upload MD inputs and submit the MD jobs.').parse_args()
    jobs = selectedJobs(args)
    print(f'Deploying to {cfg.REMOTE_HOST}...')
    ssh(f'mkdir -p {q(cfg.REMOTE_BASE + "/error_files")}')
    if cfg.UPLOAD_EXECUTABLE:
        print(f'Uploading the executable to {cfg.REMOTE_BASE}...')
        upload([cfg.LOCAL_EXECUTABLE], cfg.REMOTE_BASE)

    for mol, abbr, energy in jobs:
        inpFile = inpName(mol, energy)
        slurmFile = f'{abbr}_{energyText(energy)}ev.slurm'
        settingsFile = f'md_settings_{energyText(energy)}ev.json'
        localFiles = [os.path.join(localDir(mol), name) for name in (inpFile, slurmFile, settingsFile)]
        missing = [path for path in localFiles if not os.path.exists(path)]
        if missing and not args.dry_run:
            print(f'Skipping {mol}, {energyText(energy)} eV: {missing} not found, run make_md_inputs.py first')
            continue

        print(f'\n--- {mol} ({abbr}), {energyText(energy)} eV ---')
        ssh(f'mkdir -p {q(remoteDir(mol))}')
        upload(localFiles, remoteDir(mol))
        output = ssh(f'cd {q(remoteDir(mol))} && sbatch --parsable {q(slurmFile)}')
        jobId = parseJobId(output)
        print(f'  MD array job submitted: {jobId}')
        saveJobRecord(mol, energy, md_job_id=jobId, n_runs=cfg.N_RUNS, run_dir=f'{remoteDir(mol)}/{runDirName(energy)}')

    print('\nDone. queue_aas_jobs.py will make the analysis wait for these jobs.')


if __name__ == '__main__':
    main()
