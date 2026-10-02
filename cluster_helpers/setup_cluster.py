"""One-time setup on the cluster: a Python environment with BREAKDOWN installed in it.

The analysis jobs only activate this environment. They never install anything, so many jobs at the same time cannot
damage it. Run again (it is safe) to update BREAKDOWN, and pin BREAKDOWN_INSTALL to a tag or commit in
cluster_config.py so that results can be reproduced.
"""
import os

import cluster_config as cfg
from cluster_common import makeParser, q, run, selectedJobs, ssh, upload


def main():
    parser = makeParser('Create the environment on the cluster.', energies=False, molecules=False)
    parser.add_argument('--local-source', metavar='FOLDER',
                        help='install from this local BREAKDOWN folder (built here and uploaded) instead of from '
                             'BREAKDOWN_INSTALL, for clusters that cannot reach GitHub')
    args = parser.parse_args()
    selectedJobs(args)  # sets the dry run mode
    modules = '; '.join(f'module load {m}' for m in cfg.MODULES)
    env = f'{cfg.REMOTE_BASE}/env'
    if args.local_source:
        wheelDir = os.path.join(args.local_source, 'dist_for_cluster')
        run(['python', '-m', 'pip', 'wheel', '--no-deps', '-w', wheelDir, args.local_source])
        ssh(f'mkdir -p {q(cfg.REMOTE_BASE + "/wheels")}')
        wheels = [os.path.join(wheelDir, w) for w in os.listdir(wheelDir)] if os.path.isdir(wheelDir) else []
        upload(wheels, f'{cfg.REMOTE_BASE}/wheels')
        install = f'pip install --upgrade {q(cfg.REMOTE_BASE + "/wheels")}/breakdown-*.whl'
    else:
        install = f'pip install --upgrade {q(cfg.BREAKDOWN_INSTALL)}'
    commands = f"""set -e
module purge
{modules}
mkdir -p {q(cfg.REMOTE_BASE)}
cd {q(cfg.REMOTE_BASE)}
[ -d env ] || python -m venv env
source env/bin/activate
pip install --upgrade pip
{install}
echo "Installed:"
breakdown --version"""
    print(ssh(commands))


if __name__ == '__main__':
    main()
