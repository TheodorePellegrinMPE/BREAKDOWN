"""Settings shared by all cluster helper scripts.

Copy this file to cluster_config.py and edit it. cluster_config.py is in .gitignore, so your own paths and host name stay
on your machine. Every script reads its settings from there, nothing is hard-coded in the scripts.
"""

# --- Where things are -------------------------------------------------------------------------------------------
LOCAL_BASE = '/path/to/local/md_stuff'              # one folder per molecule, with <molecule>.out of its optimization
REMOTE_HOST = 'my_cluster'                          # ssh host (alias from ~/.ssh/config)
REMOTE_BASE = '/path/on/cluster/molecular_dynamics'
LOCAL_EXECUTABLE = LOCAL_BASE + '/deMon.static-website.x'
LOCAL_BASIS = LOCAL_BASE + '/basis'                 # used by submit_opt_jobs.py, on the cluster it is REMOTE_BASE/basis
UPLOAD_EXECUTABLE = False                           # True to copy LOCAL_EXECUTABLE to REMOTE_BASE when queueing MD jobs

# --- What to run --------------------------------------------------------------------------------------------------
MOLECULES = {'fluorene': 'flr'}                     # folder name -> short name used in job names and plot names
FUNCTIONAL_GROUPS = {'flr': None}                   # short name -> group for makechartsandstats -f, or None
ENERGIES_EV = [25]                                  # kinetic energies (eV) given to the molecule, one set of runs each
N_RUNS = 100                                        # trajectories per molecule and energy
MAX_PARALLEL_RUNS = 20                              # at most this many MD array tasks at the same time (%N in sbatch)

# --- deMon-Nano input of the MD runs --------------------------------------------------------------------------------
MD_TIMESTEP_FS = 0.5
MD_STEPS = 100000
MD_OUTPUT_INTERVAL = 10                             # steps between frames in deMon.mol, deMon-Nano's default
CHARGE = 1
MULTIPLICITY = 2
MD_FERMI_K = 1500
OPT_FERMI_K = 9670                                  # used by submit_opt_jobs.py

# --- Slurm ---------------------------------------------------------------------------------------------------------
SLURM_NODE = None                                   # e.g. 'node03' to pin jobs to one node, None to let Slurm choose
MODULES = ['python-venv']                           # `module load` these in the analysis jobs
MD_TIME = '10:00:00'                                # per run
MD_MEMORY = '1gb'
ANALYSIS_TIME = '01:00:00'
ANALYSIS_MEMORY = '16gb'
ANALYSIS_CPUS = 8

# --- BREAKDOWN on the cluster --------------------------------------------------------------------------------------
# What setup_cluster.py installs. Pin a tag or commit (...BREAKDOWN.git@v0.1.0) so that results can be reproduced.
BREAKDOWN_INSTALL = 'git+https://github.com/TheodorePellegrinMPE/BREAKDOWN.git'
# Options of `breakdown analysedemonfolder`, they are recorded in every log and in full_log.txt.
BOND_OPTIONS = {'formation_distance': 1.63, 'break_distance': 2.5, 'min_lifetime': 1}
