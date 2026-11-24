"""Write real MD geometries of a state, or of the steps around a transition, as .xyz files.

The `--coordinates` option of makechartsandstats builds force-field geometries from the SMILES strings. The MD frames
are better starting points for DFT: the geometry the molecule really had. For a transition A -> B the frames from
inside A to inside B are a ready made first guess for a NEB, with the first frame as reactant and the last as product.

The step number in a run log is the frame number in its .mol file. The bonds of every written frame are found again
with the settings the run was analysed with and compared with the state in the log, so a wrong .mol file or changed
settings are noticed.
"""
import os
import warnings
import networkx as nx
import numpy as np
from breakdown.auto_detect_isomerization import AutoDetectIsomerization
from breakdown.bond_adder import BondAdder
from breakdown.create_smiles_from_dataframe import makeSmiles
from breakdown.run_log import LOG_SUFFIX, listRunLogs, parseRun
from breakdown.run_summary import loadBlock, resolveTimePerStep

FRAMES_FOLDER = 'frames'


def getVisits(run):
    """Split a run into visits: stretches of steps in the same state.

    Args:
        run (Run): A run read from its log.

    Returns:
        list: (smiles, first step, last step) of every visit, in order.
    """
    visits = []
    for index, (step, _, smiles, _) in enumerate(run.steps):
        if visits and visits[-1][0] == smiles:
            continue
        if visits:
            visits[-1][2] = run.steps[index - 1].step
        visits.append([smiles, step, step])
    if visits:
        visits[-1][2] = run.steps[-1].step
    return [tuple(visit) for visit in visits]


def selectRunLogs(folder, names):
    """Names of the run logs that can contain the states, using the runSequences block when there is one.

    Args:
        folder (str): Folder with the run logs and full_log.txt.
        names (list): State names that must all appear in a run.

    Returns:
        list: Log file names.
    """
    sequences = loadBlock(folder, 'runSequences')
    logs = listRunLogs(folder)
    if sequences is None:
        return logs
    wanted = {f'{os.path.splitext(run)[0]}{LOG_SUFFIX}' for run, sequence in sequences.items()
              if all(name in sequence.values() for name in names)}
    return [log for log in logs if log in wanted]


def largestFragmentAtoms(nAtoms, atomTypes, pairs):
    """Atoms of the largest fragment: most heavy atoms, then most atoms.

    Args:
        nAtoms (int): Number of atoms.
        atomTypes (np.ndarray): Element of every atom.
        pairs (list): Bonded (first, second) atom pairs.

    Returns:
        set: Atom indices.
    """
    graph = nx.Graph()
    graph.add_nodes_from(range(nAtoms))
    graph.add_edges_from(pairs)
    return max(nx.connected_components(graph),
               key=lambda atoms: (sum(atomTypes[a] != 'H' for a in atoms), len(atoms), -min(atoms)))


class FrameExtractor:
    def __init__(self, folder, state, fromState=None, before=20, after=20, maxRuns=5, largestFragment=False,
                 outputFolder=None):
        """Initialize FrameExtractor.

        Args:
            folder (str): Folder with the .mol files, their run logs and full_log.txt.
            state (str): Name of the state, as in full_log.txt.
            fromState (str, optional): Write the transitions from this state into state instead of state alone.
            before (int): Steps before a transition to write, at most back to the start of the visit of fromState.
            after (int): Steps after a transition to write, at most to the end of the visit of state.
            maxRuns (int): Number of runs to take frames from, the runs are taken in file name order.
            largestFragment (bool): Only write the atoms of the largest fragment. For a transition, the atoms of the
                largest fragment at its first or at its last frame, so every frame has the same atoms.
            outputFolder (str, optional): Where to write, folder/frames by default.
        """
        self.folder = folder
        self.state = state
        self.fromState = fromState
        self.before = before
        self.after = after
        self.maxRuns = maxRuns
        self.largestFragment = largestFragment
        self.outputFolder = outputFolder or os.path.join(folder, FRAMES_FOLDER)

    def extract(self):
        """Write the frames.

        Returns:
            list: Paths of the files written.

        Raises:
            ValueError: Raised if the folder has no statistics or a state name is unknown.
        """
        isomerNames = loadBlock(self.folder, 'isomerNames')
        if isomerNames is None:
            raise ValueError(f'{self.folder} has no full_log.txt with statistics, run makechartsandstats first.')
        self.smilesOf = {name: smiles for smiles, name in isomerNames.items()}
        self.nameOf = isomerNames
        for name in filter(None, (self.state, self.fromState)):
            if name not in self.smilesOf:
                raise ValueError(f'{name} is not a state name in {self.folder}/full_log.txt, the names are listed '
                                 'under "Automatically generated molecule names."')
        self.timePerStepFs = resolveTimePerStep([self.folder])
        os.makedirs(self.outputFolder, exist_ok=True)
        written = []
        usedRuns = 0
        for logName in selectRunLogs(self.folder, [n for n in (self.state, self.fromState) if n]):
            if usedRuns >= self.maxRuns:
                break
            with open(os.path.join(self.folder, logName)) as file:
                run = parseRun(file.read())
            paths = self.extractFromRun(run)
            if paths:
                written += paths
                usedRuns += 1
        return written

    def findWanted(self, run):
        """The frames wanted from one run.

        Args:
            run (Run): The run.

        Returns:
            Tuple[list, int] or None: Steps to write, and the step of the transition (or the chosen frame of the
            state), None if the run does not have the state or transition.
        """
        visits = getVisits(run)
        target = self.smilesOf[self.state]
        if self.fromState is None:
            stays = [visit for visit in visits if visit[0] == target]
            if not stays:
                return None
            _, first, last = max(stays, key=lambda visit: visit[2] - visit[1])
            middle = (first + last) // 2
            return list(range(first, last + 1)), middle
        source = self.smilesOf[self.fromState]
        for previous, current in zip(visits, visits[1:]):
            if previous[0] == source and current[0] == target:
                start = max(previous[1], current[1] - self.before)
                end = min(current[2], current[1] + self.after)
                return list(range(start, end + 1)), current[1]
        return None

    def extractFromRun(self, run):
        """Write the frames of one run.

        Args:
            run (Run): The run.

        Returns:
            list: Paths written, empty if the run does not have the state or its .mol file is missing.
        """
        wanted = self.findWanted(run)
        if wanted is None:
            return []
        steps, eventStep = wanted
        molPath = os.path.join(self.folder, run.name)
        if not os.path.isfile(molPath):
            warnings.warn(f'{molPath} is missing, the frames of {run.name} are not written.')
            return []
        settings = run.settings or {}
        reader = AutoDetectIsomerization(self.folder, run.name)
        atomTypes, coords = reader.readTrajectory(molPath)
        coords = coords[:steps[-1] + 1]
        bondAdder = BondAdder(len(atomTypes), settings.get('formation_distance', 1.63),
                              settings.get('break_distance', 2.5))
        firstAtoms, secondAtoms, bonded = bondAdder.findBonds(atomTypes, coords)
        loggedSmiles = {step.step: step.smiles for step in run.steps}

        def smilesAt(step):
            pairs = np.flatnonzero(bonded[step])
            return makeSmiles(atomTypes, firstAtoms[pairs], secondAtoms[pairs])

        def pairsAt(step):
            pairs = np.flatnonzero(bonded[step])
            return list(zip(firstAtoms[pairs], secondAtoms[pairs]))

        runStem = os.path.splitext(run.name)[0]
        if self.fromState is None:
            # The middle of the longest visit, or the closest frame to it whose bonds give the logged state
            # (they can differ for min_lifetime above 1)
            frame = next((s for s in sorted(steps, key=lambda s: abs(s - eventStep))
                          if smilesAt(s) == loggedSmiles[s]), None)
            if frame is None:
                warnings.warn(f'No frame of {self.state} in {run.name} has the bonds of the logged state, check '
                              'that the run was analysed with the settings in its log. Writing the middle frame.')
                frame = eventStep
            atoms = self.selectAtoms(atomTypes, [pairsAt(frame)])
            path = os.path.join(self.outputFolder, f'{self.state}__{runStem}_step{frame}.xyz')
            self.writeXyz(path, atomTypes, coords, [frame], atoms, run.name, loggedSmiles)
            return [path]

        for end in (steps[0], steps[-1]):
            if smilesAt(end) != loggedSmiles[end]:
                warnings.warn(f'The bonds at step {end} of {run.name} do not give the logged state, check that the '
                              'run was analysed with the settings in its log.')
        atoms = self.selectAtoms(atomTypes, [pairsAt(steps[0]), pairsAt(steps[-1])])
        stem = os.path.join(self.outputFolder, f'{self.fromState}_to_{self.state}__{runStem}_step{eventStep}')
        paths = [f'{stem}.xyz', f'{stem}_reactant.xyz', f'{stem}_product.xyz']
        for path, frames in zip(paths, (steps, [steps[0]], [steps[-1]])):
            self.writeXyz(path, atomTypes, coords, frames, atoms, run.name, loggedSmiles)
        return paths

    def selectAtoms(self, atomTypes, pairLists):
        """Atoms to write: all, or the union of the largest fragments of the given frames.

        Args:
            atomTypes (np.ndarray): Element of every atom.
            pairLists (list): Bonded pairs of each frame.

        Returns:
            list: Sorted atom indices.
        """
        if not self.largestFragment:
            return list(range(len(atomTypes)))
        atoms = set()
        for pairs in pairLists:
            atoms |= largestFragmentAtoms(len(atomTypes), atomTypes, pairs)
        return sorted(atoms)

    def writeXyz(self, path, atomTypes, coords, frames, atoms, runName, loggedSmiles):
        """Write frames as a (multi-frame) xyz file, the comment line says the run, step, time and state.

        Args:
            path (str): File to write.
            atomTypes (np.ndarray): Element of every atom.
            coords (np.ndarray): Coordinates (steps, atoms, 3).
            frames (list): Steps to write.
            atoms (list): Atom indices to write, in this order.
            runName (str): .mol file of the run.
            loggedSmiles (dict): Step to the SMILES string in the log.
        """
        lines = []
        for frame in frames:
            state = self.nameOf.get(loggedSmiles[frame], '?')
            timePs = frame * self.timePerStepFs / 1000
            lines.append(str(len(atoms)))
            lines.append(f'{runName} step {frame} t={timePs:.3f} ps state {state}')
            for atom in atoms:
                x, y, z = coords[frame, atom]
                lines.append(f'{atomTypes[atom]:2} {x:14.6f} {y:14.6f} {z:14.6f}')
        with open(path, 'w') as file:
            file.write('\n'.join(lines) + '\n')
