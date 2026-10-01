import json
import re
from collections import Counter
from breakdown.isomer_drawer import IsomerDrawer
from breakdown.coordinates_file_maker import CoordinatesFileMaker
from breakdown.auto_detect_utils import (EVENT_LABELS, calculateNominalMass, countAtoms,
                                             getFragmentNames, processSmileToFragment)

OXYGEN_GROUPS = ['hydroxyl', 'ketone', 'ketene', 'ether', 'aldehyde', 'carboxyl', 'epoxide']
SULFUR_GROUPS = ['thiol', 'thioketone', 'thioaldehyde', 'thioether', 'thiirane', 'disulfide']

class IsomerDicts:
    def __init__(self, allLogString, moleculeName, saveFolder, groupName, heteroatom='O'):
        """Initialize IsomerDicts.

        Args:
            allLogString (str): Complete log string.
            moleculeName (str): Name of the molecule.
            saveFolder (str): Path to folder where files will be saved.
            groupName (str): Molecule functional group to track.
            heteroatom (str): Element of the groups that are tracked at the end, O or S.
        """
        self.allLogString = allLogString
        self.moleculeName = moleculeName
        self.saveFolder = saveFolder
        self.groupName = groupName
        self.runs = self.parseRuns(allLogString)
        self.heteroatom = heteroatom

    @staticmethod
    def parseRuns(allLogString):
        """Split the log string into runs, once, so that every statistic can use the result.

        Args:
            allLogString (str): Complete log string.

        Returns:
            list: One dictionary per run, with its file name under 'name', the list of its
            steps as tuples of (step, event, smiles, analysis line) under 'steps', and the
            list of chemical formulas present at the end under 'final' (None if missing).
        """
        runs = []
        run = None
        lines = allLogString.splitlines()
        for index, line in enumerate(lines):
            if line.startswith('In file'):
                run = {'name': line.split(': ', 1)[1].strip(), 'steps': [], 'final': None}
                runs.append(run)
            elif line.startswith('At step'):
                if run is None or run['final'] is not None:
                    run = {'name': f'run_{len(runs) + 1}', 'steps': [], 'final': None}
                    runs.append(run)
                parts = line.split()
                event = parts[3] if len(parts) == 5 and parts[3] in EVENT_LABELS else ''
                analysis = lines[index + 1] if event and index + 1 < len(lines) else ''
                run['steps'].append((int(parts[2][:-1]), event, parts[-1], analysis))
            elif line.startswith('Final molecules present: ') and run is not None:
                run['final'] = line.split(': ', 1)[1].split(', ')
        return runs

    def createIsomerDicts(self):
        """Create isomer statistics, isomerization statistics, and isomer names."""
        self.isomerStats = self.createIsomerStatistics()
        self.isomerizationStats = self.createIsomerizationStats()
        self.netIsomerizationStats = self.createNetIsomerizationStats()
        self.fragmentationStats = self.createFragmenationStats()
        self.isomerNames = self.createIsomerNames()
        self.isomerRuns = self.createIsomerRunsStats()
        self.endFragments = self.createEndFragmentsStats()
        self.isomerChemicalNames = self.createIsomerChemicalNames()
        self.isomerMasses = self.createIsomerMassDict()
        self.firstEvents = self.createFirstEventStats()
        if self.groupName is not None:
            self.groupSteps = self.createFunctionalGroupStats()
            self.endingGroups = self.createGroupEndingStats()

    def createIsomerStatistics(self):
        """Make a dictionary with all SMILES strings in the file as keys and
        the count of their number of occurrences in the file as values.

        Returns:
            dict: Dictionary with SMILES strings as keys and their
            counts in the file as values.
        """
        isomerStats = Counter(smiles for run in self.runs for _, _, smiles, _ in run['steps'])
        isomerStats = dict(sorted(isomerStats.items(), key=lambda item: item[1], reverse=True))
        return isomerStats

    def createIsomerizationStats(self):
        """Create isomerization statistics. A dictionary where the key is the SMILES string of
        an isomer, and the values are dictionaries where the the keys are isomers that transition
        into the isomer, and the values are the number of times this transition happens.

        Returns:
            dict: Dictionary representing isomerization statistics.
        """
        isomerizationStats = {}
        for run in self.runs:
            lastSmiles = ''
            for _, event, smiles, _ in run['steps']:
                if event:
                    transitions = isomerizationStats.setdefault(smiles, {})
                    transitions[lastSmiles] = transitions.get(lastSmiles, 0) + 1
                lastSmiles = smiles
        return isomerizationStats

    def createNetIsomerizationStats(self):
        """Create a dictionary with isomer names as keys and dictionaries of net number of transitions
        as values.

        Returns:
            dict: Dictionary with isomer names as keys and dictionaries of net net number of transitions
            as values.
        """
        netIsomerizationStats = {}
        for smiles, values in self.isomerizationStats.items():
            netValues = {}
            for otherSmiles, count in values.items():
                if otherSmiles in self.isomerizationStats and smiles in self.isomerizationStats[otherSmiles]:
                    netCount = count - self.isomerizationStats[otherSmiles][smiles]
                else:
                    netCount = count
                netValues[otherSmiles] = netCount
            netValues = {k: v for k, v in netValues.items() if v > 0}
            if netValues:
                netIsomerizationStats[smiles] = netValues
        return netIsomerizationStats
    
    def createFragmenationStats(self):
        """Create a dictionary with isomer names as keys and the net number of times the isomer transitioned
        into a fragment with a different number of atoms as values.

        Returns:
            dict: Dictionary with isomer names as keys and the net number of times the isomer transitioned into
            a fragment with a different number of atoms as values.
        """
        fragmentationStats = {}
        for smiles, values in self.netIsomerizationStats.items():
            for otherSmiles, count in values.items():
                if smiles.count('.') != otherSmiles.count('.'):
                    if smiles not in fragmentationStats:
                        fragmentationStats[smiles] = {}
                    fragmentationStats[smiles][otherSmiles] = count
        return fragmentationStats
    
    def createIsomerNames(self):
        """Create isomer names.

        Returns:
            dict: Dictionary mapping SMILES strings to isomer names.
        """
        isomerNames = {}
        nameCounts = {}
        for smiles in self.isomerStats.keys():
            isomerNames[smiles] = self.makeIsomerName(smiles, nameCounts)
        return isomerNames

    def makeIsomerName(self, smiles, nameCounts):
        """Make isomer name. The name is generated to include a letter indicating the number of lost
        carbon atoms, and a number n indicating that the isomer is the nth most common isomer of that type.

        Args:
            smiles (str): SMILES string.
            nameCounts (dict): Number of names already made for each name without a number.

        Returns:
            str: Isomer name.
        """
        letter = self.determineCarbonLossLetter(smiles)
        name = f'{self.moleculeName}_{letter}'
        nameCounts[name] = nameCounts.get(name, 0) + 1
        return f'{name}{nameCounts[name]}'

    def determineCarbonLossLetter(self, smiles):
        """Determine the letter for carbon loss.

        Args:
            smiles (str): SMILES string.

        Returns:
            str: Letter representing carbon loss, A for none, B for one carbon
            and so on. After Z it is Z followed by the number of carbons lost.
        """
        carbonCounts = [countAtoms(fragment).get('C', 0) for fragment in smiles.split('.')]
        carbonLoss = sum(carbonCounts) - max(carbonCounts)
        if carbonLoss > 25:
            return f'Z{carbonLoss}'
        return chr(ord('A') + carbonLoss)

    def createFirstEventStats(self):
        """For every run, the first step of each type of event, None if it never happened.

        Returns:
            dict: File name of the run to a dictionary with the keys isomerization,
            fragmentation, recombination, anyChange, and lastStep, the last step of the run.
        """
        keys = {'Isomerization': 'isomerization', 'Fragmentation': 'fragmentation',
                'Recombination': 'recombination'}
        firstEvents = {}
        for run in self.runs:
            if not run['steps']:
                continue
            record = dict.fromkeys(list(keys.values()) + ['anyChange'])
            for step, event, _, _ in run['steps']:
                if event:
                    for key in (keys[event], 'anyChange'):
                        if record[key] is None:
                            record[key] = step
            record['lastStep'] = run['steps'][-1][0]
            firstEvents[run['name']] = record
        return firstEvents

    def createIsomerRunsStats(self):
        """Create a dictionary with all SMILES strings in the file as keys and
        the number of runs in which the isomer occurs as values.

        Returns:
            dict: Dictionary with SMILES strings as keys and their
            counts in the file as values.
        """
        isomerRuns = Counter()
        for run in self.runs:
            isomerRuns.update(dict.fromkeys(smiles for _, _, smiles, _ in run['steps']).keys())
        isomerRuns = dict(sorted(isomerRuns.items(), key=lambda item: item[1], reverse=True))
        return isomerRuns

    def createEndFragmentsStats(self):
        """Create a dictionary with all fragments present at the last step of the run as keys
        and the number of runs in which the fragment is present at the end as values.

        Returns:
            dict: Dictionary with fragment names as keys and their counts in the file as values."""
        endFragmentStats = Counter()
        for run in self.runs:
            if run['final'] is not None:
                endFragmentStats.update(run['final'])
        return dict(endFragmentStats)

    def createIsomerChemicalNames(self):
        """Create a dictionary with isomer names as keys and lists of isomer chemical names as values.

        Returns:
            dict: Dictionary with isomer names as keys and lists of isomer chemical names as values.
        """
        isomerChemicalNames = {}
        for smiles, name in self.isomerNames.items():
            for fragmentSmiles in smiles.split('.'):
                fragment = processSmileToFragment(fragmentSmiles)
                if name not in isomerChemicalNames:
                    isomerChemicalNames[name] = [fragment]
                else:
                    if fragment not in isomerChemicalNames[name]:
                        isomerChemicalNames[name].append(fragment)
        return isomerChemicalNames
    
    def createIsomerMassDict(self):
        """Create a dictionary with molecular mass as keys and isomer names as values.

        Returns:
            dict: Dictionary with molecular mass as keys and isomer chemical names as values.
        """
        massIsomerDict = {}
        for smiles in self.isomerNames:
            for fragmentSmiles in smiles.split('.'):
                fragment = processSmileToFragment(fragmentSmiles)
                formulas = massIsomerDict.setdefault(calculateNominalMass(fragmentSmiles), [])
                if fragment not in formulas:
                    formulas.append(fragment)
        massIsomerDict = dict(sorted(massIsomerDict.items(), key=lambda item: item[0], reverse=True))
        return massIsomerDict

    def getIsomerNames(self):
        """Get isomer names as a formatted string.

        Returns:
            str: Formatted string of isomer names.
        """
        isomerNames = json.dumps(self.isomerNames, indent=4) + '\n'
        return isomerNames

    def getFormattedStats(self, dictionaryName):
        """Get formatted statistics as a string, with SMILES strings replaced by isomer names.

        Args:
            dictionaryName (str): Name of the attribute holding the dictionary of statistics.

        Returns:
            str: Formatted string of the statistics.
        """
        names = {'endFragmentStats': 'endFragments'}
        stats = getattr(self, names.get(dictionaryName, dictionaryName))
        statsString = json.dumps(self.replaceWithMolNames(stats), indent=4) + '\n'
        return statsString

    def replaceWithMolNames(self, stats):
        """Replace SMILES strings that are dictionary keys with isomer names, in the
        dictionary and in the dictionaries inside it.

        Args:
            stats (dict): Dictionary of statistics.

        Returns:
            dict: Copy of the dictionary using isomer names as keys.
        """
        renamed = {}
        for key, value in stats.items():
            if isinstance(value, dict):
                value = self.replaceWithMolNames(value)
            renamed[self.isomerNames.get(key, key)] = value
        return renamed

    def makeAllMoleculeDrawings(self):
        """Make drawings of all molecules."""
        drawingsFolder = f'{self.saveFolder}/drawings'
        isomerDrawer = IsomerDrawer(drawingsFolder)
        for smiles, name in self.isomerNames.items():
            try:
                isomerDrawer.drawAndSaveSmiles(smiles, name)
            except Exception as e:
                print(f'Error drawing {name}: {e}')
    
    def makeAllIsomerCoordinates(self):
        """Generate coordinates for all isomers and save them to files.
        """
        coordinatesFolder = f'{self.saveFolder}/coordinates'
        fileMaker = CoordinatesFileMaker(coordinatesFolder)

        for smiles, name in self.isomerNames.items():
            self.makeIsomerCoordinates(fileMaker, smiles, name)

    def makeIsomerCoordinates(self, fileMaker, smiles, name):
        """Generate coordinates for a specific isomer and save them to a file.

        Args:
            fileMaker (CoordinatesFileMaker): An instance of CoordinatesFileMaker
                responsible for creating and saving coordinates files.
            smiles (str): The SMILES representation of the isomer.
            name (str): The name associated with the isomer.
        """
        fragments, isomerNames = getFragmentNames(smiles, name)
        for i in range(len(fragments)):
            fileMaker.makeXYZForSmiles(fragments[i], isomerNames[i])

    def createFunctionalGroupStats(self):
        """
        Creates a dictionary mapping each timestep to the number of runs
        containing the specified functional group, supporting count thresholds
        via the group_name_# syntax.

        A group is present from the start of a run, until the analysis of a molecule
        after an isomerization or fragmentation does not have it in enough numbers.

        Returns:
            dict: {step: number_of_runs_with_group}
        """
        targetGroup = self.groupName.lower()
        minCount = 1
        if "_" in targetGroup:
            parts = targetGroup.split("_")
            if parts[-1].isdigit():
                minCount = int(parts[-1])
                targetGroup = "_".join(parts[:-1])
        targetPattern = re.compile(rf"\b{re.escape(targetGroup)}")
        countPattern = re.compile(rf"\b{re.escape(targetGroup)} count: (\d+)")

        globalStepStats = {}
        for run in self.runs:
            isPresent = True
            for step, event, _, analysis in run['steps']:
                if event:
                    analysisLine = analysis.lower()
                    if targetPattern.search(analysisLine):
                        countMatch = countPattern.search(analysisLine)
                        if countMatch:
                            isPresent = int(countMatch.group(1)) >= minCount
                        else:
                            isPresent = minCount <= 1
                    else:
                        isPresent = False
                globalStepStats[step] = globalStepStats.get(step, 0) + int(isPresent)
        return globalStepStats

    def createGroupEndingStats(self):
        """
        Calculates the number of runs ending with specific combinations of
        oxygen (or sulfur) containing groups and fragments.

        The groups are the ones counted in the analysis of the molecule after the last
        isomerization or fragmentation of a run. A run without either keeps its
        original group. A run whose last analysis has no groups and no fragments with
        the heteroatom is counted as 'other'.

        Returns:
            dict: {heteroatom_combination: number_of_runs}
        """
        targetGroups = SULFUR_GROUPS if self.heteroatom == 'S' else OXYGEN_GROUPS
        endingStats = {}
        for run in self.runs:
            if not run['steps']:
                continue
            lastAnalysis = next((analysis for _, event, _, analysis in reversed(run['steps'])
                                 if event), None)
            finalGroups = []
            for group in targetGroups:
                pattern = rf"\b{group}(?:\s+group)?\s+count:\s+(\d+)"
                countMatch = re.search(pattern, lastAnalysis or '', re.IGNORECASE)
                if countMatch:
                    finalGroups.extend([group] * int(countMatch.group(1)))
            allFeatures = sorted(finalGroups + self.getHeteroatomFragments(run['final'] or []))
            if lastAnalysis is None:
                combinedKey = self.groupName
            elif not allFeatures:
                combinedKey = 'other'
            else:
                combinedKey = " / ".join(allFeatures)
            endingStats[combinedKey] = endingStats.get(combinedKey, 0) + 1
        return endingStats

    def getHeteroatomFragments(self, fragments):
        """Get the fragments, apart from the parent, that contain the heteroatom.

        Args:
            fragments (list): Chemical formulas of the molecules present at the end of a run.

        Returns:
            list: Formulas of the fragments with the heteroatom, the parent is
            the one with the most carbon atoms.
        """
        carbonCounts = [self.getElementCount(formula, 'C') for formula in fragments]
        parentIndex = carbonCounts.index(max(carbonCounts)) if carbonCounts and max(carbonCounts) > 0 else -1
        return [formula for index, formula in enumerate(fragments)
                if index != parentIndex and self.getElementCount(formula, self.heteroatom) > 0]

    @staticmethod
    def getElementCount(formula, element):
        """Get the number of atoms of an element in a chemical formula like C6H5O2.

        Args:
            formula (str): Chemical formula.
            element (str): Element symbol.

        Returns:
            int: Number of atoms of the element.
        """
        match = re.search(rf'{element}(?![a-z])(\d*)', formula)
        if not match:
            return 0
        return int(match.group(1)) if match.group(1) else 1
