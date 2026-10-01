import json
import os
import re
import warnings
from functools import lru_cache
import numpy as np
from breakdown.auto_detect_utils import classifyRun, countAtoms, countLetters, orderElements
from breakdown.isomer_frag_pie_chart import IsomerFragPieChart
from breakdown.fragments_pie_chart import FragmentsPieChart
from breakdown.fragments_bar_chart import FragmensBarChart
from breakdown.dissociation_time_chart import DissociationTimeChart
from breakdown.isomer_dicts import IsomerDicts
from breakdown.elemental_stability_chart import checkElementalDissociation
from breakdown.event_survival_chart import EVENT_NAMES, FIRST_EVENTS_HEADING, EventSurvivalChart
from breakdown.run_summary import getElementalIntactSeries, saveSummary

@lru_cache(maxsize=None)
def getMajorFragmentCounts(smiles):
    """Get the numbers of atoms of each element in the largest fragment of a molecule.
    The largest fragment is the one containing the most atoms.

    Args:
        smiles (str): SMILES string of the molecule, fragments are separated by periods.

    Returns:
        dict: Element symbols as keys and atom counts as values.
    """
    fragments = sorted(smiles.split('.'), key=countLetters, reverse=True)
    return countAtoms(fragments[0])


def inferMoleculeFromLogs(folder):
    """Get the number of carbons and the elements of the starting molecule from the
    first step of the first log in a folder.

    Args:
        folder (str): Path to folder with _log.txt files.

    Returns:
        Tuple[int, list]: Number of carbon atoms and list of element symbols, C and H first.

    Raises:
        ValueError: Raised if no log with a first step is found.
    """
    logFiles = sorted(file for file in os.listdir(folder)
                      if file.endswith('_log.txt') and not file.endswith('full_log.txt'))
    for logName in logFiles:
        with open(f'{folder}/{logName}') as file:
            for line in file:
                if line.startswith('At step'):
                    counts = {}
                    for fragment in line.split()[-1].split('.'):
                        for element, count in countAtoms(fragment).items():
                            counts[element] = counts.get(element, 0) + count
                    return counts.get('C', 0), orderElements(counts)
    raise ValueError(f'No _log.txt files with steps found in {folder}.')


class GenerateChartsStats:
    def __init__(self, folder, chartName, moleculeName, numCarbons, moleculeElements,
                 functionalGroup, drawIsomers, makeCoordinates, legendCarbons=None,
                 timePerStepFs=5.0):
        """Initialize the GenerateChartsStats class.

        Args:
            folder (str): Path to folder with log.txt files.
            chartName (str): Name to be used in the 3 plots, such as the molecule
                name, a suffix to specify which plot will be appended onto this name.
            moleculeName (str): Name of the molecule.
            numCarbons (int): Number of carbon atoms in the starting molecule.
            moleculeElements (list): Elements present in the starting molecule.
            functionalGroup (str or None): Original functional group to track.
            drawIsomers (bool): Draw every isomer.
            makeCoordinates (bool): Make a coordinates file for every isomer.
            legendCarbons (list, optional): Carbon counts to show in the dissociation legend.
            timePerStepFs (float): Simulated time between steps of the logs, in femtoseconds.
        """
        self.moleculeElements = moleculeElements
        self.heteroatom = 'S' if 'S' in moleculeElements else 'O'
        self.folder = folder
        self.chartName = chartName
        self.moleculeName = moleculeName
        self.allLogString = ''
        self.isoFragPieChart = IsomerFragPieChart(chartName, folder)
        self.fragmentPieChart = FragmentsPieChart(chartName, folder)
        self.fragmentBarChart = FragmensBarChart(chartName, folder)
        self.disTimeChart = DissociationTimeChart(chartName, folder, numCarbons, legendCarbons,
                                                  has_sulfur=self.heteroatom == 'S',
                                                  timePerStepFs=timePerStepFs,
                                                  heteroatoms=[e for e in moleculeElements if e not in ('C', 'H')])
        self.timePerStepFs = timePerStepFs
        self.drawIsomers = drawIsomers
        self.makeCoordinates = makeCoordinates
        self.functionalGroup = functionalGroup

    def autoAnalyseInFolder(self):
        """Analyses all log.txt files in a folder, outputs information in
        fullLog.txt and plots statistics in 3 plots.
        """
        saveName = 'full_log.txt'
        saveFile = f'{self.folder}/{saveName}'
        
        logFiles = sorted(file for file in os.listdir(self.folder)
                          if file.endswith('_log.txt') and not file.endswith('full_log.txt'))
        if not logFiles:
            print(f'No _log.txt files found in {self.folder}.')
            return
        logStrings = []
        for logName in logFiles:
            logString = self.addLogToCharts(logName)
            if logString is not None:
                logStrings.append(logString + '\n\n')
        self.allLogString = ''.join(logStrings)

        self.addStatsToString()
        self.saveStringAsFile(saveFile)
        self.saveSummary()
        self.createCharts()

    def addLogToCharts(self, logName):
        """Add log information to charts.

        Args:
            logName (str): Name of the log file.

        Returns:
            str or None: Log string, None if the run did not finish and was skipped.
        """
        logString = self.readLogFile(logName)
        if 'Final molecules present: ' not in logString or 'At step' not in logString:
            warnings.warn(f'{logName} is incomplete and is ignored.')
            return None
        self.addToDissociationChart(logString)
        self.addToIsomerizationCharts(logString)
        self.addToFragmentationCharts(logString)
        return logString

    def readLogFile(self, logName):
        """Read the content of a log file.

        Args:
            logName (str): Name of the log file.

        Returns:
            str: Log content
        """
        filePath = f'{self.folder}/{logName}'
        with open(filePath, 'r') as file:
            logString = file.read()
        return logString

    def addStatsToString(self):
        """Add statistical inforation about the amount of amount of time spent in each isomer, and
        the number of transitions between each pair of isomers to the log string. The blocks are
        also kept in self.sections, to be saved in the summary.
        """
        isomerDicts = self.handleIsomerDicts()
        self.isomerDicts = isomerDicts
        heteroatom = 'sulfur' if self.heteroatom == 'S' else 'oxygen'
        blocks = [
            ('Automatically generated molecule names.', 'isomerNames'),
            ('Total number of steps at which an isomer is present.', 'isomerStats'),
            ('Total number of runs in which an isomer appeared.', 'isomerRuns'),
            ('Number of runs at which each chemical fragment was present at the end.', 'endFragmentStats'),
            ('Total number of transitions to an isomer from an isomer.', 'isomerizationStats'),
            ('Net number of transitions to an isomer from an isomer.', 'netIsomerizationStats'),
            ('Number of transitions only including fragmentation.', 'fragmentationStats'),
            ('Chemical names of fragments of this isomerization state.', 'isomerChemicalNames'),
            ('Molecular masses of different observed isomers.', 'isomerMasses'),
        ]
        if self.functionalGroup is not None:
            blocks += [
                (f'Runs containing original functional group {self.functionalGroup} at step.', 'groupSteps'),
                (f'Runs with each {heteroatom} in each final group and/or fragment.', 'endingGroups'),
            ]
        blocks.append((FIRST_EVENTS_HEADING, 'firstEvents'))
        self.sections = {}
        for heading, attribute in blocks:
            if attribute == 'isomerNames':
                text = isomerDicts.getIsomerNames()
            else:
                text = isomerDicts.getFormattedStats(attribute)
            self.allLogString += heading + '\n' + text
            self.sections[heading] = json.loads(text)

    def saveSummary(self):
        """Save run_summary.json, which the comparison commands read instead of the huge full_log.txt."""
        trackedElements = [e for e in self.moleculeElements if e not in ('C', 'H')]
        summary = {
            'molecule': self.moleculeName,
            'heteroatom': self.heteroatom,
            'timePerStepFs': self.timePerStepFs,
            'numRuns': len(self.isomerDicts.runs),
            'sections': self.sections,
            'firstEvents': self.sections[FIRST_EVENTS_HEADING],
            'elementalIntact': getElementalIntactSeries(self.isomerDicts.runs, trackedElements,
                                                        checkElementalDissociation),
        }
        saveSummary(self.folder, summary)

    def createEventSurvivalChart(self):
        """Plot the time to the first change, isomerization and fragmentation across the runs."""
        firstEvents = self.sections[FIRST_EVENTS_HEADING]
        chart = EventSurvivalChart(self.chartName, f'{self.folder}/plots')
        for event in ('anyChange', 'isomerization', 'fragmentation'):
            chart.addCurve(EVENT_NAMES[event], firstEvents, event, self.timePerStepFs)
        chart.createPlot()

    def handleIsomerDicts(self):
        """Handle isomer statistical information dictionaries.

        Returns:
            IsomerDicts: IsomerDicts object
        """
        isomerDicts = IsomerDicts(self.allLogString, self.moleculeName, self.folder, self.functionalGroup,
                                  self.heteroatom)
        isomerDicts.createIsomerDicts()

        if self.drawIsomers:         
            isomerDicts.makeAllMoleculeDrawings()
        if self.makeCoordinates:
            isomerDicts.makeAllIsomerCoordinates()
        return isomerDicts

    def saveStringAsFile(self, saveFile):
        """Save a string to a file.

        Args:
            saveFile (str): Full path URL for the file to be saved at.
        """
        with open(saveFile, 'w') as file:
            file.write(self.allLogString)

    def addToIsomerizationCharts(self, logString):
        """Add to isoFragPieChart if a run included isomerization, was fragmented at the
        end, both or neither. A run that fragments and recombines is not counted as fragmented.

        Args:
            logString (str): String of SMILES and isomerization details for a run.
        """
        isomerized, _, fragmented = classifyRun(logString)
        if isomerized and fragmented:
            self.isoFragPieChart.addIsomerizedAndFragmented()
        elif isomerized:
            self.isoFragPieChart.addIsomerizedNotFragmented()
        elif fragmented:
            self.isoFragPieChart.addNotIsomerizedButFragmented()
        else:
            self.isoFragPieChart.addNotIsomerizedOrFragmented()

    def addToFragmentationCharts(self, logString):
        """Get fragments present at the end of the run and add them to
        fragmentBarChart and fragmentPieChart.

        Args:
            logString (str): String of SMILES and isomerization details for a run.
        """
        fragments = logString.split('Final molecules present: ')[-1].strip()
        fragments = fragments.split(', ')
        self.fragmentBarChart.addFragments(fragments)
        self.addToFragmentationPieChart(fragments)

    def addToFragmentationPieChart(self, fragments):
        """Add a fragment to fragmentPieChart based on the fragments list.

        If there is only 1 string in fragments, then the molecule has not
        fragmented. If there are 3 or more than there are multiple fragments.
        If there are 2 strings in fragments, the second, smaller fragment is
        recorded, fragments of C and H or C, H and O are stacked, this may be
        changed in the future.

        Args:
            fragments (list): List of strings, for molecules as chemical formulas.
        """
        if len(fragments) < 2:
            self.fragmentPieChart.addFragment('No fragmentation')
        elif len(fragments) > 2:
            self.fragmentPieChart.addFragment('Multiple fragments')
        else:
            fragment = fragments[1]
            if re.match(r'^C\d*H\d*O\d*$', fragment):
                self.fragmentPieChart.addFragment('CmHnOk')
            elif re.match(r'^C\d*H\d*S\d*$', fragment):
                self.fragmentPieChart.addFragment('CmHnSk')
            elif re.match(r'^C\d*H\d*$', fragment) and (fragment != 'C2H2'):
                self.fragmentPieChart.addFragment('CmHn')
            else:
                self.fragmentPieChart.addFragment(fragment)

    def addToDissociationChart(self, logString):
        """Add dissociation data to the dissociation time chart.

        Args:
            logString (str): String of SMILES and isomerization details for a run.
        """
        carbonCounts = []
        for line in logString.split('\n'):
            if line.startswith('At step'):
                majorFragment = getMajorFragmentCounts(line.split()[-1])
                numCarbon = majorFragment.get('C', 0)
                carbonCounts.append(numCarbon)
                self.disTimeChart.addToHydrogenArray(numCarbon, majorFragment.get('H', 0))
                for element in self.disTimeChart.heteroatoms:
                    self.disTimeChart.addToHeteroatomArray(element, numCarbon, majorFragment.get(element, 0))
        self.disTimeChart.addRun(np.array(carbonCounts))

    def createCharts(self):
        """Run method to make and save the plots from each chart object.
        """
        self.isoFragPieChart.createIsomerFragPieChart()
        self.fragmentBarChart.createFragmentsBarChart()
        self.fragmentPieChart.createFragmentsPieChart()
        self.disTimeChart.createDissociationTimeChart()
        self.createEventSurvivalChart()
