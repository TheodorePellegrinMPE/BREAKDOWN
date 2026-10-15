import os
import warnings
from breakdown.auto_detect_utils import countAtoms, orderElements
from breakdown.isomer_frag_pie_chart import IsomerFragPieChart
from breakdown.fragments_pie_chart import FragmentsPieChart
from breakdown.fragments_bar_chart import FragmensBarChart
from breakdown.dissociation_time_chart import DissociationTimeChart
from breakdown.isomer_dicts import IsomerDicts
from breakdown.elemental_stability_chart import checkElementalDissociation
from breakdown.event_survival_chart import EVENT_NAMES, EventSurvivalChart
from breakdown.run_log import FULL_LOG_FILE, isComplete, iterSteps, listRunLogs, parseRuns
from breakdown.run_settings import collectSettings
from breakdown.run_summary import formatBlocks, getElementalIntactSeries, saveSummary


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
    for logName in listRunLogs(folder):
        with open(f'{folder}/{logName}') as file:
            for _, _, smiles in iterSteps(file):
                counts = {}
                for fragment in smiles.split('.'):
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
                                                  timePerStepFs=timePerStepFs,
                                                  heteroatoms=[e for e in moleculeElements if e not in ('C', 'H')])
        self.timePerStepFs = timePerStepFs
        self.statisticsOptions = {'chart_name': chartName, 'molecule_name': moleculeName, 'num_carbons': numCarbons,
                                  'elements': list(moleculeElements), 'func_group': functionalGroup,
                                  'legend_carbons': legendCarbons, 'time_per_step_fs': timePerStepFs,
                                  'draw': drawIsomers, 'coordinates': makeCoordinates}
        self.drawIsomers = drawIsomers
        self.makeCoordinates = makeCoordinates
        self.functionalGroup = functionalGroup

    def autoAnalyseInFolder(self):
        """Analyses all log.txt files in a folder, outputs information in
        fullLog.txt and plots statistics in 3 plots.
        """
        saveFile = f'{self.folder}/{FULL_LOG_FILE}'
        logFiles = listRunLogs(self.folder)
        if not logFiles:
            print(f'No _log.txt files found in {self.folder}.')
            return
        logStrings = []
        for logName in logFiles:
            logString = self.readLogFile(logName)
            if isComplete(logString):
                logStrings.append(logString + '\n\n')
            else:
                warnings.warn(f'{logName} is incomplete and is ignored.')
        self.allLogString = ''.join(logStrings)
        self.runs = parseRuns(self.allLogString)

        self.addStatsToString()
        self.saveStringAsFile(saveFile)
        self.saveSummary()
        self.createCharts()

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
        self.isomerDicts = self.handleIsomerDicts()
        self.settings = collectSettings([run.settings for run in self.runs], self.statisticsOptions)
        self.reportSettings()
        self.blocks = self.isomerDicts.getNamedStatistics()
        self.blocks['settings'] = self.settings
        text, self.sections = formatBlocks(self.blocks, self.functionalGroup, self.heteroatom)
        self.allLogString += text

    def reportSettings(self):
        """Warn if the runs were made with different settings or without recorded settings."""
        settings = self.settings
        if len(settings['analysis']) > 1:
            warnings.warn(f"The runs in {self.folder} were analysed with {len(settings['analysis'])} different sets of "
                          'settings (bond distances, minimum lifetime), the statistics mix them. See the settings '
                          'block in full_log.txt.')
        if settings['runsWithoutRecordedSettings']:
            warnings.warn(f"{settings['runsWithoutRecordedSettings']} logs have no recorded settings "
                          '(made by an older version).')

    def saveSummary(self):
        """Save run_summary.json, which the comparison commands read instead of the huge full_log.txt."""
        trackedElements = [e for e in self.moleculeElements if e not in ('C', 'H')]
        summary = {
            'molecule': self.moleculeName,
            'heteroatom': self.heteroatom,
            'timePerStepFs': self.timePerStepFs,
            'numRuns': len(self.runs),
            'sections': self.sections,
            'settings': self.settings,
            'firstEvents': self.blocks['firstEvents'],
            'elementalIntact': getElementalIntactSeries(self.runs, trackedElements,
                                                        checkElementalDissociation),
        }
        saveSummary(self.folder, summary)

    def createEventSurvivalChart(self):
        """Plot the time to the first change, isomerization and fragmentation across the runs."""
        firstEvents = self.blocks['firstEvents']
        chart = EventSurvivalChart(self.chartName, f'{self.folder}/plots')
        for event in ('anyChange', 'isomerization', 'fragmentation'):
            chart.addCurve(EVENT_NAMES[event], firstEvents, event, self.timePerStepFs)
        chart.createPlot()

    def handleIsomerDicts(self):
        """Handle isomer statistical information dictionaries.

        Returns:
            IsomerDicts: IsomerDicts object
        """
        isomerDicts = IsomerDicts(self.runs, self.moleculeName, self.folder, self.functionalGroup,
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

    def createCharts(self):
        """Run method to make and save the plots from each chart object.
        """
        for run in self.runs:
            for chart in (self.disTimeChart, self.isoFragPieChart, self.fragmentBarChart, self.fragmentPieChart):
                chart.addRun(run)
        self.isoFragPieChart.createIsomerFragPieChart()
        self.fragmentBarChart.createFragmentsBarChart()
        self.fragmentPieChart.createFragmentsPieChart()
        self.disTimeChart.createDissociationTimeChart()
        self.createEventSurvivalChart()
