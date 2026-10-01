from functools import lru_cache
import os
import matplotlib.pyplot as plt
from rdkit import Chem
from breakdown.plot_style import savePdfAndPng, textStyle
from breakdown.run_summary import loadSummary, resolveTimePerStep


@lru_cache(maxsize=None)
def checkElementalDissociation(smiles, heteroatom='O'):
    """Checks if any C or heteroatom atoms are in a fragment that is not the parent.
    The parent is the fragment with the most heavy atoms.

    Args:
        smiles (str): The SMILES string from the log.
        heteroatom (str): Element symbol of the heteroatom that is tracked.

    Returns:
        tuple: (boolCarbonInFragment, boolHeteroatomInFragment)
    """
    fragMols = [Chem.MolFromSmiles(f, sanitize=False) for f in smiles.split('.')]
    fragMols = [m for m in fragMols if m]
    if len(fragMols) < 2:
        return False, False
    parentIdx = max(range(len(fragMols)), key=lambda i: fragMols[i].GetNumHeavyAtoms())
    symbols = {a.GetSymbol() for i, m in enumerate(fragMols) if i != parentIdx for a in m.GetAtoms()}
    return 'C' in symbols, heteroatom in symbols


class ElementalStabilityChart:
    def __init__(self, chartName, saveFolder, heteroatom='O', timePerStepFs=5.0):
        """Initializes the ElementalStabilityChart object.

        Args:
            chartName (str): The name to be used for the output file.
            saveFolder (str): The directory where the plot will be saved.
            heteroatom (str): Element symbol of the heteroatom that is tracked, O or S.
            timePerStepFs (float): Simulated time between steps of the logs, in femtoseconds.
        """
        self.heteroatom = heteroatom
        self.chartName = chartName
        self.saveFolder = saveFolder
        self.timeStepToPs = timePerStepFs / 1000
        self.dataSeries = []

    def addDataFromFolders(self, folderPaths, abbreviations):
        """Loads atom count data from multiple simulation directories.

        Args:
            folderPaths (list): Paths to the directories containing full_log.txt.
            abbreviations (list): Display names for the molecules.
        """
        for folder, abbr in zip(folderPaths, abbreviations):
            data = self.loadDissociationData(folder)
            if data:
                self.dataSeries.append({
                    'abbr': abbr,
                    'cSeries': data['cIntact'],
                    'oSeries': data['oIntact']
                })

    def loadDissociationData(self, folder):
        """Gets, for every step, the number of runs in which the carbon atoms and the heteroatoms are
        all still in the parent. Read from run_summary.json if there is one, otherwise from full_log.txt,
        which is read line by line, only the runs at the start of it are used.

        Args:
            folder (str): The directory containing full_log.txt.

        Returns:
            dict: Lists of (time, count) tuples for Carbon and the heteroatom.
        """
        summary = loadSummary(folder)
        if summary and self.heteroatom in summary['elementalIntact']:
            intact = summary['elementalIntact']
            steps = sorted(int(step) for step in intact['runsAtStep'])
            return {'cIntact': [(s * self.timeStepToPs, intact['C'][str(s)]) for s in steps],
                    'oIntact': [(s * self.timeStepToPs, intact[self.heteroatom][str(s)]) for s in steps]}

        logPath = os.path.join(folder, 'full_log.txt')
        if not os.path.exists(logPath):
            return None

        carbonStats = {}
        heteroatomStats = {}
        totalRunsAtStep = {}

        with open(logPath, 'r', encoding='utf-8') as file:
            for line in file:
                if line.startswith('Automatically generated molecule names.'):
                    break
                elif line.startswith('At step'):
                    parts = line.split()
                    step = int(parts[2][:-1])
                    cLostNow, hLostNow = checkElementalDissociation(parts[-1], self.heteroatom)
                    carbonStats[step] = carbonStats.get(step, 0) + int(not cLostNow)
                    heteroatomStats[step] = heteroatomStats.get(step, 0) + int(not hLostNow)
                    totalRunsAtStep[step] = totalRunsAtStep.get(step, 0) + 1

        sortedSteps = sorted(totalRunsAtStep.keys())
        cSeries = [(s * self.timeStepToPs, carbonStats[s]) for s in sortedSteps]
        oSeries = [(s * self.timeStepToPs, heteroatomStats[s]) for s in sortedSteps]

        return {'cIntact': cSeries, 'oIntact': oSeries}

    def createStabilityPlot(self):
        """Generates a two-subplot figure for Carbon and heteroatom stability.

        Returns:
            str or None: Path of the pdf, None if there was no data.
        """
        if not self.dataSeries:
            return None
        with textStyle(), plt.rc_context({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'}):
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True, dpi=300)

            colors = plt.cm.Set1(range(9))
            lineStyles = ['-', '--', ':', '-.', (0, (3, 5, 1, 5)), (0, (5, 10))]

            for i, data in enumerate(self.dataSeries):
                time = [p[0] for p in data['cSeries']]

                ax1.plot(time, [p[1] for p in data['cSeries']],
                         label=data['abbr'], color=colors[i % 9],
                         linestyle=lineStyles[i % len(lineStyles)], linewidth=2)

                ax2.plot(time, [p[1] for p in data['oSeries']],
                         color=colors[i % 9], linestyle=lineStyles[i % len(lineStyles)],
                         linewidth=2)

            heteroatomName = 'Sulfur' if self.heteroatom == 'S' else 'Oxygen'
            self.applyFormatting(ax1, "Carbon Skeleton Intact ($N$)")
            self.applyFormatting(ax2, f"{heteroatomName} Content Intact ($N$)")
            ax2.set_xlabel('Time (ps)', fontsize=12, fontweight='bold')

            ax1.legend(fontsize=9, loc='lower left', frameon=True, edgecolor='black')

            plt.tight_layout()
            return self.savePlotToDisk()

    def applyFormatting(self, ax, ylabel):
        """Applies consistent professional styling."""
        ax.set_ylabel(ylabel, fontsize=11, fontweight='bold')
        ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.3)
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)

    def savePlotToDisk(self):
        """Saves the figure as a pdf and a png, returns the path of the pdf."""
        return savePdfAndPng(self.saveFolder, f"{self.chartName}_elemental_stability")


def generateElementalComparison(folders, names, title, outputDir, heteroatom='O', timePerStepFs=None):
    """Convenience function to generate the elemental comparison chart.

    Args:
        folders (list): Directories to analyze.
        names (list): Abbreviations for the legend.
        title (str): Output filename.
        outputDir (str): Output directory.
        heteroatom (str): Element symbol of the heteroatom that is tracked, O or S.
        timePerStepFs (float): Simulated time between steps of the logs, in femtoseconds.

    Returns:
        str or None: Path of the pdf.
    """
    chart = ElementalStabilityChart(title, outputDir, heteroatom, resolveTimePerStep(folders, timePerStepFs))
    chart.addDataFromFolders(folders, names)
    return chart.createStabilityPlot()
