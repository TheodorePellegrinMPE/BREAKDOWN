from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from breakdown.plot_style import rangeFormulaLabel, savePdfAndPng, textStyle
from breakdown.auto_detect_utils import ensureFolderExists


class DissociationTimeChart:
    def __init__(self, fileName, folder, initialCarbon, legendCarbons=None, has_sulfur=False,
                 timePerStepFs=5.0, heteroatoms=None):
        """Constructs all necessary attributes of the DissociationTimeChart object.

        Args:
            fileName (str): The name of the file the chart will be saved as.
            folder (str): The name of the folder in which the chart will be saved.
            initialCarbon (int): The initial number of carbon atoms.
            legendCarbons (list, optional): List of integers for carbons to show in legend.
            has_sulfur (bool, optional): If True, track Sulfur ('S') instead of Oxygen ('O') in the legend,
                used if heteroatoms is not given.
            timePerStepFs (float): Simulated time between the steps of the runs, in femtoseconds.
                The default is 0.5 fs time steps with output every 10 steps.
            heteroatoms (list, optional): Elements other than C and H to show in the legend formulas.
        """
        self.fileName = fileName
        self.folder = folder
        self.chartTitle = fileName.capitalize().replace('26_5', '26.5').replace('_', ' ')
        self.initialCarbon = initialCarbon
        self.stepSize = timePerStepFs / 1000
        self.runCarbons = []
        self.has_sulfur = has_sulfur
        self.heteroatoms = list(heteroatoms) if heteroatoms else (['S'] if has_sulfur else ['O'])
        # For every number of carbon atoms in the largest fragment, the numbers of H and of the other atoms seen
        self.hydrogenArray = defaultdict(list)
        self.heteroatomArrays = {element: defaultdict(list) for element in self.heteroatoms}
        self.legendCarbons = legendCarbons

    def addRun(self, carbonCounts):
        """Add the number of carbon atoms in the largest fragment at every step of a run.

        Args:
            carbonCounts (np.ndarray): Number of carbons in the largest fragment at each step.
        """
        self.runCarbons.append(np.asarray(carbonCounts, dtype=int))

    def getDissociationArray(self):
        """Count, at every step, the percentage of runs with each carbon count in the
        largest fragment. A run that has stopped is assumed to stay in its last state.

        Returns:
            np.ndarray: Array with shape (steps, largest carbon count + 1).
        """
        numSteps = max(len(run) for run in self.runCarbons)
        numCarbons = max(self.initialCarbon, max(run.max() for run in self.runCarbons)) + 1
        counts = np.zeros([numSteps, numCarbons])
        stepIndices = np.arange(numSteps)
        for run in self.runCarbons:
            padded = np.pad(run, (0, numSteps - len(run)), mode='edge')
            counts[stepIndices, padded] += 1
        return 100 * counts / len(self.runCarbons)

    def addToHydrogenArray(self, numCarbon, numHydrogen):
        """Adds a number of hydrogens to the list for a number of carbons in hydrogenArray if that number is not already present.

        Args:
            numCarbon (int): Integer number of carbons.
            numHydrogen (int): Integer number of hydrogens.
        """
        if not numHydrogen in self.hydrogenArray[numCarbon]:
            self.hydrogenArray[numCarbon].append(numHydrogen)

    def addToHeteroatomArray(self, element, numCarbon, numAtoms):
        """Adds a number of atoms of an element to the list for a number of carbons if that number is not already present.

        Args:
            element (str): Element symbol, one of the heteroatoms of the chart.
            numCarbon (int): Integer number of carbons.
            numAtoms (int): Integer number of atoms of the element.
        """
        if element in self.heteroatomArrays and numAtoms not in self.heteroatomArrays[element][numCarbon]:
            self.heteroatomArrays[element][numCarbon].append(numAtoms)

    def addToOxygenArray(self, numCarbon, numOxygen):
        """Adds a number of oxygens to the list for a number of carbons if that number is not already present."""
        self.addToHeteroatomArray('O', numCarbon, numOxygen)

    def addToSulfurArray(self, numCarbon, numSulfur):
        """Adds a number of sulfurs to the list for a number of carbons if that number is not already present."""
        self.addToHeteroatomArray('S', numCarbon, numSulfur)

    def createDissociationTimeChart(self):
        """Creates a line graph based on the size of parent atom over time.

        Returns:
            str: Path of the pdf, the png is saved next to it.
        """
        dissociationArray = self.getDissociationArray()
        steps, numCarbons = dissociationArray.shape
        with textStyle():
            cmap = plt.get_cmap('tab10')
            lineStyles = ['-', '--', '-.', ':']
            plt.figure(figsize=(9, 5))
            legend_lines = []
            legend_labels = []

            for m in reversed(range(numCarbons)):
                carbons = dissociationArray[:, m]
                if np.any(carbons != 0):
                    color = cmap(m % cmap.N)
                    line_style = lineStyles[m % len(lineStyles)]
                    line, = plt.plot(np.multiply(range(steps), self.stepSize),
                                     carbons, marker='', label=self.makeLegendString(m),
                                     color=color, linestyle=line_style, linewidth=3)
                    if self.legendCarbons is None or m in self.legendCarbons:
                        legend_lines.append(line)
                        legend_labels.append(self.makeLegendString(m))

            plt.xlabel('Time, t (ps)', fontsize=13)
            plt.ylabel('% ions', fontsize=13)
            plt.legend(legend_lines, legend_labels, loc=(0.6, 0.5), fontsize=11, handlelength=2.5)
            plt.title(self.chartTitle, fontsize=15)
            plt.xticks(fontsize=11)
            plt.yticks(fontsize=11)
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15)
            return savePdfAndPng(f'{self.folder}/plots', f'{self.fileName}_dissociation')

    def makeLegendString(self, numCarbon):
        """Make a string to be used in the dissociation chart's legend, a formula such as C_6H_{6-7}O_{0-1} where the
        subscripts are the numbers (or the range of numbers) of atoms seen in the largest fragment of that size.

        Args:
            numCarbon (int): Index of the line to be given a legend.

        Returns:
            str: String to be used in the legend.
        """
        others = {element: arrays[numCarbon] for element, arrays in self.heteroatomArrays.items()}
        return rangeFormulaLabel(numCarbon, self.hydrogenArray[numCarbon], others)
