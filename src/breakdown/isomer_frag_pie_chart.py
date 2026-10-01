import numpy as np
import matplotlib.pyplot as plt
from breakdown.plot_style import savePdfAndPng, textStyle

class IsomerFragPieChart:
    """A class to create a double nested pie chart showing the number of
    isomerizations on the outer ring and fragmentations on the inner.
    
    Attributes:
        fileName (str): The name of the file the chart will be saved as.
        folder (str): The name of the folder in which the chart will be saved.
        chartTitle (str): The title to be printed at the top of the chart.
        numNotIsomerizedOrFragmented (int): Number of molecules that did not isomerize or fragment.
        numIsomerizedNotFragmented (int): Number of molecules that only isomerized.
        numNotIsomerizedButFragmented (int): Number of molecules that only fragmented.
        numIsomerizedAndFragmented (int): Number of molecules that both isomerized and ended in more
            than one piece. A molecule that fragments and recombines later is not counted as fragmented.
    """

    def __init__(self, fileName, folder):
        """Constructs all necessary attributes of the IsomerFragPieChart object.

        Args:
            fileName (str): The name of the file the chart will be saved as.
            folder (str): The name of the folder in which the chart will be saved.
        """
        self.fileName = fileName
        self.folder = folder
        self.chartTitle = fileName.capitalize().replace('26_5', '26.5').replace('_', ' ')
        self.numNotIsomerizedOrFragmented = 0
        self.numIsomerizedNotFragmented = 0
        self.numNotIsomerizedButFragmented = 0
        self.numIsomerizedAndFragmented = 0

    def addNotIsomerizedOrFragmented(self):
        """Increment the number that were not isomerized or fragmented by one."""
        self.numNotIsomerizedOrFragmented += 1

    def addIsomerizedNotFragmented(self):
        """Increment the number that were only isomerized by one."""
        self.numIsomerizedNotFragmented += 1

    def addNotIsomerizedButFragmented(self):
        """Increment the number that were only fragmented by one."""
        self.numNotIsomerizedButFragmented += 1

    def addIsomerizedAndFragmented(self):
        """Increment the number that were isomerized and fragmented by one."""
        self.numIsomerizedAndFragmented += 1

    def createIsomerFragPieChart(self):
        """Creates a double nested pie chart based on the number of molecules
        in each of the 4 categories of the attributes. The chart is titled
        based on chartTitle and then saved at the fileName in the given folder.
        """
        total = (self.numNotIsomerizedOrFragmented + self.numIsomerizedNotFragmented +
                 self.numNotIsomerizedButFragmented + self.numIsomerizedAndFragmented)
        if total == 0:
            return None
        with textStyle():
            return self.drawChart()

    def drawChart(self):
        """Draws the chart, see createIsomerFragPieChart."""
        fig, ax = plt.subplots()
        size = 0.3
        fig.set_size_inches(7, 6)
        cmap = plt.get_cmap('tab20c')
        outer_colors = cmap([1, 2])
        inner_colors = cmap([6, 5, 5, 6])

        vals = self.getChartVals()

        ax.pie(vals.sum(axis=1), radius=1, colors=outer_colors,
           labels=[self.numIsomerizedNotFragmented +
               self.numIsomerizedAndFragmented,
               self.numNotIsomerizedOrFragmented +
               self.numNotIsomerizedButFragmented],
           autopct='', wedgeprops=dict(width=size, edgecolor='w'))

        ax.pie(vals.flatten(), radius=1-size, colors=inner_colors,
           labels=[self.numIsomerizedAndFragmented,
               self.numIsomerizedNotFragmented,
               self.numNotIsomerizedOrFragmented,
               self.numNotIsomerizedButFragmented],
           autopct='', wedgeprops=dict(width=size, edgecolor='w'))

        plt.title(self.chartTitle, fontsize=14)
        ax.legend(['Isomerization', 'No isomerization',
               'Fragmented at the end', 'Not fragmented at the end'], loc=(1, 0.4), fontsize=10)
        plt.subplots_adjust(left=0.1, right=0.7, top=0.9, bottom=0.1)
        return savePdfAndPng(f'{self.folder}/plots', f'{self.fileName}_frags_isos_pie')

    def getChartVals(self):
        """Makes a 2D numpy array of the values to be used to create the pie chart.

        Returns:
            np.ndarray: 2D numpy array.
        """
        vals = np.array([[self.numIsomerizedAndFragmented,
                          self.numIsomerizedNotFragmented],
                         [self.numNotIsomerizedOrFragmented,
                          self.numNotIsomerizedButFragmented]])
        return vals
