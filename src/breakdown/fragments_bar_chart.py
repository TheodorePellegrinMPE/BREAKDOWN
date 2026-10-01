import matplotlib.pyplot as plt
from breakdown.plot_style import savePdfAndPng, textStyle

class FragmensBarChart:
    """A class to create a bar chart showing the number of instances of each
    fragment seen at the end of each simulation. All pieces including the
    parent molecule are considered if they are present at the end.

    Attributes:
        fileName (str): The name of the file the chart will be saved as.
        folder (str): The name of the folder in which the chart will be saved.
        chartTitle (str): The title to be printed at the top of the chart.
        fragmentCounts (dict): A dictionary with the names of fragments as its
            keys and the count of each fragment as its values.
    """
    
    def __init__(self, fileName, folder):
        """Constructs all necessary attributes of the FragmentsBarChart object.

        Args:
            fileName (str): The name of the file the chart will be saved as.
            folder (str): The name of the folder in which the chart will be saved.
        """
        self.fileName = fileName
        self.folder = folder
        self.chartTitle = fileName.capitalize().replace('26_5', '26.5').replace('_', ' ')
        self.fragmentCounts = {}

    def addFragments(self, fragments):
        """Increments the value of a fragment type in the fragmentCounts
        dictionary, if the type is not a key yet it is added and set to one.

        Args:
            fragments (list): List of names of fragments present at the last
                step of an individual simulation.
        """
        for fragment in fragments:
            if fragment in self.fragmentCounts:
                self.fragmentCounts[fragment] += 1
            else:
                self.fragmentCounts[fragment] = 1

    def createFragmentsBarChart(self):
        """Creates a bar chart based on the count of each fragment type.
        The chart is titled based on chartTitle and then saved at the
        fileName in the given folder.
        """
        if not self.fragmentCounts:
            return None
        fragmentNameList = sorted(list(self.fragmentCounts.keys()))
        fragmentCountList = [self.fragmentCounts[fragment] for fragment in fragmentNameList]
        with textStyle():
            plt.figure(figsize=(10, 5))
            plt.bar(fragmentNameList, fragmentCountList)
            plt.xlabel('Fragment', fontsize=12)
            plt.ylabel('Count', fontsize=12)
            plt.title(self.chartTitle, fontsize=14)
            plt.xticks(rotation=90, fontsize=10)
            plt.yticks(fontsize=10)
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.25)
            return savePdfAndPng(f'{self.folder}/plots', f'{self.fileName}_fragmentations_bar')
