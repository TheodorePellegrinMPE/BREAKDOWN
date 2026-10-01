import matplotlib.pyplot as plt
from breakdown.plot_style import savePdfAndPng, textStyle

class FragmentsPieChart:
    """A class to create a pie chart showing the number of instances of each
    fragment seen in the simulations. The largest piece is not considered
    a fragment for this.

    Attributes:
        fileName (str): The name of the file the chart will be saved as.
        folder (str): The name of the folder in which the chart will be saved.
        chartTitle (str): The title to be printed at the top of the chart.
        fragmentCounts (dict): A dictionary with the names of fragments as its
            keys and the count of each fragment as its values.
    """

    def __init__(self, fileName, folder):
        """Constructs all necessary attributes of the FragmentsPieChart object.

        Args:
            fileName (str): The name of the file the chart will be saved as.
            folder (str): The name of the folder in which the chart will be saved.
        """
        self.fileName = fileName
        self.folder = folder
        self.chartTitle = fileName.capitalize().replace('26_5', '26.5').replace('_', ' ')
        self.fragmentCounts = {}

    def addFragment(self, fragType):
        """Increments the value of a fragment type in the fragmentCounts
        dictionary, if the type is not a key yet it is added and set to one.

        Args:
            fragType (str): Name of the fragment to be added to the dictionary.
        """
        if fragType in self.fragmentCounts:
            self.fragmentCounts[fragType] += 1
        else:
            self.fragmentCounts[fragType] = 1

    def createFragmentsPieChart(self):
        """Creates a pie chart based on the count of each fragment type.
        The chart is titled based on chartTitle and then saved at the
        fileName in the given folder.
        """
        if not self.fragmentCounts:
            return None
        with textStyle():
            plt.figure(figsize=(7, 5))
            plt.pie(self.fragmentCounts.values(),
                    labels=self.fragmentCounts.keys(), textprops={'fontsize': 10})
            plt.title(self.chartTitle, fontsize=14)
            plt.legend(self.fragmentCounts.keys(), loc=(1.15, 0.4), fontsize=10)
            plt.subplots_adjust(left=0.1, right=0.65, top=0.9, bottom=0.1)
            return savePdfAndPng(f'{self.folder}/plots', f'{self.fileName}_fragmentations_pie')
