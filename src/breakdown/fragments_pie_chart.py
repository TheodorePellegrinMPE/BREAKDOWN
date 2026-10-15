import re
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

    def addRun(self, run):
        """Count the fragment that a run ended with.

        If there is only 1 formula at the end, then the molecule has not
        fragmented. If there are 3 or more than there are multiple fragments.
        If there are 2, the second, smaller fragment is
        recorded, fragments of C and H or C, H and O are stacked, this may be
        changed in the future.

        Args:
            run (Run): The run, see run_log.
        """
        fragments = run.final
        if len(fragments) < 2:
            self.addFragment('No fragmentation')
        elif len(fragments) > 2:
            self.addFragment('Multiple fragments')
        else:
            fragment = fragments[1]
            if re.match(r'^C\d*H\d*O\d*$', fragment):
                self.addFragment('CmHnOk')
            elif re.match(r'^C\d*H\d*S\d*$', fragment):
                self.addFragment('CmHnSk')
            elif re.match(r'^C\d*H\d*$', fragment) and (fragment != 'C2H2'):
                self.addFragment('CmHn')
            else:
                self.addFragment(fragment)

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
