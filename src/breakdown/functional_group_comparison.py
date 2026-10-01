import matplotlib.pyplot as plt
from breakdown.plot_style import savePdfAndPng, textStyle
from breakdown.run_summary import loadSection, resolveTimePerStep

class FunctionalGroupComparisonChart:
    def __init__(self, chartName, saveFolder, timePerStepFs=5.0):
        """Initializes the FunctionalGroupComparisonChart object.

        Args:
            chartName (str): The name to be used for the output file and title.
            saveFolder (str): The directory where the resulting plot will be saved.
            timePerStepFs (float): Simulated time between steps of the logs, in femtoseconds.
        """
        self.chartName = chartName
        self.saveFolder = saveFolder
        self.timeStepToPs = timePerStepFs / 1000
        # T50 is the time at which half of the runs have lost the group
        self.thresholdFraction = 0.5
        self.dataSeries = []

    def addDataFromFolders(self, folderPaths, abbreviations, groupNames):
        """Loads data from multiple simulation directories.

        Args:
            folderPaths (list): Paths to the directories containing full_log.txt.
            abbreviations (list): Display names for the molecules.
            groupNames (list): The specific functional groups to track for each.
        """
        for folder, abbr, group in zip(folderPaths, abbreviations, groupNames):
            countDict = loadSection(folder, f'Runs containing original functional group {group} at step.')
            if countDict:
                self.processSeries(abbr, group, countDict)
            else:
                print(f'No data for group {group} in {folder}, was makechartsandstats run with -f {group}?')

    def processSeries(self, abbreviation, groupName, countDict):
        """Calculates statistics and stores the data for plotting.

        Args:
            abbreviation (str): Molecule abbreviation.
            groupName (str): Functional group name.
            countDict (dict): The raw step-to-count dictionary.
        """
        timeSeries = self.convertStepsToTime(countDict)
        threshold = self.thresholdFraction * max(count for _, count in timeSeries)
        t50 = self.calculateT50(timeSeries, threshold)
        self.dataSeries.append({
            'threshold': threshold,
            'abbr': abbreviation,
            'group': groupName,
            'points': timeSeries,
            't50': t50
        })

    def convertStepsToTime(self, countDict):
        """Converts step strings to integer steps and scales to picoseconds.

        Args:
            countDict (dict): Dictionary with string keys (steps).

        Returns:
            list: Sorted list of tuples (time_ps, count).
        """
        sortedKeys = sorted([int(k) for k in countDict.keys()])
        return [(step * self.timeStepToPs, countDict[str(step)]) for step in sortedKeys]

    def calculateT50(self, timeSeries, threshold):
        """Identifies the first time value where the count is half the starting count.

        Args:
            timeSeries (list): List of (time, count) tuples.
            threshold (float): Count at or under which half of the runs have lost the group.

        Returns:
            float or None: The time in ps, or None if it never drops that far.
        """
        for time, count in timeSeries:
            if count <= threshold:
                return time
        return None

    def createComparisonPlot(self):
        """Generates the final comparison plot with professional styling.

        Returns:
            str or None: Path of the pdf, None if there was no data.
        """
        if not self.dataSeries:
            return None
        with textStyle(), plt.rc_context({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'}):
            fig, ax = plt.subplots(figsize=(8, 5), dpi=300)

            colors = plt.cm.Set1(range(9))
            lineStyles = ['-', '--', ':', '-.', (0, (3, 5, 1, 5)), (0, (5, 10))]

            for i, data in enumerate(self.dataSeries):
                self.plotSingleLine(ax, data, colors[i % 9], lineStyles[i % len(lineStyles)])

            self.applyPlotFormatting(ax)
            return self.savePlotToDisk()

    def plotSingleLine(self, ax, data, color, lineStyle):
        """Plots a single data series onto the provided axis.

        Args:
            ax (matplotlib.axes.Axes): The axis to plot on.
            data (dict): Dictionary containing points and metadata.
            color (tuple): RGBA color tuple.
            lineStyle (str/tuple): Matplotlib linestyle.
        """
        times = [p[0] for p in data['points']]
        counts = [p[1] for p in data['points']]
        label = self.generateLegendString(data)
        ax.plot(times, counts, label=label, color=color,
                linestyle=lineStyle, linewidth=2, alpha=0.9)

    def generateLegendString(self, data):
        """Creates a legend label with an approximate T50 if the count is dropping.

        Args:
            data (dict): Series metadata including points and calculated t50.

        Returns:
            str: Formatted legend string with LaTeX subscripts and approximate T50.
        """
        t50Val = data['t50']
        points = data['points']
        lastTime = points[-1][0]

        if t50Val is not None:
            t50Str = f"{t50Val:.1f}"
        else:
            if len(points) > 1:
                firstCount = points[0][1]
                lastCount = points[-1][1]

                if lastCount < firstCount and lastCount > data['threshold']:
                    slope = (lastCount - firstCount) / lastTime
                    approxT50 = (data['threshold'] - firstCount) / slope
                    t50Str = f"~{approxT50:.0f}"
                else:
                    t50Str = f">{lastTime:.0f}"
            else:
                t50Str = "n/a"

        return f"{data['abbr']} ({data['group']}): $T_{{50}}$={t50Str} ps"

    def applyPlotFormatting(self, ax):
        """Sets the aesthetic parameters for publication standards.

        Args:
            ax (matplotlib.axes.Axes): The axis to format.
        """
        ax.set_xlabel('Time (ps)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Active Simulations ($N$)', fontsize=12, fontweight='bold')

        ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.3)
        ax.tick_params(axis='both', which='major', labelsize=10, width=1.5, length=6)

        legend = ax.legend(fontsize=9, loc='best', frameon=True, edgecolor='black')
        legend.get_frame().set_linewidth(0.8)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)

        plt.tight_layout()

    def savePlotToDisk(self):
        """Saves the plot as a pdf and a png, returns the path of the pdf."""
        return savePdfAndPng(self.saveFolder, f"{self.chartName}_comparison")

def generateComparison(folders, names, groups, title, outputDir, timePerStepFs=None):
    """Convenience function to generate a comparison chart.

    Args:
        folders (list): Directories to analyze.
        names (list): Abbreviations for the legend.
        groups (list): Functional groups to extract.
        title (str): Output filename.
        outputDir (str): Output directory.
        timePerStepFs (float or None): Simulated time between steps of the logs, in femtoseconds,
            read from the summaries if None.

    Returns:
        str or None: Path of the pdf.
    """
    chart = FunctionalGroupComparisonChart(title, outputDir, resolveTimePerStep(folders, timePerStepFs))
    chart.addDataFromFolders(folders, names, groups)
    return chart.createComparisonPlot()
