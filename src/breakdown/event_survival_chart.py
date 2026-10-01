from breakdown.plot_style import (COLORS, MAX_CURVES_WITH_BANDS, fractionYLimits, mixedLegend, newFigure,
                                      plotStyle, saveFigure, styleFor)
from breakdown.run_summary import getSurvivalCurve, loadSection, loadSummary

EVENT_NAMES = {'anyChange': 'Any change', 'isomerization': 'Isomerization',
               'fragmentation': 'Fragmentation', 'recombination': 'Recombination'}
FIRST_EVENTS_HEADING = 'First step of each event type in each run.'


def loadFirstEvents(folder):
    """Get the first event of each type in each run of a folder, and the time per step.

    Args:
        folder (str): Folder with full_log.txt, and run_summary.json if it was made.

    Returns:
        Tuple[dict or None, float or None]: The first events (None if the folder has none, as for
        logs made by older versions), and the time per step in fs stored with them.
    """
    summary = loadSummary(folder)
    timePerStepFs = summary['timePerStepFs'] if summary else None
    return loadSection(folder, FIRST_EVENTS_HEADING), timePerStepFs


# Fixed colours for the event types, the same in every figure and in the run outcome chart
EVENT_COLORS = {'anyChange': '#444444', 'isomerization': COLORS[0], 'fragmentation': COLORS[1],
                'recombination': COLORS[2]}
EVENT_LINE_STYLES = {'anyChange': '-', 'isomerization': '--', 'fragmentation': '-.', 'recombination': ':'}


class EventSurvivalChart:
    """Kaplan-Meier curves: the percentage of runs in which the event has not happened yet, against time.
    A run without the event is counted as still waiting until its last step. With two curves or fewer, the
    95 % Greenwood bands are drawn. The legend gives the median (time at which half of the runs had the event,
    n.r. if fewer than half did within the simulated time) and the number of runs with the event."""

    def __init__(self, chartName, saveFolder):
        """
        Args:
            chartName (str): Name used for the file.
            saveFolder (str): Folder for the plot.
        """
        self.chartName = chartName
        self.saveFolder = saveFolder
        self.series = []

    def addCurve(self, label, firstEvents, event, timePerStepFs, color=None, linestyle=None):
        """Add the survival curve of one event type of one molecule.

        Args:
            label (str): Legend label.
            firstEvents (dict): First events per run.
            event (str): anyChange, isomerization, fragmentation or recombination.
            timePerStepFs (float): Simulated time between steps, in femtoseconds.
            color (str or None): Line colour, the colour of the event type if None.
            linestyle (str or tuple or None): Line style, the style of the event type if None.
        """
        curve = getSurvivalCurve(firstEvents, event, timePerStepFs)
        endTime = max(record['lastStep'] for record in firstEvents.values()) * timePerStepFs / 1000
        self.series.append({'label': label, 'curve': curve, 'endTime': endTime,
                            'color': color or EVENT_COLORS[event], 'linestyle': linestyle or EVENT_LINE_STYLES[event]})

    def createPlot(self, suffix='first_event_survival'):
        """Draw and save the chart.

        Returns:
            str: Path of the pdf.
        """
        drawBands = len(self.series) <= MAX_CURVES_WITH_BANDS
        with plotStyle():
            fig, ax = newFigure('medium', 3.0)
            # Legend in the order of the curves at the end of the plot, so that it can be read along the lines
            ordered = sorted(self.series, key=lambda item: -item['curve']['survival'][-1])
            for item in ordered:
                curve, endTime = item['curve'], item['endTime']
                times = curve['times'] + [endTime]
                median = f"{curve['median']:.1f} ps" if curve['median'] is not None else 'n.r.'
                ax.step(times, [100 * v for v in curve['survival'] + [curve['survival'][-1]]], where='post',
                        color=item['color'], linestyle=item['linestyle'],
                        label=f"{item['label']}: {median} ({curve['numEvents']}/{curve['numRuns']})")
                if drawBands:
                    ax.fill_between(times, [100 * v for v in curve['lower'] + [curve['lower'][-1]]],
                                    [100 * v for v in curve['upper'] + [curve['upper'][-1]]], step='post',
                                    alpha=0.15, color=item['color'], linewidth=0)
            bottom = fractionYLimits(ax, margin=6)
            if bottom == 0:
                ax.axhline(50, color='#777777', linewidth=0.6, linestyle=(0, (2, 2)))
            ax.set_xlabel('Time (ps)')
            ax.set_ylabel('Runs without the event (%)')
            ax.set_xlim(left=0)
            legend = mixedLegend(ax, None, drawBands, loc='center left', bbox_to_anchor=(1.01, 0.5),
                                 title='Median (events/runs)', title_fontsize=7.5)
            return saveFigure(fig, self.saveFolder, f'{self.chartName}_{suffix}')


def generateEventComparison(folders, names, event, title, outputDir, timePerStepFs=None):
    """Compare the time to the first event of one type between molecules.

    Args:
        folders (list): Folders with full_log.txt.
        names (list): Legend labels of the molecules.
        event (str): anyChange, isomerization, fragmentation or recombination.
        title (str): Name of the chart and the file.
        outputDir (str): Folder for the pdf.
        timePerStepFs (float or None): Time between steps in fs, taken from the summaries if None.

    Returns:
        str: Path of the pdf.

    Raises:
        ValueError: Raised if a folder has no first event data, run makechartsandstats again for it.
    """
    chart = EventSurvivalChart(title, outputDir)
    styles = styleFor(names)
    for folder, name in zip(folders, names):
        firstEvents, storedTime = loadFirstEvents(folder)
        if not firstEvents:
            raise ValueError(f'{folder} has no first event data, run makechartsandstats for it again.')
        chart.addCurve(name, firstEvents, event, timePerStepFs or storedTime or 5.0, *styles[name])
    return chart.createPlot(suffix=f'{event}_survival')
