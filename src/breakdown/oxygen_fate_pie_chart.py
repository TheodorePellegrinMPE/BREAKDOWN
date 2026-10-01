import matplotlib.pyplot as plt
from breakdown.plot_style import savePdfAndPng, textStyle
from breakdown.run_summary import loadSection

class OxygenFatePieChart:
    def __init__(self, chartName, saveFolder):
        self.chartName = chartName
        self.saveFolder = saveFolder
        self.dataSeries = []

    def addDataFromFolders(self, folderPaths, abbreviations, groupNames):
        for folder, abbr, group in zip(folderPaths, abbreviations, groupNames):
            endingStats = loadSection(folder,
                                      'Runs with each oxygen in each final group and/or fragment.',
                                      'Runs with each sulfur in each final group and/or fragment.')
            if endingStats:
                self.processFateData(abbr, group, endingStats)
            else:
                print(f'No data on final groups in {folder}, was makechartsandstats run with -f {group}?')

    def processFateData(self, abbr, startGroup, endingStats):
        startGroupClean = startGroup.lower()
        fates = {}
        totalDissociated = 0

        for key, count in endingStats.items():
            if key.lower() == startGroupClean:
                continue
            
            fates[key] = count
            totalDissociated += count
        
        if totalRuns := sum(endingStats.values()):
            self.dataSeries.append({
                'abbr': abbr,
                'group': startGroup,
                'totalRuns': totalRuns,
                'totalDissociated': totalDissociated,
                'fates': fates
            })

    def createPiePlots(self):
        n = len(self.dataSeries)
        if n == 0:
            return None
        with textStyle(), plt.rc_context({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'}):
            return self.drawPiePlots(n)

    def drawPiePlots(self, n):
        """Draws one pie chart for every molecule, see createPiePlots."""

        all_fates = sorted(list(set(label for d in self.dataSeries for label in d['fates'].keys())))
        cmap = plt.cm.Set3
        color_map = {label: cmap(i / (len(all_fates) if len(all_fates) > 0 else 1)) for i, label in enumerate(all_fates)}

        cols = 2
        rows = (n + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(4.5 * cols, 4.5 * rows), dpi=300, squeeze=False)
        axes = axes.flatten()

        for i, data in enumerate(self.dataSeries):
            ax = axes[i]
            if not data['fates']:
                ax.text(0.5, 0.5, "100% Parent Intact", ha='center', fontsize=12)
                ax.set_title(f"{data['abbr']} ({data['group']})", fontweight='bold')
                ax.axis('off')
                continue

            labels = list(data['fates'].keys())
            values = list(data['fates'].values())
            current_colors = [color_map[l] for l in labels]
            
            wedges, texts, autotexts = ax.pie(
                values, 
                labels=labels, 
                autopct=lambda p: '{:.0f}'.format(p * sum(values) / 100),
                startangle=140, 
                colors=current_colors,
                pctdistance=0.75,
                wedgeprops={'linewidth': 1.2, 'edgecolor': 'white'}
            )

            plt.setp(autotexts, size=8, weight="bold")
            plt.setp(texts, size=9)
            
            ax.set_title(f"{data['abbr']} ({data['group']})\n($N_{{changed}}$ = {data['totalDissociated']}/{data['totalRuns']})", 
                         fontsize=11, fontweight='bold', pad=15)

        for j in range(i + 1, len(axes)):
            axes[j].axis('off')

        plt.tight_layout()
        return self.savePlotToDisk()

    def savePlotToDisk(self):
        return savePdfAndPng(self.saveFolder, f"{self.chartName}_oxygen_fate_pies")

def generateOxygenFates(folders, names, groups, title, outputDir):
    chart = OxygenFatePieChart(title, outputDir)
    chart.addDataFromFolders(folders, names, groups)
    return chart.createPiePlots()