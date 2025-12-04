# A suite of plotting functions for the BPDE analysis.
import pandas
from enum import Enum
from typing import List, Dict
from plotnine import *
from benbiohelpers.Plotting.PlotnineHelpers import defaultTextScaling, blankBackground, flushXAxis
from mutperiodpy.GeneratePlotnineFigures import NucleosomeColorPalette


class Plottable:
    """
    Wraps together the necessary information for plotting a data set, so it easier to pass to the generalPlotting function.
    """

    def __init__(self, dataFrame: pandas.DataFrame, dataColumn: str, label,
                 color, underlaidColor = None, colorColumn = None, underlaidColorColumn = None,
                 colorColumnBreaksAndLabels: Dict[str,str] = None,
                 smoothData = False, overlaySmoothedAndNormal = False):
        self.dataFrame = dataFrame.copy()
        self.dataColumn = dataColumn
        self.smoothedDataColumn = dataColumn + "_Smoothed"
        self.label = label
        self.color = color
        self.underlaidColor = underlaidColor
        self.colorColumn = colorColumn
        self.underlaidColorColumn = underlaidColorColumn
        self.columnColorBreaksAndLabels = colorColumnBreaksAndLabels
        self.smoothData = smoothData
        self.overlaySmoothedAndNormal = overlaySmoothedAndNormal
        

def generalPlotting(plottables: List[Plottable],
                     title = " ", xAxisLabel = "Position Relative to Feature (bp)", yAxisLabel = "Normalized Activity",
                     xlim = None, ylim = None, dropZeroRows = True, smoothData = True, overlaySmoothedAndNormal = True):
    """
    Generalized plotting function for plotting an arbitrary number of features together (e.g., damage and repair).
    """

    # Get rid of any rows with no counts.
    if dropZeroRows:
        for plottable in plottables:
            plottable = plottable.dataFrame.loc[plottable.dataFrame[plottable.dataColumn] != 0].copy()

    # If requested, smooth the data.
    if smoothData:
        for plottable in plottables:
            plottable.dataFrame[plottable.smoothedDataColumn] = plottable.dataFrame[plottable.dataColumn].rolling(11, center = True, min_periods=1).mean()

    # Generate the plot.
    plot = ggplot(mapping = aes(x = "Dyad_Position"))

    if overlaySmoothedAndNormal:
        for plottable in plottables:
            if plottable.underlaidColorColumn is not None: underlaidColorString = plottable.underlaidColorColumn
            else: underlaidColorString = f"'{plottable.underlaidColor}'"

            plot = plot + geom_path(aes(y = plottable.dataColumn, color = f"{underlaidColorString}"), plottable.dataFrame, size = 1, group = 1, lineend = "round")

    for plottable in plottables:
        if plottable.colorColumn is not None: colorString = plottable.colorColumn
        else: colorString = f"'{plottable.color}'"

        if smoothData:
            plot = plot + geom_path(aes(y = plottable.smoothedDataColumn, color = f"{colorString}"), plottable.dataFrame, size = 1, group = 1, lineend = "round")
        else:
            plot = plot + geom_path(aes(y = plottable.dataColumn, color = f"{colorString}"), plottable.dataFrame, size = 1, group = 1, lineend = "round")

    if len(plottables) > 1:
        breaks = list()
        labels = list()
        for plottable in plottables:
            if plottable.columnColorBreaksAndLabels is not None:
                breaks += plottable.columnColorBreaksAndLabels.keys()
                labels += plottable.columnColorBreaksAndLabels.values()
            else:
                breaks.append(plottable.color)
                labels.append(plottable.label)
        plot = (plot +
            scale_color_identity(guide = "legend", name = " ", breaks = breaks, labels = labels) +
            theme(figure_size = (12,6))
        )
    else: plot = plot + scale_color_identity(guide = None, name = " ") + theme(figure_size=(9,6))

    plot = (plot + 
        coord_cartesian(xlim = xlim, ylim = ylim) +
        xlab(xAxisLabel) + ylab(yAxisLabel) + ggtitle(title) +
        defaultTextScaling + blankBackground
    )

    return plot


class PlottingColors(Enum):

    BLACK = "#000000"
    BLACK_UNDERLAY = "#BFBFBF"
    DAMAGE = "#707070"
    DAMAGE_UNDERLAY = "#CCCCCC"
    REPAIR = "#F29500"
    REPAIR_UNDERLAY = "#F2DBB6"
    MUTATION = "#CC0000"
    MUTATION_UNDERLAY = "#CC9999"


def plotDamageAndRepair(damageData: pandas.DataFrame, repairData: pandas.DataFrame, dataColumn = "Normalized_Both_Strands",
                        title = "Damage vs. Repair", xAxisLabel = "Position Relative to Feature (bp)", yAxisLabel = "Normalized Damage/Repair Activity",
                        xlim = None, ylim = None, dropZeroRows = True, smoothData = True, overlaySmoothedAndNormal = True,
                        damageColor = PlottingColors.DAMAGE.value, repairColor = PlottingColors.REPAIR.value,
                        underlaidDamageColor = PlottingColors.DAMAGE_UNDERLAY.value, underlaidRepairColor = PlottingColors.REPAIR_UNDERLAY.value):
    """
    Plot damage and repair together using the generalizable plotting function
    """
    damagePlottable = Plottable(damageData, dataColumn, "Damage", damageColor, underlaidDamageColor)
    repairPlottable = Plottable(repairData, dataColumn, "Repair", repairColor, underlaidRepairColor)
    return generalPlotting([damagePlottable, repairPlottable],
                           title, xAxisLabel, yAxisLabel, xlim, ylim,
                           dropZeroRows, smoothData, overlaySmoothedAndNormal)


def plotDamageRepairAndMutagenesis(
        damageData: pandas.DataFrame, repairData: pandas.DataFrame, mutationData: pandas.DataFrame, dataColumn = "Normalized_Both_Strands",
        title = "Damage vs. Repair", xAxisLabel = "Position Relative to Feature (bp)", yAxisLabel = "Enrichment",
        xlim = None, ylim = None, dropZeroRows = True, smoothData = True, overlaySmoothedAndNormal = True,
        damageColor = PlottingColors.DAMAGE.value, underlaidDamageColor = PlottingColors.DAMAGE_UNDERLAY.value,
        repairColor = PlottingColors.REPAIR.value, underlaidRepairColor = PlottingColors.REPAIR_UNDERLAY.value,
        mutationColor = PlottingColors.MUTATION.value, underlaidMutationColor = PlottingColors.MUTATION_UNDERLAY.value
    ):
    """
    Plot damage, repair, and mutations together using the generalizable plotting function
    """
    damagePlottable = Plottable(damageData, dataColumn, "Damage", damageColor, underlaidDamageColor)
    repairPlottable = Plottable(repairData, dataColumn, "Repair", repairColor, underlaidRepairColor)
    mutationPlottable = Plottable(mutationData, dataColumn, "Mutations", mutationColor, underlaidMutationColor)
    return generalPlotting([damagePlottable, repairPlottable, mutationPlottable],
                           title, xAxisLabel, yAxisLabel, xlim, ylim,
                           dropZeroRows, smoothData, overlaySmoothedAndNormal)


def plotSingleFeature(data: pandas.DataFrame, dataColumn = "Normalized_Both_Strands",
                      title = "Misc Data", xAxisLabel = "Position Relative to Feature (bp)", yAxisLabel = "Counts",
                      xlim = None, ylim = None, dropZeroRows = True, smoothData = True, overlaySmoothedAndNormal = True,
                      color = PlottingColors.BLACK.value, underlaidColor = PlottingColors.BLACK_UNDERLAY.value):
    """
    Plot a single arbitrary feature using the generalizable plotting function.
    """
    return generalPlotting([Plottable(data, dataColumn, " ", color, underlaidColor)],
                           title, xAxisLabel, yAxisLabel, xlim, ylim,
                           dropZeroRows, smoothData, overlaySmoothedAndNormal)


def plotFeatureWithCustomNucleosomes(data: pandas.DataFrame, nucleosomePositions: List[int], dataColumn = "Normalized_Both_Strands",
                                     title = "Misc Data", xAxisLabel = "Position Relative to Feature (bp)", yAxisLabel = "Counts",
                                     xlim = None, ylim = None, dropZeroRows = True, smoothData = True, overlaySmoothedAndNormal = True,
                                     nucleosomeColorPalette = NucleosomeColorPalette(), blacklistedRegions = None):
    """
    Plots a single arbitrary feature with coloring for given nucleosome positions.
    Uses the generalizable plotting function.
    """
    data = data.copy()
    data["Color"] = nucleosomeColorPalette.linker
    data["Underlaid_Color"] = nucleosomeColorPalette.linkerUnderlaid
    if blacklistedRegions is not None:
        data.loc[data["Dyad_Position"].isin(blacklistedRegions),"Color"] = PlottingColors.BLACK.value
        data.loc[data["Dyad_Position"].isin(blacklistedRegions),"Underlaid_Color"] = PlottingColors.BLACK_UNDERLAY.value
    for nucleosomePosition in nucleosomePositions:
        data.loc[data["Dyad_Position"].isin(range(nucleosomePosition-73,nucleosomePosition+73+1)),"Color"] = nucleosomeColorPalette.nucleosomal
        data.loc[data["Dyad_Position"].isin(range(nucleosomePosition-73,nucleosomePosition+73+1)),"Underlaid_Color"] = nucleosomeColorPalette.nucleosomalUnderlaid

    return generalPlotting([Plottable(data, dataColumn, " ", None, None, "Color", "Underlaid_Color",
                                      {nucleosomeColorPalette.nucleosomal:"Nucleosomal", nucleosomeColorPalette.linker:"Linker"})],
                           title, xAxisLabel, yAxisLabel, xlim, ylim,
                           dropZeroRows, smoothData, overlaySmoothedAndNormal)


def plotReadLengthDistribution(data: pandas.DataFrame, readLengthDataCol = "Sequence_Length", readCountsDataCol = "Read_Count",
                               title = "Read Length Distribution", xAxisLabel = "Read Length", yAxisLabel = "Number of Reads",
                               xlim = None, xAxisBreaks = True, ylim = None, yAxisBreaks = True) -> ggplot:
    """
    Plots read length distribution using an enrinched indices data frame.
    """
    readLengthPlot = (
        ggplot(data, aes(x = readLengthDataCol, y = readCountsDataCol, fill = '"black"')) +
        geom_bar(stat = "identity") + scale_fill_identity() +
        coord_cartesian(xlim = xlim, ylim = ylim) +
        scale_x_continuous(breaks = xAxisBreaks) + scale_y_continuous(breaks = yAxisBreaks, expand = flushXAxis) + 
        xlab(xAxisLabel) + ylab(yAxisLabel) + ggtitle(title) +
        theme(figure_size=(9,6)) + defaultTextScaling + blankBackground
    )

    return readLengthPlot


def plotEnrichedEndices(data: pandas.DataFrame, readLengthDataCol = "Sequence_Length", maxFrequencyDataCol = "G_Max_Frequency",
                        maxFrequenceyPosDataCol = "G_Max_Frequency_Position",
                        title = "Enriched indices", xAxisLabel = "Read Length", yAxisLabel = "Maximum Guanine Frequency",
                        xlim = None, xAxisBreaks = True, ylim = None, yAxisBreaks = True) -> ggplot:
    """
    Plots relative enrichment of sequences at each read length using an enriched indices data frame.
    """

    enrichedIndicesPlot = (
        ggplot(data, aes(x = readLengthDataCol, y = maxFrequencyDataCol, fill = '"black"')) +
        geom_bar(stat = "identity") + scale_fill_identity() +
        geom_text(aes(label = maxFrequenceyPosDataCol), nudge_y = 0.02, format_string = "{:.0f}") +
        coord_cartesian(xlim = xlim, ylim = ylim) +
        scale_x_continuous(breaks = xAxisBreaks) + scale_y_continuous(breaks = yAxisBreaks, expand = flushXAxis) + 
        xlab(xAxisLabel) + ylab(yAxisLabel) + ggtitle(title) +
        theme(figure_size=(9,6)) + defaultTextScaling + blankBackground
    )

    return enrichedIndicesPlot


def plotIndividualFrequencies(data: pandas.DataFrame, readLengthDataCol = "Sequence_Length", positionDataCol = "Position",
                              nucFreqDataCol = "G_Frequency", validReadLengths = range(22,30),
                              title = "Guanine Frequency By Read Length", xAxisLabel = "Position Relative to 3' End", yAxisLabel = "Guanine Frequency",
                              xlim = None, ylim = None) -> ggplot:
    """
    Plots Individual nucleotide frequencies across a given set of read lengths.
    """

    data = data.loc[(readLength in validReadLengths for readLength in data[readLengthDataCol])]

    individualFrequenciesPlot = (
        ggplot(data, aes(x = positionDataCol, y = nucFreqDataCol)) +
        geom_bar(aes(fill = '"black"'), stat = "identity") + scale_fill_identity() +
        facet_grid(readLengthDataCol) +
        coord_cartesian(xlim = xlim, ylim = ylim) + scale_y_continuous(expand = flushXAxis) +
        xlab(xAxisLabel) + ylab(yAxisLabel) + ggtitle(title) +
        theme(figure_size=(9,1+1.5*len(validReadLengths))) + defaultTextScaling + blankBackground
    )

    return individualFrequenciesPlot