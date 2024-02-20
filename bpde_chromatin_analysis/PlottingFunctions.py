# A suite of plotting functions for the BPDE analysis.
import pandas
from plotnine import *
from benbiohelpers.Plotting.PlotnineHelpers import defaultTextScaling, blankBackground


def plotTFBSData(damageData: pandas.DataFrame, repairData: pandas.DataFrame, dataColumn = "Normalized_Both_Strands",
                 title = "TFBS Damage vs. Repair", xAxisLabel = "Position Relative to Midpoint (bp)", yAxisLabel = "Normalized Damage/Repair Activity",
                 xlim = None, ylim = None, dropZeroRows = True, damageColor = "#707070", repairColor = "#F1A226"):
    """
    Plots repair and damage data relative to transcription factor binding site midpoints.
    """

    assert damageData is not None or repairData is not None

    # Get rid of any rows with no counts.
    if dropZeroRows:
        if damageData is not None:
            damageData = damageData.loc[damageData[dataColumn] != 0].copy()
        if repairData is not None:
            repairData = repairData.loc[repairData[dataColumn] != 0].copy()

    # Generate the plot.
    TFBSPlot = ggplot()

    if damageData is not None:
        TFBSPlot = TFBSPlot + geom_path(aes(x = "Dyad_Position", y = dataColumn, group = 1, color = 'damageColor'),
                                        data = damageData, size = 1, lineend = "round")
    if repairData is not None:
        TFBSPlot = TFBSPlot + geom_path(aes(x = "Dyad_Position", y = dataColumn, group = 1, color = 'repairColor'),
                                        data = repairData, size = 1, lineend = "round")

    if repairData is not None and damageData is not None:
        TFBSPlot = (TFBSPlot +
            scale_color_identity(guide = "legend", name = " ", breaks = (damageColor,repairColor), labels = ("damage","repair")) +
            theme(figure_size = (12,6))
        )
    else: TFBSPlot = TFBSPlot + theme(figure_size=(9,6))

    TFBSPlot = (TFBSPlot + 
        coord_cartesian(xlim = xlim, ylim = ylim) +
        xlab(xAxisLabel) + ylab(yAxisLabel) + ggtitle(title) +
        defaultTextScaling + blankBackground
    )

    return TFBSPlot


def plotTSSData(damageData: pandas.DataFrame, repairData: pandas.DataFrame, dataColumn = "Normalized_Both_Strands",
                title = "TSS Damage vs. Repair", xAxisLabel = "Position Relative to TSS (bp)", yAxisLabel = "Normalized Damage/Repair Activity",
                xlim = None, ylim = None, damageColor = "#707070", repairColor = "#F1A226"):
    """
    Turns out, this works pretty much the same as the TFBS data!
    """
    return plotTFBSData(damageData, repairData, dataColumn, title, xAxisLabel, yAxisLabel, xlim, ylim, False, damageColor, repairColor)


def plotMiscData(data: pandas.DataFrame, dataColumn = "Normalized_Both_Strands",
                 title = "Misc Data", xAxisLabel = "Position Relative to Feature (bp)", yAxisLabel = "Counts",
                 xlim = None, ylim = None, color = "black"):
    """
    Yup. This one just uses the TFBS plotting function again.
    """
    return plotTFBSData(data, None, dataColumn, title, xAxisLabel, yAxisLabel, xlim, ylim, False, color)