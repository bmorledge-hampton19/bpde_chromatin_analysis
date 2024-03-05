# This script uses a specialized peak calling algorithm to call dyad positions from a table of counts.
# (E.g., for determining nucleosome positions around CTCF binding sites).
from pandas import DataFrame


def callDyadPeaks(dyadCounts: DataFrame, countsColumn = "Both_Strands_Counts",
                  smoothingWindow = None, blacklistedRegions = None,
                  peakExclusionRadius = 146, boundaryExclusion = None, maxPeaks = None, meanRatioCutoff = 1):
    """
    Given a pandas data frame of dyad counts, returns a list of dyad positions that can be used to define concrete nucleosome positions.
    - countsColumn defines the column in the data table to use.
    - smoothingWindow should be None for no smoothing or an integer value for the rolling mean window
    - blacklistedRegion should be None or some iterable to define regions to remove prior to peak calling.
    - peakExclusionRadius describes how many positions on either side of each peak are invalidated for further peak calling.
    - boundaryExclusion describes the number of positions at each end that will be used for smoothing but will not be considered reliable for peak calling.
      If these positions are called as peaks, they still exclude positions around them according to the peakExclusion radius but will not be added to the list of peaks.
      If set to None, the value defaults to half the smoothing window, rounded up (or 1 if the smoothing window is None as well).
    - maxPeaks and meanRatioCutoff are termination conditions. If maxPeaks peaks are found, or
      if none of the remaining values are greater than the original mean (after blacklisting) times meanRatioCutoff, peak calling is finished.
    """

    # Start by copying the given data frame, so as not to alter the original.
    dyadCounts = dyadCounts.copy()

    # If requested, set all blacklisted regions to None. This ensures that they do not impact smoothing.
    if blacklistedRegions is not None:
        dyadCounts.loc[dyadCounts["Dyad_Position"].isin(blacklistedRegions), countsColumn] = None

    # Compute the cutoff value.
    cutoff = dyadCounts[countsColumn].mean()*meanRatioCutoff

    # If requested, smooth the data.
    if smoothingWindow is not None:
        dyadCounts[countsColumn+"_Smoothed"] = dyadCounts[countsColumn].rolling(11, center = True, min_periods=1).mean()
        countsColumn = countsColumn+"_Smoothed"

    # Determine the range for boundary exclusion if it is still None.
    if boundaryExclusion is None:
        if smoothingWindow is None: boundaryExclusion = 1
        else: boundaryExclusion = round(smoothingWindow/2)

    # Begin the core peak calling algorithm. Maximum positions are called as peaks and
    # surrounding values are zeroed out based on the peak exclusion radius until either termination condition is reached.
    peaks = list()
    maxValue = dyadCounts[countsColumn].max()
    maxPosition = dyadCounts.loc[dyadCounts[countsColumn]==maxValue, "Dyad_Position"].iloc[0]
    while (maxPeaks is None or len(peaks) < maxPeaks) and maxValue > cutoff:

        if abs(maxPosition) <= dyadCounts["Dyad_Position"].max() - boundaryExclusion: peaks.append(maxPosition)

        dyadCounts.loc[dyadCounts["Dyad_Position"].isin(range(maxPosition-peakExclusionRadius, maxPosition+peakExclusionRadius+1)),
                       countsColumn] = 0

        maxValue = dyadCounts[countsColumn].max()
        maxPosition = dyadCounts.loc[dyadCounts[countsColumn]==maxValue, "Dyad_Position"].iloc[0]

    return sorted(peaks)