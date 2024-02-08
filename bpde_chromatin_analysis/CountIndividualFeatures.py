# This script preforms the usual analysis of counting mutations, repair, damage, etc. around features like nucleosomes but
# with additional stratification by each feature, in order to identify those with unusual counts.
import os
from typing import List
from helper_scripts.BPDE_DataDir import getDataDir
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.InputParsing.CheckForNumber import checkForNonNegativeInteger
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.InputDataStructures import EncompassingDataDefaultStrand, EncompassingData, ENCOMPASSING_DATA
from benbiohelpers.CountThisInThat.CounterOutputDataHandler import CounterOutputDataHandler, AmbiguityHandling, OutputDataWriter


# Define the necessary functions/classes for counting.
def getCountDerivatives(outputDataWriter: OutputDataWriter, getHeaders):
    if getHeaders: return ["Both_Strands_Counts", "Aligned_Strands_Counts"]
    else:
        thisPlusCounts = outputDataWriter.outputDataStructure[outputDataWriter.previousKeys[0]][outputDataWriter.previousKeys[1]][True]
        thisMinusCounts = outputDataWriter.outputDataStructure[outputDataWriter.previousKeys[0]][outputDataWriter.previousKeys[1]][False]
        oppositeMinusCounts = outputDataWriter.outputDataStructure[outputDataWriter.previousKeys[0]][-outputDataWriter.previousKeys[1]][False]
        return [str(thisPlusCounts+thisMinusCounts),str(thisPlusCounts+oppositeMinusCounts)]


class MutationsInNucleosomesCounter(ThisInThatCounter):

    def setUpOutputDataHandler(self):
        self.outputDataHandler = CounterOutputDataHandler(self.writeIncrementally, trackAllEncompassing=True)
        self.outputDataHandler.addEncompassingFeatureStratifier()
        self.outputDataHandler.addRelativePositionStratifier(self.currentEncompassingFeature, extraRangeRadius = self.encompassingFeatureExtraRadius,
                                                             outputName = "Dyad_Position")
        self.outputDataHandler.addStrandComparisonStratifier(strandAmbiguityHandling = AmbiguityHandling.tolerate)
        self.outputDataHandler.createOutputDataWriter(self.outputFilePath, customStratifyingNames=(None, None, {True:"Plus_Strand_Counts", False:"Minus_Strand_Counts"}),
                                                      getCountDerivatives = getCountDerivatives, writeHeadersImmediately = True)

    def constructEncompassingFeature(self, line) -> EncompassingDataDefaultStrand:
        return EncompassingDataDefaultStrand(line, self.acceptableChromosomes)


class MutationsInStrandedNucleosomesCounter(ThisInThatCounter):

    def setUpOutputDataHandler(self):
        self.outputDataHandler = CounterOutputDataHandler(self.writeIncrementally, trackAllEncompassing=True)
        self.outputDataHandler.addEncompassingFeatureStratifier()
        self.outputDataHandler.addRelativePositionStratifier(self.currentEncompassingFeature, extraRangeRadius = self.encompassingFeatureExtraRadius,
                                                             outputName = "Dyad_Position", strandSpecificPos = True)
        self.outputDataHandler.addStrandComparisonStratifier(strandAmbiguityHandling = AmbiguityHandling.tolerate)
        self.outputDataHandler.createOutputDataWriter(self.outputFilePath, customStratifyingNames=(None, None, {True:"Plus_Strand_Counts", False:"Minus_Strand_Counts"}),
                                                      getCountDerivatives = getCountDerivatives, writeHeadersImmediately = True)

    def constructEncompassingFeature(self, line) -> EncompassingData:
        return EncompassingData(line, self.acceptableChromosomes)
    

def countIndividualFeatures(encompassedFeaturesFilePaths: List[str], encompassingFeaturesFilePaths: List[str],
                            encompassingFeaturesAreStranded = False, encompassingFeatureExtraRadius = 0):
    """
    Modifies the normal nucleosome counting classes used in mutperiod to also stratify by individual encompassing features.
    Two lists are passed the function, one for encompassed features and one for encompassing features. Each encompassed feature is paired with each encompassing feature.
    If the encompassing features are stranded (e.g., transcription start sites), the encompassingFeaturesAreStranded variable should be set to True.
    Returns a list of the output file paths, which are stored relative to the encompassed features file paths.
    """
    
    # Set the relevant counter based on whether or not the encompassing features are stranded.
    if encompassingFeaturesAreStranded: Counter = MutationsInStrandedNucleosomesCounter
    else: Counter = MutationsInNucleosomesCounter

    outputFilePaths: List[str] = list()

    # Count in each unique pair of encompassed and encompassing features.
    for encompassedFeaturesFilePath in encompassedFeaturesFilePaths:
        for encompassingFeaturesFilePath in encompassingFeaturesFilePaths:

            print(f"\nWorking with {os.path.basename(encompassedFeaturesFilePath)} and {os.path.basename(encompassingFeaturesFilePath)}...")

            outputFilePath = os.path.join(os.path.dirname(encompassedFeaturesFilePath),
                                          os.path.basename(encompassingFeaturesFilePath).rsplit('.',1)[0]+'_'+
                                          os.path.basename(encompassedFeaturesFilePath).rsplit('.',1)[0]+"_individual_feature_counts.tsv")
            outputFilePaths.append(outputFilePath)

            Counter(encompassedFeaturesFilePath, encompassingFeaturesFilePath, outputFilePath,
                    encompassingFeatureExtraRadius = encompassingFeatureExtraRadius, writeIncrementally = ENCOMPASSING_DATA).count()

    return outputFilePaths


def main():
    
    with TkinterDialog(workingDirectory = getDataDir(), title = "Analyze Individual Features") as dialog:
        dialog.createMultipleFileSelector("Encompassed Features:", 0, "context_mutations.bed", ("Bed files", ".bed"))
        dialog.createMultipleFileSelector("Encompassing Features:", 1, "nuc_map.bed", ("Bed files", ".bed"))
        dialog.createCheckbox("Encompassing features are stranded", 2, 0)
        dialog.createTextField("Encompassing feature extra radius", 3, 0, defaultText="100")

    countIndividualFeatures(dialog.selections.getFilePathGroups()[0], dialog.selections.getFilePathGroups()[1],
                            dialog.selections.getToggleStates()[0], checkForNonNegativeInteger(dialog.selections.getTextEntries()[0]))


if (__name__ == "__main__"): main()