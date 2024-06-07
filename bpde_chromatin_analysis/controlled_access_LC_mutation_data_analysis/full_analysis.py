import os
from bpde_chromatin_analysis.helper_scripts.BPDE_DataDir import getDataDir
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getExternalDataDirectory as getMutperiodExternalDataDirectory, DataTypeStr
from mutperiodpy.input_parsing.ParseCustomBed import parseCustomBed
from mutperiodpy.RunAnalysisSuite import runAnalysisSuite
from mutperiodpy.RunNucleosomeMutationAnalysis import runNucleosomeMutationAnalysis
from benbiohelpers.FileSystemHandling.DirectoryHandling import getFilesInDirectory
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.OutputDataStratifiers import AmbiguityHandling


# Prepare paths to relevant data files.
mutperiodExternalDataDirectory = getMutperiodExternalDataDirectory()
mutperiodHg19Directory = os.path.join(mutperiodExternalDataDirectory, "hg19")
hg19FastaFilePath = os.path.join(mutperiodHg19Directory, "hg19.fa")

lcMutationDataDir = os.path.join(getDataDir(), "controlled_access_LC")
lcMutationInputDataFilePath = os.path.join(lcMutationDataDir, "controlled_access_LC_input_data.bed")

lcMutationMutperiodInputFilePaths = parseCustomBed([lcMutationInputDataFilePath], hg19FastaFilePath, onlySingleBaseSubs = True)


# Run nucleosome periodicity analysis
nucleosomeMapFilePaths = list()
nucleosomeMapFilePaths.append(os.path.join(mutperiodHg19Directory, "hg19_hybrid_nuc_map", "hg19_hybrid_nuc_map.bed"))
nucleosomeMapFilePaths.append(os.path.join(mutperiodHg19Directory, "hg19_LCL_MNase_nuc_map", "hg19_LCL_MNase_nuc_map.bed"))
nucleosomeMapFilePaths.append(os.path.join(mutperiodHg19Directory, "hg19_NHF1_MNase_nuc_map", "hg19_NHF1_MNase_nuc_map.bed"))

nucleosomeMapNames = [os.path.basename(nucleosomeMapFilePath).rsplit('.', 1)[0] for nucleosomeMapFilePath in nucleosomeMapFilePaths]

runAnalysisSuite(lcMutationMutperiodInputFilePaths, nucleosomeMapNames, normalizationMethod = "Trinuc/Quadrunuc", customBackgroundDir = None,
                 useSingleNucRadius = True, includeLinker = False, useNucGroupRadius = True)
lcMutationMutperiodInputFilePaths = getFilesInDirectory(lcMutationDataDir, DataTypeStr.mutations + ".bed")

rawNucleosomeCountsFilePaths = list()
normalizedNucleosomeCountsFilePaths = list()
for nucleosomeMapName in nucleosomeMapNames:
    normalizedNucleosomeCountsFilePaths += getFilesInDirectory(os.path.join(lcMutationDataDir, nucleosomeMapName), DataTypeStr.normNucCounts + ".tsv")
    rawNucleosomeCountsFilePaths += getFilesInDirectory(os.path.join(lcMutationDataDir, nucleosomeMapName), DataTypeStr.rawNucCounts + ".tsv")

runNucleosomeMutationAnalysis(normalizedNucleosomeCountsFilePaths,
                              outputFilePath = os.path.join(lcMutationDataDir, "controlled_access_LC_mutations_normalized_periodicity_data.tsv"),
                              overridePeakPeriodicityWithExpected = False, alignStrands = False)
runNucleosomeMutationAnalysis(rawNucleosomeCountsFilePaths,
                              outputFilePath = os.path.join(lcMutationDataDir, "controlled_access_LC_mutations_raw_periodicity_data.tsv"),
                              overridePeakPeriodicityWithExpected = False, alignStrands = False)


# Run TFBS analysis
TFBS_FilePaths = list()
TFBS_FilePaths.append(os.path.join(mutperiodHg19Directory, "hg19_CTCF_known", "hg19_CTCF_known.bed"))
TFBS_FilePaths.append(os.path.join(mutperiodHg19Directory, "hg19_SP1_known", "hg19_SP1_known.bed"))

TFBS_Names = [os.path.basename(TFBS_FilePath).rsplit('.', 1)[0] for TFBS_FilePath in TFBS_FilePaths]

runAnalysisSuite(lcMutationMutperiodInputFilePaths, TFBS_Names, normalizationMethod = "Trinuc/Quadrunuc", customBackgroundDir = None,
                 useSingleNucRadius = True, includeLinker = False, useNucGroupRadius = True, useNucStrand = True)


# Run TSS analysis
TSS_FilePaths = list()
TSS_FilePaths.append(os.path.join(mutperiodHg19Directory, "hg19_protein_coding_genes_TSSs", "hg19_protein_coding_genes_TSSs.bed"))

TSS_Names = [os.path.basename(TSS_FilePath).rsplit('.', 1)[0] for TSS_FilePath in TSS_FilePaths]

runAnalysisSuite(lcMutationMutperiodInputFilePaths, TSS_Names, normalizationMethod = "Trinuc/Quadrunuc", customBackgroundDir = None,
                 useSingleNucRadius = True, includeLinker = False, useNucGroupRadius = True, useNucStrand = True)


# Run transcriptional asymmetry analysis
def binInGenes(featureFilePaths, geneDesignationsFilePath, flankingBinSize = 0, flankingBinNum = 0, 
               filePathSuffix = "", colorColIndex = None, geneFractionNum = 6):
    """
    Count features (e.g., mutations) on the transcribed and nontranscribed strands of genes and bin them across 6 gene fractions.
    The flankingBinSize and flankingBinNum parameters add additional bins of a constant length on the regions flanking genes (on each side). Importantly,
    these regions must already be a part of the regions given in the gene designations file.
    """

    class BinInGenesCounter(ThisInThatCounter):

        def setUpOutputDataHandler(self):
            super().setUpOutputDataHandler()
            if colorColIndex is not None:
                self.outputDataHandler.addSimpleEncompassingColStratifier(outputName = "Color_Domain", colIndex = colorColIndex)
            self.outputDataHandler.addFeatureFractionStratifier(outputName = "Gene_Fraction", fractionNum = geneFractionNum,
                                                                flankingBinSize = flankingBinSize, flankingBinNum = flankingBinNum)
            self.outputDataHandler.addStrandComparisonStratifier(strandAmbiguityHandling = AmbiguityHandling.tolerate)
            if colorColIndex is None: customStratifyingNames=(None, {True:"Coding_Strand_Counts", False:"Noncoding_Strand_Counts"})
            else: customStratifyingNames=(None, None, {True:"Coding_Strand_Counts", False:"Noncoding_Strand_Counts"})
            self.outputDataHandler.createOutputDataWriter(self.outputFilePath, customStratifyingNames = customStratifyingNames)

    for featureFilePath in featureFilePaths:

        print("\nWorking in", os.path.basename(featureFilePath))

        outputFilePath = featureFilePath.rsplit('.', 1)[0] + "_gene_bins"
        if filePathSuffix: outputFilePath += '_' + filePathSuffix
        outputFilePath += ".tsv"
        metadataFilePath = outputFilePath.rsplit('.',1)[0] + ".metadata"

        counter = BinInGenesCounter(featureFilePath, geneDesignationsFilePath, outputFilePath)
        counter.count()

        with open(metadataFilePath, 'w') as metadataFile:
            metadataFile.write(f"Flanking_Bin_Size:\t{flankingBinSize}\n")
            metadataFile.write(f"Flanking_Bins_Each_Side:\t{flankingBinNum}\n")
            metadataFile.write(f"Gene_Designations_File_Path:\t{geneDesignationsFilePath}\n")

hg19GeneDesignationsFilePath = os.path.join(mutperiodHg19Directory, "hg19_protein_coding_genes",
                                            "hg19_protein_coding_genes_32658_bp_extended_ENCODE_blacklist_filtered.bed")

binInGenes(lcMutationMutperiodInputFilePaths, hg19GeneDesignationsFilePath, filePathSuffix = "TS_vs_NTS", flankingBinNum = 3, flankingBinSize = 10886)
binInGenes(lcMutationMutperiodInputFilePaths, hg19GeneDesignationsFilePath, filePathSuffix = "TS_vs_NTS_fine-grain", flankingBinNum = 300, flankingBinSize = 108,
           geneFractionNum = 600)