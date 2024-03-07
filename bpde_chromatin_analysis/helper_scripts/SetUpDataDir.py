# This script is meant to be used by third parties interested in reproducing this analysis.
# It creates the necessary directory trees for the data used in the Jupyter notebooks.
# Note that the actual data still needs to be supplied to the directories.

import os
from bpde_chromatin_analysis.helper_scripts.BPDE_DataDir import getDataDir
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getExternalDataDirectory as getMutperiodExternalDataDirectory
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs

def setUpDataDir():

    # Create BPDE data directories
    BPDE_DataDir = getDataDir()

    mutationDataDir = os.path.join(BPDE_DataDir, "Alexandrov_LUAD")
    damageDataDirs = [os.path.join(BPDE_DataDir, "Jiang_BPDE_damage_maps", dataSet) for dataSet in ("BEAS-2B_2uM_BPDE_cell_24h", "BEAS-2B_2uM_BPDE_nDNA_24h")]
    repairDataDir = os.path.join(BPDE_DataDir, "Li_tXR-seq")
    nucleosomeRelativeDataDirs = [os.path.join(BPDE_DataDir, "relative_nucleosome_patterns", nucMap) for nucMap in ("hybrid", "LCL_MNase")]
    TFBS_RelativeDataDirs = [os.path.join(BPDE_DataDir, "relative_TFBS_patterns", tf) for tf in ("CTCF", "SP1")]

    checkDirs(mutationDataDir, *damageDataDirs, repairDataDir, *nucleosomeRelativeDataDirs, *TFBS_RelativeDataDirs)

    # Create mutperiod data directories
    mutperiodHg19Dir = os.path.join(getMutperiodExternalDataDirectory(), "hg19")

    mutperiodHybridNucMapDir = os.path.join(mutperiodHg19Dir, "hg19_hybrid_nuc_map")
    mutperiodLCL_MNaseNucMapDir = os.path.join(mutperiodHg19Dir, "hg19_LCL_MNase_nuc_map")
    mutperiodCTCFDir = os.path.join(mutperiodHg19Dir, "hg19_CTCF_known")
    mutperiodSP1Dir = os.path.join(mutperiodHg19Dir, "hg19_SP1_known")
    mutperiodTSSDir = os.path.join(mutperiodHg19Dir, "hg19_protein_coding_genes_TSSs")

    checkDirs(mutperiodHybridNucMapDir, mutperiodLCL_MNaseNucMapDir, mutperiodCTCFDir, mutperiodSP1Dir, mutperiodTSSDir)

    def createAppendToDataNameFile(directory, dataNameToAppend):
        with open(os.path.join(directory, "append_to_data_name.txt"), 'w') as appendToDataNameFile:
            appendToDataNameFile.write(dataNameToAppend)
    
    createAppendToDataNameFile(mutperiodHybridNucMapDir, "_hybrid_nuc_map")
    createAppendToDataNameFile(mutperiodLCL_MNaseNucMapDir, "_LCL_MNase_map")
    createAppendToDataNameFile(mutperiodCTCFDir, "_CTCF_known")
    createAppendToDataNameFile(mutperiodSP1Dir, "_SP1_known")
    createAppendToDataNameFile(mutperiodTSSDir, "_TSSs")


if __name__ == "__main__": setUpDataDir()