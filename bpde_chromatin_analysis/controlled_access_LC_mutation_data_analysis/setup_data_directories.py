import os
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory as getMutperiodDataDir
from bpde_chromatin_analysis.helper_scripts.BPDE_DataDir import getDataDir as getBPDEDataDir
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs

mutperiodDataDir = getMutperiodDataDir()
bpdeDataDir = getBPDEDataDir()

checkDirs(os.path.join(bpdeDataDir, "controlled_access_LC"))