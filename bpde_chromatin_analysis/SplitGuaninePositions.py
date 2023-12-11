# This script takes a file of BPDE damage positions and splits the guanine positions away from non-guanine bases.

import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from bpde_chromatin_analysis.helper_scripts.BPDE_DataDir import getDataDir


def splitGuaninePositions(damagePositionFilePaths: List[str]) -> List[str]:
    """
    This function takes a list of bed file paths for BPDE-damage positions and splits guanine positions away from non-guanine positions.
    Returns a list of the file paths for guanine positions.
    """

    guanineDamagePositionFilePaths = list()

    # Iterate through the given file paths.
    for damagePositionFilePath in damagePositionFilePaths:

        print(f"\nWorking with {os.path.basename(damagePositionFilePath)}")

        # Create the paths to the output files.
        basename = os.path.basename(damagePositionFilePath).rsplit(".bed",1)[0]
        guanineDamagePositionFilePath = os.path.join(os.path.dirname(damagePositionFilePath), basename + "_G_only.bed")
        nonGuanineDamagePositionFilePath = os.path.join(os.path.dirname(damagePositionFilePath), basename + "_not_G.bed")
        guanineDamagePositionFilePaths.append(guanineDamagePositionFilePath)

        # Split the input file across the output files.
        with open(damagePositionFilePath, 'r') as damagePositionFilePath, open(guanineDamagePositionFilePath, 'w') as guanineDamagePositionFile, \
             open(nonGuanineDamagePositionFilePath, 'w') as nonGuanineDamagePositionFile:
            for line in damagePositionFilePath:

                if line.split()[6] == "G": guanineDamagePositionFile.write(line)
                else: nonGuanineDamagePositionFile.write(line)

    return guanineDamagePositionFilePaths


def main():

    # Create the Tkinter UI
    with TkinterDialog(workingDirectory = getDataDir(), title = "Split Guanine Positions") as dialog:
        dialog.createMultipleFileSelector("Damage Position Files:", 0, "damage_pos.bed", ("Bed Files",".bed"))

    splitGuaninePositions(dialog.selections.getFilePathGroups()[0])


if __name__ == "__main__": main()