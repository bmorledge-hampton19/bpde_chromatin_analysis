# This script takes a file of BPDE damage positions and splits the guanine positions away from non-guanine bases.

import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.InputParsing.ParseToIterable import parseToIterable
from bpde_chromatin_analysis.helper_scripts.BPDE_DataDir import getDataDir


def filterOnSpecificNucleotides(damagePositionFilePaths: List[str], nucleotides: List[str] = ['G']) -> List[str]:
    """
    This function takes a list of bed file paths for BPDE-damage positions and splits entries containing specific nucleotides from those that don't.
    Returns a list of the file paths for the files containing entries with the given nucleotides.
    """

    guanineDamagePositionFilePaths = list()

    # Iterate through the given file paths.
    for damagePositionFilePath in damagePositionFilePaths:

        print(f"\nWorking with {os.path.basename(damagePositionFilePath)}")

        # Create the paths to the output files.
        basename = os.path.basename(damagePositionFilePath).rsplit(".bed",1)[0]
        guanineDamagePositionFilePath = os.path.join(os.path.dirname(damagePositionFilePath), basename + "_" + '_'.join(nucleotides) + "_only.bed")
        nonGuanineDamagePositionFilePath = os.path.join(os.path.dirname(damagePositionFilePath), basename + "_not_" + '_'.join(nucleotides) + ".bed")
        guanineDamagePositionFilePaths.append(guanineDamagePositionFilePath)

        # Split the input file across the output files.
        with open(damagePositionFilePath, 'r') as damagePositionFilePath, open(guanineDamagePositionFilePath, 'w') as guanineDamagePositionFile, \
             open(nonGuanineDamagePositionFilePath, 'w') as nonGuanineDamagePositionFile:
            for line in damagePositionFilePath:

                if line.split()[6] in nucleotides: guanineDamagePositionFile.write(line)
                else: nonGuanineDamagePositionFile.write(line)

    return guanineDamagePositionFilePaths


def main():

    # Create the Tkinter UI
    with TkinterDialog(workingDirectory = getDataDir(), title = "Split Guanine Positions") as dialog:
        dialog.createMultipleFileSelector("Damage Position Files:", 0, "damage_pos.bed", ("Bed Files",".bed"))
        dialog.createTextField("Nucleotide(s): ", 1, 0, defaultText = "G, C")

    filterOnSpecificNucleotides(dialog.selections.getFilePathGroups()[0],
                                parseToIterable(dialog.selections.getTextEntries()[0]))


if __name__ == "__main__": main()