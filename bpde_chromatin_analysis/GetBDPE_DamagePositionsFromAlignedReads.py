# This script derives BPDE damage positions from aligned sequencing reads in bed format.
# It also records the nucleotide at the damaged position.
# This script expects data prepared as in the 2023 Jiang et al. paper: https://pubs.acs.org/doi/10.1021/acscentsci.2c01100

import os
from typing import List
from benbiohelpers.FileSystemHandling.AddSequenceToBed import addSequenceToBed
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from bpde_chromatin_analysis.helper_scripts.BPDE_DataDir import getDataDir


def getBPDE_DamagePositionsFromAlignedReads(alignedReadsFilePaths: List[str], genomeFastaFilePath, isSingleEnd = False) -> List[str]:
    """
    This function takes a list of bed file paths for formatted damage reads and
    converts them to bed files of the expected BPDE lesion positions (-1 relative to 5' end)
    """

    BPDE_PositionsFilePaths = list()

    # Iterate through the given file paths.
    for alignedReadsFilePath in alignedReadsFilePaths:

        print(f"\nWorking with {os.path.basename(alignedReadsFilePath)}")

        # Create the path to the output file.
        basename = os.path.basename(alignedReadsFilePath).rsplit(".bed",1)[0]
        BPDE_PositionsFilePath = os.path.join(os.path.dirname(alignedReadsFilePath), basename + "_BPDE_damage_pos.bed")
        BPDE_PositionsFilePaths.append(BPDE_PositionsFilePath)

        # Create the intermediate file with the expected BPDE positions
        print("Retrieving expected BPDE positions...")
        with open(alignedReadsFilePath, 'r') as alignedReadsFile, open(BPDE_PositionsFilePath, 'w') as BPDE_PositionsFile:
            for line in alignedReadsFile:

                splitLine = line.split()

                # Make sure we only work with the first read in each pair
                if not isSingleEnd and splitLine[3].endswith("/2"):
                    print("WARNING: Found read 2 in paired-end read. Are you sure the reads have been formatted?")
                    continue

                # Adjust the positions to fit the expected BPDE position
                if splitLine[5] == '+':
                    splitLine[2] = splitLine[1]
                    splitLine[1] = str(int(splitLine[1]) - 1)
                else:
                    splitLine[1] = splitLine[2]
                    splitLine[2] = str(int(splitLine[2]) + 1)

                # Drop unnecessary columns
                splitLine[3] = '.'
                splitLine[4] = '.'

                # Write the new line to the output file.
                BPDE_PositionsFile.write('\t'.join(splitLine)+'\n')

        # Add information on the nucleotide base present at the expected BPDE position.
        print("Retrieving bases at expected BPDE positions...")
        addSequenceToBed(BPDE_PositionsFilePath, genomeFastaFilePath)

    return BPDE_PositionsFilePaths


def main():

    # Create the Tkinter UI
    with TkinterDialog(workingDirectory = getDataDir(), title = "Get BPDE Positions From Aligned Reads") as dialog:
        dialog.createMultipleFileSelector("Aligned Reads Files:", 0, "formatted.bed", ("Bed Files",".bed"))    
        dialog.createGenomeSelector(1, 0)
        dialog.createCheckbox("Input is single-end", 2, 0)

    getBPDE_DamagePositionsFromAlignedReads(dialog.selections.getFilePathGroups()[0], dialog.selections.getGenomes(returnType="fasta")[0],
                                            dialog.selections.getToggleStates()[0])


if __name__ == "__main__": main()