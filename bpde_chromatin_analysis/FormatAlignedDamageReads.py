# This script formats aligned damage reads for later steps in the bpde_chromatin_analysis pipeline.
# (See function docstring for further info)
import os, subprocess, shutil
from typing import List, Dict
from benbiohelpers.Alignment.CombinePairedBedReads import combinePairedBedReads
from benbiohelpers.FileSystemHandling.RemoveDuplicates import removeDuplicates


def formatAlignedDamageReads(alignedDamageReadsFilePaths: List[str]):
    """
    Formats aligned damage reads for later steps in the bpde_chromatin_analysis pipeline.
    Formatting includes the following (in order):
    - Paired-end reads are combined where possible
    - Lone "2nd pair" reads are removed,
    - Duplicate reads are removed,
    - Repetitions are combined
    
    NOTE: Do not combine multiple data sets in the same directory or repetition combining will fail.
    
    Returns a list of formatted files.
    """

    formattedDamageReadsFilePaths = list()
    noDupsFilePathsByDataDirectory: Dict[str, List[str]] = dict()

    for alignedDamageReadsFilePath in alignedDamageReadsFilePaths:

        print(f"\nWorking in {os.path.basename(alignedDamageReadsFilePath)}...")

        # Combine the paired-end reads.
        print("Combining paired end reads...")
        combinedReadsFilePath = combinePairedBedReads([alignedDamageReadsFilePath], outputToTmpDir = True, verbose = False)[0]

        # Remove lone "2nd pair" reads.
        print("Removing lone \"2nd pair\" reads...")
        read2FilteredFilePath = combinedReadsFilePath.rsplit(".bed",1)[0] + "_read_2_filtered.bed"
        with open(combinedReadsFilePath, 'r') as combinedReadsFile, open(read2FilteredFilePath, 'w') as read2FilteredFile:
            for line in combinedReadsFile:
                if not line.split()[3].endswith('2'): read2FilteredFile.write(line)

        # Sort file and remove duplicates (sorting is necessary for duplicate removal)
        print("Sorting...")
        subprocess.check_call(("sort", "-k1,1", "-k2,2n", "-k3,3n", "-k6,6", "-o", read2FilteredFilePath, read2FilteredFilePath))
        print("Removing duplicates...")
        noDupsFilePath = removeDuplicates([read2FilteredFilePath], verbose = False)[0]
        if os.path.dirname(alignedDamageReadsFilePath) not in noDupsFilePathsByDataDirectory:
            noDupsFilePathsByDataDirectory[os.path.dirname(alignedDamageReadsFilePath)] = list()
        noDupsFilePathsByDataDirectory[os.path.dirname(alignedDamageReadsFilePath)].append(noDupsFilePath)

    # Combine replicates
    print()
    for dataDirectory in noDupsFilePathsByDataDirectory:

        print(f"Combining replicates in directory {os.path.basename(dataDirectory)}...")
        # If there was only one replicate, just rename it.
        if len(noDupsFilePathsByDataDirectory[dataDirectory]) == 1:
            formattedDamageReadsFilePath = os.path.join(dataDirectory,
                                                        os.path.basename(noDupsFilePathsByDataDirectory[dataDirectory][0]).replace(
                                                            "combined_reads_read_2_filtered_no_dups.bed", "formatted.bed"
                                                        ))
            shutil.copyfile(noDupsFilePathsByDataDirectory[dataDirectory][0], formattedDamageReadsFilePath)
            formattedDamageReadsFilePaths.append(formattedDamageReadsFilePath)
        # Otherwise, combine reps. (Assumes that the input file paths end in "rep*")
        else:
            formattedDamageReadsFilePath = os.path.join(dataDirectory,
                                                        os.path.basename(noDupsFilePathsByDataDirectory[dataDirectory][0]).rsplit(
                                                            "combined_reads_read_2_filtered_no_dups.bed",1
                                                        )[0].rsplit("rep",1)[0] + "all_reps_formatted.bed")
            with open(formattedDamageReadsFilePath, 'w') as formattedDamageReadsFile:
                for noDupsFilePath in noDupsFilePathsByDataDirectory[dataDirectory]:
                    with open(noDupsFilePath, 'r') as noDupsFile:
                        shutil.copyfileobj(noDupsFile, formattedDamageReadsFile)
            formattedDamageReadsFilePaths.append(formattedDamageReadsFilePath)

    return formattedDamageReadsFilePaths