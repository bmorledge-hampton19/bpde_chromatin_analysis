# bpde_chromatin_analysis
Analyzing lung cancer mutations, BPDE damage, and BPDE repair with respect to genomic features like nucleosomes and transcription factors
***

## Table of Contents
1. [Project Overview](#project-overview)
2. [Dependencies](#dependencies)
3. [Preparing Data for Analysis](#preparing-data-for-analysis)
4. [Analysis and Figure Generation](#analysis-and-figure-generation)
5. [Data Availability](#data-availability)
***

## Project Overview

### Accompanying Literature
The code in this repository was used to produce the findings in [this publication](NO_LINK_YET). This paper contains additional insight into the analysis using this code.

### Background
Lung cancer mutations have an unusual [preference for linker DNA over nucleosomal DNA](https://doi.org/10.1016/j.cell.2018.10.004). While previously established [repair patterns](https://doi.org/10.1073/pnas.1706021114) do not explain this preference, more recent [damage data](https://doi.org/10.1021/acscentsci.2c01100) display a very similar periodic trend.

### Overview of Methodology
The [mutperiod package](https://github.com/bmorledge-hampton19/mutperiod) is used heavily to count mutations, damage, and repair relative to chromatin features (and even chromatin features relative to other chromatin features). Because of this, the bulk of this repository simply contains jupyter notebooks to format data for mutperiod and leverage its analysis pipelines. The rest of the repository primarily contains scripts for curating the data, calling nucleosome dyads in irregular settings, and plotting results.

### A Disclaimer...
Admittedly, the software in this repository was largely not designed with other users in mind. Regardless, the core analysis is reproducible after taking a few steps to set up a specific environment and obtain/format the necessary input files. In the event that the following steps are insufficient to reproduce the desired results, you are more than welcome to [post an issue](../../issues) or email me directly at b.morledge-hampton@wsu.edu.
***

## Dependencies
The most up-to-date versions of the following packages are recommended:
- [mutperiod](https://github.com/bmorledge-hampton19/mutperiod) (This analysis requires the most up-to-date version of mutperiod, which is not currently available through the related ppa and requires manual installation of the associated dependencies. If there's genuine interest in reproducing this analysis, I may be willing to update the ppa for an easier installation.)
- [chromatinfeaturesanalysis](https://github.com/bmorledge-hampton19/Chromatin_Features_Analysis) (Only the Python files are required and can be installed using the setup.py file at the root of the repository)
- The following tools must also be available through the command line (Recommend installing through apt, where possible):
  - [bedtools](https://bedtools.readthedocs.io/en/latest/)
  - [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - [samtools](http://www.htslib.org/)
  - [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) or [bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/)

## Preparing Data for Analysis

### Setting up data directories
The necessary data directories can be set up automatically using the [SetUpDataDir.py](/bpde_chromatin_analysis/helper_scripts/SetUpDataDir.py) script. Running the script will first prompt you to choose locations for this package's and mutperiod's parent data directories (`BPDE_data` and `mutperiod_data`), provided they do not already exist. These directories will then be automatically populated with the necessary sub-directories.

### Formatting input data.
All inputs need to be provided manually, with the proper names and in the proper directories. See the [Data Availability](#data-availability) section for more information on obtaining the data.

#### Lung cancer mutation data
- Mutation data should be copied in its original format to `BPDE_data/Alexandrov_LUAD/Lung_Adeno_clean_somatic_mutations_for_signature_analysis.txt`.

#### BPDE damage data
- Damage-seq reads will need to be aligned to the hg19 genome and converted to bed format. Adapter trimming is recommended prior to alignment. See [this paper](https://doi.org/10.1021/acscentsci.2c01100) for information on adapters. I recommend using the [AlignReads.py](https://github.com/bmorledge-hampton19/benbiohelpers/blob/main/python/benbiohelpers/Alignment/AlignReads.py) script in benbiohelpers to accomplish this.
- The 3 cellular and 3 naked-DNA samples with 2-micromolar exposure concentrations are the only samples that are analyzed by this pipeline.
- The 3 bed files from the aligned cellular samples should be copied to `BPDE_data/Jiang_BPDE_damage_maps/BEAS-2B_2uM_BPDE_cell_24h/BEAS-2B_2uM_BPDE_cell_24h_rep#.bed` where `#` is replaced by the repetition number.
- The 3 bed files from the aligned naked-DNA samples should be copied to `BPDE_data/Jiang_BPDE_damage_maps/BEAS-2B_2uM_BPDE_nDNA_24h/BEAS-2B_2uM_BPDE_nDNA_24h_rep#.bed` where `#` is replaced by the repetition number.
- Run the [prepocessing notebook](bpde_chromatin_analysis/_notebooks/Jiang_damage_maps/Jiang_damage_maps_preprocessing.ipynb) to finish formatting the data.

#### BPDE repair data
- tXR-seq reads will need to be aligned to the hg19 genome and converted to bed format. Adapter trimming is recommended prior to the alignment. See [this paper](http://www.genesdev.org/cgi/doi/10.1101/gad.261271.115) for information on adapters. I recommend using the [AlignReads.py](https://github.com/bmorledge-hampton19/benbiohelpers/blob/main/python/benbiohelpers/Alignment/AlignReads.py) script in benbiohelpers to accomplish this.
- Each of the two replicates for the data should be aligned separately, and the resulting bed files should be sorted alphabetically on the first column, numerically on the second and third columns, and alphabetically on the sixth column (prioritized in that order).
- The sorted files then need to be deduplicated based on the four sorted columns using the [RemoveDuplicates.py](https://github.com/bmorledge-hampton19/benbiohelpers/blob/main/python/benbiohelpers/FileSystemHandling/RemoveDuplicates.py) script from benbiohelpers (default parameters should suffice).
- The deduplicated files should be combined and sorted again.
- The resulting file should be copied to `BPDE_data/Li_tXR-seq/Li_tXR-seq.bed`

#### Nucleosome maps
- The nucleosome maps used here are tricky to recreate and are not currently hosted anywhere, but I hope this will change soon. In the meantime, I am happy to provide the original maps if you reach out to me at b.morledge-hampton@wsu.edu, or you are welcome to try to recreate them from the original DNase-seq and MNase-seq data using the methods described [here](NO_LINK_YET).
- The hybrid nucleosome map should be copied to `BPDE_data/relative_nucleosome_patterns/hybrid/hg19_hybrid_nuc_map.bed` and `mutperiod_data/__external_data/hg19/hg19_hybrid_nuc_map/hg19_hybrid_nuc_map.bed`.
- The LCL MNase-only map should be copied to `BPDE_data/relative_nucleosome_patterns/LCL_MNase/hg19_LCL_MNase_nuc_map.bed` and `mutperiod_data/__external_data/hg19/hg19_LCL_MNase_nuc_map/hg19_LCL_MNase_nuc_map.bed`.

#### Transcription factor binding sites
- The CTCF binding sites should be copied to `BPDE_data/relative_TFBS_patterns/CTCF/hg19_CTCF_known.bed` and `mutperiod_data/__external_data/hg19/hg19_CTCF_known/hg19_CTCF_known.bed`
- The SP1 binding sites should be copied to `BPDE_data/relative_TFBS_patterns/SP1/hg19_SP1_known.bed` and `mutperiod_data/__external_data/hg19/hg19_SP1_known/hg19_SP1_known.bed`

#### Transcription start sites
- The transcription start sites should be copied to `mutperiod_data/__external_data/hg19/hg19_protein_coding_genes_TSSs/hg19_protein_coding_genes_TSSs.bed`

#### Removing blacklisted regions
Using the benbiohelpers [RemoveBlacklistedRegions.py](https://github.com/bmorledge-hampton19/benbiohelpers/blob/main/python/benbiohelpers/FileSystemHandling/RemoveBlacklistedRegions.py) script, all nucleosome, transcription factor, and transcription start site files should be filtered using the ENCODE blacklisted regions in a 1000 bp radius. Note that all input files must be sorted alphabetically on the first column and numerically on the second and third columns (prioritized in that order). The transcription factor binding site and transcription start site files in the [maintained_data](/maintained_data/) folder of this repository have already been filtered this way.

***

## Analysis and Figure Generation
Once all the input data has been properly formatted, the following [jupyter notebooks](/bpde_chromatin_analysis/_notebooks) can be run exactly as they are to reproduce the analyses that contributed to [publication](NO_LINK_YET). Note that while other notebooks are available in this repository, they may require additional inputs and were not part of the core analysis.

### Lung cancer mutation analyses
- [Alexandrov_LUAD_mutation_TFBS_data_generation.ipynb](bpde_chromatin_analysis/_notebooks/Alexandrov_LUAD_mutations/Alexandrov_LUAD_mutation_TFBS_data_generation.ipynb)
- [Alexandrov_LUAD_mutation_TFBS_figure_generation.ipynb](bpde_chromatin_analysis/_notebooks/Alexandrov_LUAD_mutations/Alexandrov_LUAD_mutation_TFBS_figure_generation.ipynb)
- [Alexandrov_LUAD_mutation_nucleosome_periodicity_data_generation.ipynb](bpde_chromatin_analysis/_notebooks/Alexandrov_LUAD_mutations/Alexandrov_LUAD_mutation_nucleosome_periodicity_data_generation.ipynb)
- [Alexandrov_LUAD_mutation_nucleosome_periodicity_figure_generation.ipynb](bpde_chromatin_analysis/_notebooks/Alexandrov_LUAD_mutations/Alexandrov_LUAD_mutation_nucleosome_periodicity_figure_generation.ipynb)

### BPDE damage analyses
- [Jiang_damage_maps_TFBS_data_generation.ipynb](bpde_chromatin_analysis/_notebooks/Jiang_damage_maps/Jiang_damage_maps_TFBS_data_generation.ipynb)
- [Jiang_damage_maps_TFBS_figure_generation.ipynb](bpde_chromatin_analysis/_notebooks/Jiang_damage_maps/Jiang_damage_maps_TFBS_figure_generation.ipynb)
- [Jiang_damage_maps_TSS_data_generation.ipynb](bpde_chromatin_analysis/_notebooks/Jiang_damage_maps/Jiang_damage_maps_TSS_data_generation.ipynb)
- [Jiang_damage_maps_TSS_figure_generation.ipynb](bpde_chromatin_analysis/_notebooks/Jiang_damage_maps/Jiang_damage_maps_TSS_figure_generation.ipynb)
- [Jiang_damage_maps_nucleosome_periodicity_data_generation.ipynb](bpde_chromatin_analysis/_notebooks/Jiang_damage_maps/Jiang_damage_maps_nucleosome_periodicity_data_generation.ipynb)
- [Jiang_damage_maps_nucleosome_periodicity_figure_generation.ipynb](bpde_chromatin_analysis/_notebooks/Jiang_damage_maps/Jiang_damage_maps_nucleosome_periodicity_figure_generation.ipynb)

### BPDE repair analyses
- [Li_tXR-seq_TFBS_data_generation.ipynb](bpde_chromatin_analysis/_notebooks/Li_tXR-seq/Li_tXR-seq_TFBS_data_generation.ipynb)
- [Li_tXR-seq_TFBS_figure_generation.ipynb](bpde_chromatin_analysis/_notebooks/Li_tXR-seq/Li_tXR-seq_TFBS_figure_generation.ipynb)
- [Li_tXR-seq_TSS_data_generation.ipynb](bpde_chromatin_analysis/_notebooks/Li_tXR-seq/Li_tXR-seq_TSS_data_generation.ipynb)
- [Li_tXR-seq_TSS_figure_generation.ipynb](bpde_chromatin_analysis/_notebooks/Li_tXR-seq/Li_tXR-seq_TSS_figure_generation.ipynb)
- [Li_tXR-seq_nucleosome_periodicity_data_generation.ipynb](bpde_chromatin_analysis/_notebooks/Li_tXR-seq/Li_tXR-seq_nucleosome_periodicity_data_generation.ipynb)
- [Li_tXR-seq_nucleosome_periodicity_figure_generation.ipynb](bpde_chromatin_analysis/_notebooks/Li_tXR-seq/Li_tXR-seq_nucleosome_periodicity_figure_generation.ipynb)

### Transcription-factor-relative and nucleosome-relative analyses
- [relative_TFBS_patterns_data_generation.ipynb](bpde_chromatin_analysis/_notebooks/relative_TFBS_patterns/relative_TFBS_patterns_data_generation.ipynb)
- [relative_TFBS_patterns_figure_generation.ipynb](bpde_chromatin_analysis/_notebooks/relative_TFBS_patterns/relative_TFBS_patterns_figure_generation.ipynb)
- [relative_nucleosome_patterns_data_generation.ipynb](bpde_chromatin_analysis/_notebooks/relative_nucleosome_patterns/relative_nucleosome_patterns_data_generation.ipynb)
- [relative_nucleosome_patterns_figure_generation.ipynb](bpde_chromatin_analysis/_notebooks/relative_nucleosome_patterns/relative_nucleosome_patterns_figure_generation.ipynb)

***

## Data Availability
- Lung cancer mutation data: `ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl` through [https://doi.org/10.1038/nature12477](https://doi.org/10.1038/nature12477)
- BPDE damage data: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224001](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224001)
- BPDE repair data: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97675](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97675)
- MNase-seq data: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36979)
- DNase-seq data: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31388](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31388)
- CTCF binding sites: [/maintained_data/hg19_CTCF_known.bed](/maintained_data/hg19_CTCF_known.bed)
- SP1 binding sites: [/maintained_data/hg19_SP1_known.bed](/maintained_data/hg19_SP1_known.bed)
- Transcription start sites: [/maintained_data/hg19_protein_coding_genes_TSSs.bed](/maintained_data/hg19_protein_coding_genes_TSSs.bed)
- ENCODE blacklist: [https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz)
