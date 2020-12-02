## Data description
Enzyme activity data is in the titular “Enzyme activity” data folder. Enzyme assays were conducted by diluting 7 substrates for hydrolytic enzymes and 1 substrate for oxidative enzymes. The litter from each sample is homogenized and then plated with varying concentrations of each substrates and left to incubate and then read in the old plate reader (BioTek Synergy 4). There are 2 types of plates, each for a different type of enzymes: black plates for hydrolytic enzymes, and clear plates for oxidative enzymes. Hydrolytic enzyme activity was obtained by reading the fluorescence from each well in a sample black plate while oxidative enzyme activity was obtained from reading the absorbance from a sample clear plate. Please refer to the ***2020-11-24-Plate-Layouts.xlsx*** spreadsheet to understand the plate layout. To understand how assays are conducted, please refer to the ***Enzyme Assays Description.md*** file.<br>

Enzyme activity data are .txt files. A data file with “Black-Plates” in its title only contains data of hydrolytic enzymes, while a data file with “Clear-Plates” in its title only contains data of oxidative enzymes. In addition, each of these .txt files also contain the timepoint in its title (T0, T3, T5, or T6).
## Columns description (in order from left to right in a data file)
**Well**: the identifier of a particular well; all wells are plated with the same sample, but each well has a unique concentration of a substrate.
- Unitless

**Long sample name**: the identifier of a particular sample and plate. Each plate only contains one sample, and each plate has 96 wells, so for a particular sample, there should be 96 rows with unique Well IDs. Each of these 96 rows will have an identical Plate ID as each other. This column contains the following useful pieces of information: (1) the date the assay was performed in YYMMDD, (2) a unique *sample identifier* that is described in the descriptive metadata file **2020-11-10-Litter-bag-codes.xlsx** (3) and if the file is a “Black-Plates” file, a letter denoting the plate type (B for buffer plates and no letter if the plate is a sample plate).
- Unitless

**Fluorescence or absorbance reading**: this is the raw fluorescence (for “Black-Plates” files) or absorbance (for “Clear-Plates” files) reading for each well.
- Fluorescence/absorbance units


## Association metadata files
**2020-11-10-Litter-bag-codes.xlsx**: this association metadata file describes the unique sample identifier in the Plate ID column in the raw enzyme activity and data files. This is available under the “Metadata” folder.<br>

**2020-11-10-Dry-weight.xlsx**: this file describes the month and year (YY) of each time point in which litterbags were picked up from the field. Although this file contains metadata, it is primarily a data file, hence its location under the “Dry weights” folder.

**2020-11-24-Plate-Layouts.xlsx**: this describes the plate layout of black plates (which assays for hydrolytic enzymes) and clear plates (which assays for oxidative enzymes). Included is the highest concentration of each substrate. Each substrate is diluted serially by half with the top row having the highest substrate concentration while the bottom row has the lowest substrate concentration. This file is available under the “Metadata” folder.

**Enzyme Assays Description.md**: this file describes how the enzyme assays are conducted.

## File/data owner
Steven Allison


## Who is permitted to access
All project personnel who maintains/conducts research at the Loma Ridge GLobal Change Experiment


## Permission to update/edit/use data
Please contact the data owner (Steven Allison) for permission to change or analyze the data.

## Works cited
German, Donovan P. et al. “Optimization of hydrolytic and oxidative enzyme methods for ecosystem studies.” *Soil Biology & Biogeochemistry* vol. 43, no. 7, 2011, pp. 1387-1397.
