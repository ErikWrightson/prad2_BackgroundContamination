# prad2_BackgroundContamination
This software aims to quickly find the e-e and e-p Yields from PRad-II experimental data while also being able to compare runs from the 4 separate different run types. Using different combinations of these run types, one can find the various relative contributions of different background contamination source.

![Alt Text](./Images/PRad_SpecialBackgroundRuns.png)

|Run Type|Description|Gas in Cell?|Residual gas in chamber?|Target cell in place?|
|:------:|:---------:|:----------:|:-----------------------|:--------------------|
|(a)|Production configuration.|✔|✔|✔|
|(b)|Residual gas measurement.|X|✔|✔|
|(c)|No gas at all. Cell Measurement.|X|X|✔|
|(d)|No Gas, No Cell. Downstream measurement.|X|X|X|

## Overview
This program quickly finds rough e-e and e-p yields using various cuts for different types of runs and then compares them if at least 1 run of all 4 types is present. If not all types are present only the yields are calculated for the present types. All yields are normalized to by the liveCharge collected which is the livetime times the charge collected during the run which approximates for the amount of uptime of the beam. Finer studies would require cutting of downtime periods.

Comparisons are then made to find the amount of contamination from each non-target element. Here is a list of the various subtractions that can be done:
|Subtraction|What is left?|
|:---------:|:-----------:|
|(a)-(b)|H2 Gas Contribution (used as the general signal total)|
|(b)-(c)|Residual Gas outside the cell within the target chamber|
|(c)-(d)|Contribution of the cell itself.|
|(d)|Contribution of anything outside the target chamber, notably the collimators.|

Each of these is then divided by the liveCharge normalized yield of (a)-(b) to get the level of contamintation as compared to the good signal.

## Settup on ifarm
```console
source farm_setup.csh
make
```

## Run Options
|Option|Input|Description|
|:----:|:----:|:--------:|
|-a|<TypeA_File/List>|Marks the accompanying file as a Type A file or a txt file containing a list of Type A files|
|-b|<TypeB_File/List>|Marks the accompanying file as a Type B file or a txt file containing a list of Type A files|
|-c|<TypeC_File/List>|Marks the accompanying file as a Type C file or a txt file containing a list of Type C files|
|-d|<TypeD_File/List>|Marks the accompanying file as a Type D file or a txt file containing a list of Type D files|
|-f|<FileName (NO EXTENSION)>|Give the custom filename that you would like to be used for the aggregated output ROOT file and PDF|
|-D|<liveChargeDb>|Custom location for the liveCharge .dat file (in the format described below)|
|-L|N/A|Indicates that list .txt files for handing in many run number locations at once will be used (in format described below)|
|-v|N/A|Verbosity. Fill, print and save all histograms before and after each cut is applied.|
|-m|<nThreads>|Multithreading. Number of threads to request.|
|-G|N/A|Gems. Indicates that Møller center finding should be done with the GEM planes and includea GEM matching cut as part of the expected energy cut.|
|-h|N/A|Help. Brings up options helper menu.|

**NOTE: At least one file type -a, -b, -c or -d must be used.**

## Input File Formats
### Individual File Inputs
If the -L option is not used then for the type A-D inputs a .root file is expected in the format of the PRad-II reconstructed replay.
### File List
If the -L option is used then for the A-D type inputs a .txt file is expected in the format of all entire file path including filename of each file for this type is included separated by linebreaks.
### LiveCharge Dat File
The LiveCharge .dat file is expected to come in the format of the example (./database/beam_charge.dat), but more specifically the run number must be the first item in the line and the liveCharge must be the 4th.
