# prad2_BackgroundContamination
This software aims to quickly find the e-e and e-p Yields from PRad-II experimental data while also being able to compare runs from the 4 separate different run types. Using different combinations of these run types, one can find the various relative contributions of different background contamination source.

![Alt Text](./Images/PRad_SpecialBackgroundRuns.png)

|Run Type|Description|Gas in Cell?|Residual gas in chamber?|Target cell in place?|
|:------:|:---------:|:----------:|:-----------------------|:--------------------|
|(a)|Production configuration.|✔|✔|✔|
|(b)|Residual gas measurement.|X|✔|✔|
|(c)|No gas at all. Cell Measurement.|X|X|✔|
|(d)|No Gas, No Cell. Downstream measurement.|X|X|X|

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