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
|>> source farm_setup.csh|
|>> make|
```