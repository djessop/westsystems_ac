# westsystems_ac

Utilities for the processing of data files produced by WestSystems accumulation chambers software (FluxManager: https://www.westsystems.com/instruments/download/).  These utilities provide functionality similar to FluxRevision, and allow for molar fluxes (i.e. units of mol/m2/day) to be calculated from the slope of concentration as a function of time plots.

## Usage

The main utility is a class named WestsystemsFile, which reads and processes the raw data from gas concentration measurements at a site.

```python
from westsystems_ac.read_westsystems import WestsystemsFile

wsf = WestsystemsFile('example_000_01011970_000000.txt',
		       gas_species='CO2', 
		       ac_chamber='B')
```

The WestsystemsFile class will parse the data file "example_000_01011970_000000.txt" and display a plot of the concentration data for the chosen gas species a function of time.  The user can then select the portion of the data that he wishes using clicks of the left mouse button.  Erroneous clicks can be undone using the right mouse button.  After the region is selected, it is shown in red on the plot and the user is invited to validate ('y') the data or reselect the region ('n') by typing a response in the console.  This process will continue until 'y' has been entered.  An image of the the validated data plot will then be saved in the same directory as the data file.  

Simultaneously, a pandas.DataFrame is created (one entry) as a class attribute.

## Batch processing 

The following code will batch process all files in the directory DIR_NAME, treating CO2 data for accumulation chamber type B

```bash
python read_westsystems.py DIR_NAME CO2 B
```

## Requirements

This packaged depends on, and has been tested using the following pagages/machine configuration:
```bash
numpy     : 1.26.4
scipy     : 1.16.0
matplotlib: 3.10.5
utm       : 0.7.0
pandas    : 2.2.0
sklearn   : 1.5.2

Python implementation: CPython
Python version       : 3.12.3
IPython version      : 9.6.0

Compiler    : GCC 13.3.0
OS          : Linux
Release     : 6.8.0-90-generic
Machine     : x86_64
Processor   : x86_64
CPU cores   : 16
Architecture: 64bit
```

## TO DO:
- read all metadata, including sensor data