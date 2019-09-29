# UNISim-tool
A Unified Simulation tool for the Rudjer Boskovic Institute nuclear clustering group.

Contributors: Daniele Dell'Aquila

Rudjer Boskovic Institute, Zagreb, Croatia

This software is a unified simulation tool for the nuclear clustering group of the Rudjer Boskovic Institute. It is able to reliably run Monte Carlo simulations involving different physics models (such as sequential decay reactions) and an arbitrary configuration of detectors (DSSSDs, Lamp Wedge detectors). This software contains a high-precision energy loss program extracted by LISE++ for the calculation of the energy loss of charged particles in the matter. Another interesting novelty of this software is the advanced graphical mode for the visualization of the reaction and the detection setup, this is based on the Event Visualization Environment of ROOT.

# Useful links:
  * [RUDER_BOSKOVIC](https://www.irb.hr/eng) : Rudjer Boskovic Institute
  
Table of contents
=================
<!--ts-->
* [Getting the code](#getting-the-code)
  * [Using git](#using-git)
  * [Downloading from Git Hub](#downloading-from-git-hub)
* [Setup and Configuration](#setup-and-configuration)
  * [Requirements](#requirements)
  * [Compile the Code](#compile-the-code)
  * [Configure the Program](#configure-the-program)
  * [Limitations](#limitations)
* [The UNISim-tool environment](#the-unisim-tool-environment)
  * [Run the Code](#run-the-code)
  * [Physics Models](#physics-models)
  * [Detectors](#detectors)
  * [Output Data](#output-data)
  * [Notes for Developers](#notes-for-developers)
<!--te-->

## Getting the code
### Using git
The latest version of the code can be obtained by using the git command. This is possible after installing git on a linux machine (see https://git-scm.com/download/linux for further documentation on how to install git). Use the following command to download UNISim-tool:
````
$ git clone https://github.com/dellaquilamaster/UNISim-tool.git
````
### Downloading from Git Hub
The code can be downloaded also frm the Git Hub web page at the link: https://github.com/dellaquilamaster/UNISim-tool.git, by
clicking on the "Clone or Download" button on the right side of the page and then "Download ZIP". It is possible to donwload also a previous release of the code. For a complete list of all the releases please visit: https://github.com/dellaquilamaster/UNISim-tool/releases.
## Setup and Configuration
### Requirements
The code is compiled using the g++ compiler.
In order to compile and run the code ROOT >=6 is required (the program has been tested with version 6.16.00). It is not recommended to compile the code with an installation of ROOT <=5.
### Compile the Code
To compile the code and make a clean installation use the sequence of commands:
````
$ make clean
$ make -jN
````
N will be the number of parallel compilation processes that will be launched.
The binary file exec_UNISim-tool.exe is generated in the compilation.
### Configure the program
The program allows to perform simulations using different physics models and any arbitrary configuration of detectors. To configure the program modify the following file:
````
config/UNISim.conf
````
In the special language used to read the file, the character '\*' is used to provide a comment. Following a detailed list of the fields to configure:
* set VERBOSE_MODE : this option can be "true" or "false". If true, the program will run in verbose mode
* set GRAPHICAL_MODE : this option can be "true" or "false". If true, the program runs under the Event Visualization Environment of ROOT. WARNING: in this case it is not recommended to generate a large number of events
* set OUTPUT_DIRECTORY : this sets the directory where the output ROOT file will be stored, containing the output data stored into a TTree
* set RANDOM_SEED : this sets the seed use for the random number generator. Use "0" to have a unique sequence of random number at each time the program is launched
The physics reaction can be configured by the following fields:
* set BEAM : specify Z and A of the beam with the commands -Z=xx -A=yy
* set BEAM_ENERGY : beam kinetic energy in MeV
* set BEAM_POSITION : specify X, Y and Z coordinate of the position of the beam on the target (with respect to the target center) with the commands -X=xx -Y=yy -Z=zz
* set BEAM_ANGULAR_SPREAD : FWHM of the angular spread of the beam on the target in degrees
* set BEAM_POSITION_SPREAD : FWHM of the beam cross section on the target in cm
* set TARGET : specify Z and A of the target with the commands -Z=xx -A=yy
The physics model can be configured by using the following command:
* set PHYSICS_MODEL <the_model> "<the_physics_configuration_file>"
For a list of the available models and to learn how to prepare a physics configuration file please see section "Physics Models".
The properties of the target material can be configured as follows (these affect the energy of particles prior being detected as well as the energy of the beam before reacting, the reaction is assumed at mid-target):
* set TARGET_MATERIAL : symbol of the target material
* set TARGET_THICKNESS : target thickness in um
The detection setup can be configured as follows:
* define SETUP "<name>" : this initializes a detection setup called "name". The name of the setup will be the name of the output tree
* add <detector type> <options> : see section "Detectors" for a list of available detector types and options
### Limitations

## The HiRAEVTUnpacker Program
### Run the code
To run the UNISim-tool program use the following command from the main program directory:
````
$ ./exec_UNISim -events <N>
````
N will be the number of required events.
### Physics Models

### Detectors

### Output Data
Output data is stored in a tree called EXXXXX, where XXXXX represents the experiment name (i.e. E15190). The folder where the tree is stored is configured in the config file by setting HiRAEVTMAPPER_ROOT_FILE_PATH in the configuration file (see section "Configure the Program").

The structure of the output tree is constituted by an individual branch for each detector defined in the mapping file. In the presence of an "HiRA" detector, a sub-branch will be created for each telescope. For convenience, timestamp and TDC spare channels are treated like detectors. The second can map (see the map file contained in the folder "input") individual TDC channels that are used as benchmark, containing for example experiment triggers.
Here a summary of the data structures for each individual detector:

**_HiRA_**  
|- **fDE**  
|&nbsp;&nbsp;&nbsp;+ int fmulti  
|&nbsp;&nbsp;&nbsp;+ int fnumstrip\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ UShort_t fEnergyHi\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ UShort_t fEnergyLo\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ UShort_t fTime\[fmulti\]  
|- **fEF**  
|&nbsp;&nbsp;&nbsp;+ int fmulti  
|&nbsp;&nbsp;&nbsp;+ int fnumstrip\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ UShort_t fEnergyHi\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ UShort_t fEnergyLo\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ UShort_t fTime\[fmulti\]  
|- **fEB**  
|&nbsp;&nbsp;&nbsp;+ int fmulti  
|&nbsp;&nbsp;&nbsp;+ int fnumstrip\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ UShort_t fEnergyHi\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ UShort_t fEnergyLo\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ UShort_t fTime\[fmulti\]  
|- **fCsI**  
|&nbsp;&nbsp;&nbsp;+ int fmulti  
|&nbsp;&nbsp;&nbsp;+ int fnumcsi\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Short_t fEnergy\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Double_t fTime\[fmulti\]  
**_NeutronWall_**  
|&nbsp;&nbsp;&nbsp;+ int fmulti  
|&nbsp;&nbsp;&nbsp;+ int fnumbar\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Short_t fLeft\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Short_t fRight\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Short_t fFastLeft\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Short_t fFastRight\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Double_t fTimeLeft\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Double_t fTimeRight\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Double_t fGeoMean\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Double_t fFastGeoMean\[fmulti\]  
**_VetoWall_**  
|&nbsp;&nbsp;&nbsp;+ int fmulti  
|&nbsp;&nbsp;&nbsp;+ int fnumbar\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Short_t fBottom\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Short_t fTop\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Double_t fTimeBottom\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Double_t fTimeTop\[fmulti\]  
**_ForwardArray_**  
|&nbsp;&nbsp;&nbsp;+ int fmulti  
|&nbsp;&nbsp;&nbsp;+ int fnumdet\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Short_t fEnergy\[fmulti\]  
**_Microball_**  
|&nbsp;&nbsp;&nbsp;+ int fmulti  
|&nbsp;&nbsp;&nbsp;+ int fnumring\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ int fnumdet\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Short_t fTail\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Short_t fFast\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ Short_t fTime\[fmulti\]  
**_HTSisTimestampe15190_**  
|&nbsp;&nbsp;&nbsp;+ Long64_t fTimestamp  
|&nbsp;&nbsp;&nbsp;+ Long64_t fTimestampKoreans  
**_HTTDCSpare_**  
|&nbsp;&nbsp;&nbsp;+ Double_t (...)  
## Notes for Developers
