# UNISim-tool
A Unified Simulation tool for the Rudjer Boskovic Institute nuclear clustering group.

Contributors: Daniele Dell'Aquila

Rudjer Boskovic Institute, Zagreb, Croatia

This software is a unified simulation tool for the nuclear clustering group of the Rudjer Boskovic Institute. It is able to reliably run Monte Carlo simulations involving different physics models (such as sequential decay reactions) and an arbitrary configuration of detectors (DSSSDs, Lamp Wedge detectors). This software contains a high-precision energy loss program extracted by LISE++ for the calculation of the energy loss of charged particles in the matter. Another interesting novelty of this software is the advanced graphical mode for the visualization of the reaction and the detection setup, this is based on the Event Visualization Environment of ROOT.

The program supports the use in parallel mode on the "Malanar" compting cluster at the RBI.

# Useful links:
  * [Rudjer Boskovic Institute](https://www.irb.hr/eng) : Rudjer Boskovic Institute official website
  
Table of contents
=================
<!--ts-->
* [Getting the code](#getting-the-code)
  * [Using git](#using-git)
  * [Downloading from Git Hub](#downloading-from-git-hub)
* [Setup and Configuration](#setup-and-configuration)
  * [Requirements](#requirements)
  * [Environment Variables](#environment-variables)    
  * [Compile the Code](#compile-the-code)
  * [Configure the Program](#configure-the-program)
  * [Limitations](#limitations)
* [The UNISim-tool Environment](#the-unisim-tool-environment)
  * [Run the Code](#run-the-code)
  * [Run the Code in parallel on Malnar](#run-the-code-in-parallel-on-malnar)  
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
### Environment Variables
Before compiling or running the code, it is required to export the needed environment variables and aliases. This can be done by launching the command:
````
$ source path-to-framework/UNISim-tool.sh
````
This can usually be done automatically by adding the previous line to the bashrc file.
### Compile the Code
To compile the code and make a clean installation use the sequence of commands:
````
$ make distclean
$ make -jN
$ make install
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
* define SETUP "<setup_name>" : this initializes a detection setup called "name". The name of the setup will be the name of the output tree
* add <detector type> <options> : see section "Detectors" for a list of available detector types and options
### Limitations

## The The UNISim-tool Environment
### Run the code
To run the UNISim-tool program use the following command from the main program directory:
````
$ ./exec_UNISim.exe [OPTIONS]
````
[OPTIONS]  
* -events: this can specify the number of events. This option is mandatory;  
* -o: this can specify the name of the output file;  
* -physics: specify the physics model;  
* -reaction: specify the pathname to the reaction file.  
### Run the code in parallel on Malnar
The code can be run in parallel by using the Malanar computing cluster of the RBI. When the code runs in parallel mode, the number of events to generate are divided into subgroups, and a series of individual independent processes is launched in 
batch mode. In this case, it is important that the user sets RANDOM_SEED to 0 (see Sect "Configure the program"), in order to obtain a different sequence of random numbers for each individual process.  
To run the code in parallel mode on Malnar, the user can use the script:
````
run_UNISim-tool_Malnar.sh
````
It must be marked as executable:
````
chmod +x run_UNISim-tool_Malnar.sh
````
The variable "MAXEVENTS" sets the number of events processed by each individual process (it is not recommended to change this value for non-expert users).  
To run the script use the following command:
````
./run_UNISim-tool_Malnar.sh <N>
````
Where N will be the total number of events to process. WARNING: Make sure the folders "batch_output" and "batch_errors" exist. 
### Physics Models
The program includes a number of already implemented physics models. The user can add new physics models very easily. When the user adopts a physics model, he must specify a configuration file of the specific reaction process. Reaction configuration files are contained in the folder:
````
reactions/
````
The available physics models are:
* SequentialDecay : used to simulate reactions with the sequential decay of one or more excited nuclei.  
* RutherfordScattering : used to simulate elastic scattering processes according to the Rutherford cross-section.  
Additional physics models will be implemented in the near future.  
Following, a detailed explanation of each available physics model:  
* SequentialDecay : This model can simulate sequential reactions where one or more of the products of the reaction can decay in additional decay products. One of the additional decay product can decay in other decay products and so on. The novelty of this model is that it can simulate ad arbitrary number of steps with an arbitrary number of particles and excited nuclei, each of those with an arbitrary decay pattern and an arbitrary spectroscopy.   
Usage example:  
9Li + 9Be -> 16C* + d -> 6He + 10Be + d  
We first set the particles of the primary collision: 16C* and d  
````
define particle 0 -Z=6 -A=16 : 16C excited defined as particle "0". This will decay according to the decay scheme below
define particle 1 -Z=1 -A=2 : light recoil
````
In this case 16C is excited. We define the spectroscopy of 16C  
````
set spectroscopy 0 -Ex=21 -Gamma=0.150
set spectroscopy 0 -Ex=28 -Gamma=0.300
````
I have 2 excited states one at 21 MeV and one at 28 MeV with Gamma = 150 keV and 300 keV respectively.  
If the nucleus decays I need also to define a decay pattern (e.g. 6He+10Be):  
````
set decay 0 0 -Z=2 -A=6 *** particle 0 of the decay of the fragment 0
set decay 0 1 -Z=4 -A=10 *** particle 1 of the decay of the fragment 0
````
In this way I have defined the particles 0_0 and 0_1. The exercise can continue in the same way by setting a spectroscopy to those particles and a decay pattern and so on. In a more general case one can have the following reaction:  
P=Projectile  
T=Target  
x1,x2,...,xn = light ejectiles of the primary reaction  
X* heavy residual excited that decays into Y* + y1 + ... + ym  
P+T -> X* + x1 + x2 + ... + xn -> Y* + y1 + ... +ym + x1 + ... + xn  
Equivalently, Y* can decay into other products and so on... For each step of the reaction, one might have even more than 1 product decaying, e.g. X1*, ..., Xl* in the first step of the reaction and a number of Yi* in the second and so on.
* RutherfordScattering : This model can simulate elastic scattering events according to a pure Rutherford cross-section.
Usage example:  
20Ne + 197Au -> 20Ne + 197Au  
We set exclusively the particle of the primary collision: 20Ne and 197Au  
````
define particle 0 -Z=10 -A=20 : scatteted beam (ejectile)
define particle 1 -Z=79 -A=197 : recoiling target (recoil)
````
Additionally, one can set the minimum and/or maximum polar angles for the random generator in degrees:  
````
set min_angle X
set max_angle X
````  
By default, data will be simulated by using the standard Rutherford cross-section. However, the user can also define a "user-defined" angular distribution by using the command:
````
set user_defined_distribution formula_with_._instead_of_* num_parameters par0 par1 ...
````  
In the formula, each parameter is indicated with the sintax [N] (where N is a number from N to NumParameters-1 and all the "*" are replaced by ".".  
### Detectors
In the framework of UNISim-tool, the user can define an arbitrary number of detectors with arbitrary configurations. A list of available detectors and the corresponding options follows:  
* DSSSD (Double-Sided Silicon Strip Detector facing the target perpendicularly) -> options are: -distance (perpendicular distance from the target), -theta (polar angle of the center in degrees), -phi (azimuthal angle of the center in degrees), -strips (number of front or back strips), -strip_width (width of one strip in cm), -inter_strip (width of the interstrip in cm), -frame_width (width of the ceramic frame in cm), -dead_layer (size of a dead region of silicon before the frame in cm)
* DSSSD_ROT (Double-Sided Silicon Strip Detector with a customized position and rotation) -> options are: -X0 (X-position of the center), -Y0 (Y-position of the center), -Z0 (Z-position of the center), -tilt_X (tilt angle with respect to the X-axis in degrees), -strips (number of front or back strips), -strip_width (width of one strip in cm), -inter_strip (width of the interstrip in cm), -frame_width (width of the ceramic frame in cm), -dead_layer (size of a dead region of silicon before the frame in cm)
* LAMP_WEDGE (Lamp detector wedge) -> options: -distance (distance from the inner part of the wedge to the target) -phi_pos (azimuthal position in degrees), -tilt (with respect to the orizontal axis), -frame_distance (distance from the bottom of the frame to the beam axis), -strips (number of strips), -strip_width (width of one strip in cm), -inter_strip (width of the interstrip in cm)
* FAZIA_QUARTET (A Fazia quartet) -> options: -displacement (Z-displacement from the nominal 100 cm position) -theta_pos (theta of the quartet center when placed at the nominal 100 cm position), -phi_pos (phi of the quartet center when placed at the nominal 100 cm position), -pad_width (width of each pad in the quartet), -frame_width (width of pad frame)
### Output Data
Output data is stored in a tree called as the experimental setup. The folder where the tree is stored is configured in the config file by setting OUTPUT_DIRECTORY in the configuration file (see section "Configure the Program").  
The structure of the output tree is constituted by an individual branch containing the data stucture UNISRootEvent. This is composed by a number of sub-branches listed below:  
**_UNISRootEvent_**  
|&nbsp;&nbsp;&nbsp;+ int fmulti  
|&nbsp;&nbsp;&nbsp;+ bool fIsDetected\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ int fnumdet\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ int fnumpixel\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ double fKinEnergy\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ double fXDetHit\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ double fYDetHit\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ double fZDetHit\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ double fKinEnergyAfterTarget\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ double fKinEnergyOrigin\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ double fThetaOrigin\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ double fPhiOrigin\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ double fKinEnergyOriginCms\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ double fThetaOriginCms\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ int fZ\[fmulti\]  
|&nbsp;&nbsp;&nbsp;+ int fA\[fmulti\]  
## Notes for Developers
