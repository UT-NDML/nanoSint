# nanoSint 
This repository contains scripts related to modelling the diffusion process between nanoparticles. The diffusion between particles is tracked using a Phase Field Modelling (PFM) approach. These simulations use density and crystal order parameters as the phase field variables, tracking the diffusion process through the temporal evolution of the phase field variables. Evolution of the density variable is governed with the Cahn-Hilliard equation and the order parameter with the Landau-Gizburg equation. 

## Code breakdown
The scripts used for these simulations are divided into three major categories, the pre-processing bedGeneration scripts, the sintering simulation scripts, and the post-processing imaging and analysis scripts.

### Bed Generation Scripts
The bed generation scripts 

### Sintering Simulation Scripts
Parallel scripts run on TACC and have been tested with impi on TACC and openmpi on a regular desktop computer. 

Compilation

_**mpicxx -O3 scriptName.cpp**_

Running
<on TACC with 64 processors>

_**ibrun -np 64 ./a.out noP height**_

noP is the number of particles in the simulation bed and height is the height of the simulation box. For a single layer simulation the scriptName is singleLayer.cpp and there must be a bed.txt file in the same directory which contains the information for particles in the bed in the format of:

radius  xCenter   yCenter   zCenter

for each particle in the bed. singleLayer.cpp script assumes that all particles in the bed are being heated at the same temperature (an isothermally heated bed). 

For multiple layers multiLayer.cpp simulates the sintering of more than one layer. This requires a newLayer.dat file with the information for all the particles in the new bed. This .dat file is created with an init script in the bedGeneration folder. For this script the assumption is also made that all particles are being heated isothermally within the bed. The isothermal heating assumption is relaxed in multiLayerTemperature.cpp. This script as written performs gradient heating for a multi layer bed and so requires a newLayer.dat file as well as a tempProfile.dat file. The temProfile.dat file contains the temperature information corresponding to each pixel. The temperatures in this current version are 0, 450, 500 and 550. The temperature map for these temperatures to the sintering parameters are defined in a temperature map in the script.

On TACC the simulations are run with a batchScript (indiBatch) which defines the queue to run the simulation in, as well as the numper of processors and nodes to use in the simulation. The running procedure:

_**sbatch indiBatch**_

### Post Processing Scripts
Imaging is done in parallel with the plot_2by2.py and in series with import_auto_fromdat.py the size of the bed is set in the script with the xsize, ysize and zsize parameters. The naming convention used for reading into the file is also set in the script. Analysis is done with the limanalys.py script (in parallel) or analysis_fromdat.py (single core). In TACC procdevBatch used to run the python scripts.

* plot_2by2 
  * command line inputs:
    * zsize, ysize, xsize --> simulation box dimensions
    * SNo --> serial number tag used in naming convention for files
    * minfile, maxfile --> start and end file numbers for the span of files to plot
    * filed --> step distance between files (default 50)
    * sizep --> number of processors to use in analysis
    * simzoom --> sets whether to plot whole bed (1), just center box (2) or both (3)
  * Output is the png image files
  
* limanalys 
  * command line inputs are the same as for the plots above with the exception of the simzoom input
  * additional variables used to set parallelization. List includes directories to navigate to as well as the corresponding box dimensions for the beds in each directory. Also includes the analysis box dimensions defined in the xvec and yvec lists
  * Output are text files with the results of the analysis. The output results are then sent into the calibration scripts
  
* TACCtimecalibplot_R2error_withFitting
  * performs the calibration using the analysis results. The calibration is done within fitting of temperature bands defined in the AllCoeficients files. These files should be in the same directory as python script when run.
  * In file settings *bedList* defines the list of directories containing the analysis files. Modifications should correspond to the SNo variable that sets the file naming. xbd and ybd also affect the file naming. These correspond to the file analysis boxes
