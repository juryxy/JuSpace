# JuSpace
This is the repository for the JuSpace toolbox allowing for cross-modal correlations between imaging data and positron emission tomography derived information
|===================================================================|
|                        WELCOME to JuSpace v 1.3                   |
|                                                                   |
|                                                                   |
|    Tool for spatial correlation analyses of magnetic resonance    |
|    		imaging data with positron emission                          |
|                 tomography derived receptor maps                  |
|                                                                   |
|===================================================================|
 
       			Created by Juergen Dukart
          
          Please contact me at juergen.dukart@gmail.com for reports on bugs or other questions with respect to the toolbox
 
          PET and SPECT data sources are listed in "Sources_templates_release.txt". Please cite the sources corresponding to respective PET or SPECT maps

 
--------------------------------------------------------------------
 
INTRODUCTION:
JuSpace is a software package for the integration of different imaging modalities with positron emission tomography derived neurophysiological measures.
The toolbox is written for Matlab 2017b and following (The Mathworks Insc., MA, USA). It further requires SPM12 to be installed (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
Although it has been designed on a Windows platform, it can be used on any other OS like Unix or Linux where the Matlab package and the Stats toolbox are installed.

GETTING STARTED:
Once the program is unpackaged in its directory, start Matlab. Use the ‘current directory’ tab to go to the JuSpace directory and type
>>JuSpace

The user interface should appear.

Alternative the computing function "compute_DomainGauges" can be called directly, see help for this function on its usage
To compute exact permutation based p-values for within and between-subject designs the function "compute_exact_pvalue" needs to be called. The function requires inputs provided by the "compute_DomainGauges" function

More information about generation of the neuromorphometrics atlas included in the toolbox can be found on the following pages:

 	General Segmentation: http://neuromorphometrics.com/Seg/
	
	BrainCOLOR: http://neuromorphometrics.com/ParcellationProtocol_2010-04-05.PDF

A list of anatomical labels for the for the included neuromorphometrics (asymetric - 119 regions version) atlas is available using the following link:
   
	https://github.com/neurodebian/spm12/blob/master/tpm/labels_Neuromorphometrics.xml
 
Version history:

25.05.2020 
- fixed some visualization color mismatch for linux, 
- fixed a setting for computational option 4 (list 1 each image) and option 8 (allowing to test correlation distribution for all files from list 1 against null)

25.02.2021
- fixed a bug in computation of an exact p-value (leaded to error message under certain conditions)

30.07.2021 - Version 1.2 release
- A major update adding exact spatial permutations statistics for options 3,4 and 8 (highly recommended to use instead of the parametric p-value). Random permutations are performed by creating randomly permuted PET maps and reintroducing spatial auto-correlation (if present in the original). Auto-correlation is reintroduced by smoothing the permuted data (similar to the concept introduced by Burt et al. 2018, Nature Neuroscience, 21(9):1251–1259).

01.12.2021 - Version 1.3 release
- Fixed a bug in visualization of bar plots (shifted x labels when too many PET maps in one plot)
- Added many new PET maps

03.02.2022
- Fixed a bug in loading the maps when files with different dimensions were loaded into the same list.

18.02.2022
- Fixed a bug in generate_spatial_nullMaps (while loop did not converge for some of the newly added PET maps)
