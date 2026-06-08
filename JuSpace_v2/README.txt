|===================================================================|
|                        WELCOME to JuSpace| v2.1                  |
|                                                                   |
|                                                                   |
|    Tool for spatial correlation analyses of magnetic resonance    |
|    		imaging data with positron emission                 |
|        tomography derived receptor and other cellularmaps         |
|                                                                   |
|===================================================================|
 
       			Created by Juergen Dukart
 
          PET data sources are listed in "Sources_templates_release.txt"

 
--------------------------------------------------------------------
 
INTRODUCTION:
JuSpace is a software package for the integration of different imaging modalities 
with positron emission tomography derived neurophysiological measures.
The toolbox is written for Matlab 2017b and following (The Mathworks Insc., MA, USA). 
It further requires SPM12 to be installed (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
Although it has been designed on a Windows platform, it can be used on any other OS 
like Unix or Linux where the Matlab package and the Stats toolbox are installed.

More information about generation of the neuromorphometrics atlas included in the 
toolbox can be found on the following pages:

 	General Segmentation: http://neuromorphometrics.com/Seg/
	BrainCOLOR: http://neuromorphometrics.com/ParcellationProtocol_2010-04-05.PDF

A list of regions provided in the Neuromorphometrics atlas is available using the following link:
   
	https://github.com/neurodebian/spm12/blob/master/tpm/labels_Neuromorphometrics.xml

The list of selected regions is also provided in the included "NeuromorphTemp_labels.csv" file.

For detailed Information including citations on the included PET and cellular maps please 
read the included "Sources_templates_release.txt" file.

GETTING STARTED:
Once the program is unpackaged in its directory, start Matlab. Use the  current directory  tab to go to the JuSpace 
directory (or add JuSpace to default path) and type
>>JuSpace

The user interface should appear.

Alternative the computing function "compute_DomainGauges" can be called directly, see help for this function on its usage.
To compute exact permutation based p-values for within and between-subject designs the function 
"compute_exact_pvalue" needs to be called. The function requires inputs provided by the "compute_DomainGauges" function


