|===================================================================|
|                        WELCOME to JuSpace|
|                                                                   |
|                                                                   |
|    Tool for spatial correlation analyses of magnetic resonance    |
|    		imaging data with positron emission                 |
|                 tomography derived receptor maps                  |
|                                                                   |
|===================================================================|
 
       			Created by Juergen Dukart
 
          PET data sources are listed in "Sources_templates_release.txt"

 
--------------------------------------------------------------------
 
INTRODUCTION:
JuSpace is a software package for the integration of different imaging modalities with positron emission tomography derived neurophysiological measures.
The toolbox is written for Matlab 2017b and following (The Mathworks Insc., MA, USA). It further requires SPM12 to be installed (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
Although it has been designed on a Windows platform, it can be used on any other OS like Unix or Linux where the Matlab package and the Stats toolbox are installed.

GETTING STARTED:
Once the program is unpackaged in its directory, start Matlab. Use the �current directory� tab to go to the JuSpace directory and type
>>JuSpace

The user interface should appear.

Alternative the computing function "compute_DomainGauges" can be called directly, see help for this function on its usage
To compute exact permutation based p-values for within and between-subject designs the function "compute_exact_pvalue" needs to be called. The function requires inputs provided by the "compute_DomainGauges" function