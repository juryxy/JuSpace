JuSpace Executable

1. Prerequisites for Deployment 

Verify that MATLAB Runtime(R2023a) is installed.   
If not, you can run the MATLAB Runtime installer.
To find its location, enter
  
    >>mcrinstaller
      
at the MATLAB prompt.

Alternatively, download and install the Linux version of the MATLAB Runtime for R2023a 
from the following link on the MathWorks website:

    https://www.mathworks.com/products/compiler/mcr/index.html
   
For more information about the MATLAB Runtime and the MATLAB Runtime installer, see 
"Distribute Applications" in the MATLAB Compiler documentation  
in the MathWorks Documentation Center.

2. Files to Deploy and Package

Files to Package for Standalone 
================================
-JuSpace 
-run_JuSpace.sh (shell script for temporarily setting environment variables and executing 
                 the application)
   -to run the shell script, type
   
       ./run_JuSpace.sh <mcr_directory> <argument_list>
       
    at Linux or Mac command prompt. <mcr_directory> is the directory 
    where MATLAB Runtime(R2023a) is installed or the directory where 
    MATLAB is installed on the machine. <argument_list> is all the 
    arguments you want to pass to your application. For example, 

    If you have MATLAB Runtime(R2023a) installed in 
    /mathworks/home/application/R2023a, run the shell script as:
    
       ./run_JuSpace.sh /mathworks/home/application/R2023a
       
    If you have MATLAB installed in /mathworks/devel/application/matlab, 
    run the shell script as:
    
       ./run_JuSpace.sh /mathworks/devel/application/matlab
-MCRInstaller.zip
    Note: if end users are unable to download the MATLAB Runtime using the
    instructions in the previous section, include it when building your 
    component by clicking the "Runtime included in package" link in the
    Deployment Tool.
-This readme file 



3. Definitions

For information on deployment terminology, go to
https://www.mathworks.com/help and select MATLAB Compiler >
Getting Started > About Application Deployment >
Deployment Product Terms in the MathWorks Documentation
Center.

4. Appendix 

A. Linux systems:
In the following directions, replace MR/R2023a by the directory on the target machine 
   where MATLAB is installed, or MR by the directory where the MATLAB Runtime is 
   installed.

(1) Set the environment variable XAPPLRESDIR to this value:

MR/R2023a/X11/app-defaults


(2) If the environment variable LD_LIBRARY_PATH is undefined, set it to the following:

MR/R2023a/runtime/glnxa64:MR/R2023a/bin/glnxa64:MR/R2023a/sys/os/glnxa64:MR/R2023a/sys/opengl/lib/glnxa64

If it is defined, set it to the following:

${LD_LIBRARY_PATH}:MR/R2023a/runtime/glnxa64:MR/R2023a/bin/glnxa64:MR/R2023a/sys/os/glnxa64:MR/R2023a/sys/opengl/lib/glnxa64

    For more detailed information about setting the MATLAB Runtime paths, see Package and 
   Distribute in the MATLAB Compiler documentation in the MathWorks Documentation Center.


     
        NOTE: To make these changes persistent after logout on Linux 
              or Mac machines, modify the .cshrc file to include this  
              setenv command.
        NOTE: The environment variable syntax utilizes forward 
              slashes (/), delimited by colons (:).  
        NOTE: When deploying standalone applications, you can
              run the shell script file run_JuSpace.sh 
              instead of setting environment variables. See 
              section 2 "Files to Deploy and Package".    






