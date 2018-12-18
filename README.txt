The toolbox in this repository is designed to perform steady state optimizations of delayed systems. Constraints enforce robust asymptotic (or exponential) stability in the presence of uncertain parameters, uncertain delays and state dependent delays.
The development of this toolbox was funded by Deutsche Forschungsgemeinschaft (grant MO 1086/13).



HOW CAN I TEST THIS TOOLBOX VERY QUICKLY?

There are small examples that can be executed with minimal preparation effort. These are located in the folder matlab/delay/demo/. To test them, you need a Matlab installation including a compiler for c-code/mex-files.

Clone this repository. Start matlab and select "Home" in the graphical user interface. Then click on "Set path". In the new windows, select "Add Folder ..." and navigate to the matlab/delay folder of this repository. Navigate to matlab/delay/demo/, open the matlab-file "optimPopulationTestManifoldSlice.m" and run it with the play button in the matlab graphical user interface. The script contains sample code for a typical use case of the toolbox. 

The m-files in the subfolders of "ExampleApplications" contain additional examples. 


WHERE DO I FIND DOCUMENTATION?

The matlab code documentation is based on doxygen. It is located in the subfolder matlab/delay/Documentation. You can update your documentation using the respective makefile.



WHERE CAN I FIND MORE TECHNICAL BACKGROUND INFORMATION?

The following publications describe the technical details:

https://www.sciencedirect.com/science/article/pii/S2405896316306802
https://ieeexplore.ieee.org/document/7798469
