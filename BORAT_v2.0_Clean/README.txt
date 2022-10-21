Borat 2.0

The Eikonal Solver is added to this version to reduce computational recources.  

This version of the code runs with axis-symmetrical and bidimensional solutions. 

In the files ARD.Case, ExoMars_2D.Case and ExoMars_Larsen.Case you can find the simulation
paramenters for each different case. You can copy paste the content of these files in the first section of the BORAT_main.m file. 

Run with the flag collisionmodel=0 to deactivate the use of Mutation++. 

Run with the flag flag_exportdomain=0 to deactivate the use of import-export between Matlab and Tecplot.


During the run, a folder named RunOutput is created and all the figures and solution files are saved in this folder.

The simulation solution file (in Matlab format .mat) is saved at the end of the simulation to get access to the solution data. 

All the main functions of the code are saved in their corresponding function definition files. You can easily test each of them by loading the solution.mat file of the simulation, and call each of them separately (you can use the file test.m for this)

Figure axis are set on size of domain of bidimensional solutions, change them in case of axisymmetric solution (reduce axis limits).


