# 2armbanditvariant
data analysis and plotting of the data from the behavioral experiments.

some of these functions/scripts are for analyzing the sessions after all sessions are done and are specifically for the 2Armbanditvariant with the Risk type FixedwithHolding.
the main function lets you choose a directory from where to take the session files for the analysis and also whether session-per-session should be analyzed or all sessions pooled together (not yet implemented)
the code for the data analysis and the plotting is in different scripts. In order to run the code, see file: main

the Analysis.mat file is a function to implement in the bpod protocol of the 2armbanditvariant task.
it has a switch statement to to analysis specifically for the different Risktypes defined in the Gui
