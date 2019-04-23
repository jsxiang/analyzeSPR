# Analyze SPR binding affinity for data from Biacore X100

This depository contains the the MatLab scripts that perform binding affinity and kinetics analysis exported from a Biacore X100 instrument.

To perform the analysis for an example Biacore experiment that used the Single Cycle Kinetics protocol, run the script analyzeSCK.m. 
analyzeSCK.m uses input files FC2-1_SCK.csv as the sensorgram traces and the SCK_layout.xlsx detailing the run information.
FC2-1_SCK.csv is exported as a .txt file from the Biacore evaluation software, the first window that shows all the traces.
Microsoft Excel is used to convert the .txt file to .csv file for easy import into MatLab.
analyzeSCK.m requires other MatLab scripts: donlinmultifit.m, nlinmultifit.m, setfig.m.
analyzeSCK.m determines RNA capture, binding kinetics and affinity. 
analyzeSCK.m creates a .mat files containing a MatLab structure with the estimated parameters.

