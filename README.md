# Analyze SPR binding affinity for data from Biacore X100

This depository contains the the MatLab scripts that perform binding affinity and kinetics analysis exported from a Biacore X100 instrument.

To perform the analysis for an example Biacore experiment that used the Single Cycle Kinetics protocol, run the script analyzeSCK.m. 
analyzeSCK.m uses input files FC2-1_SCK.csv as the sensorgram traces and the SCK_layout.xlsx detailing the run information.
FC2-1_SCK.csv is exported as a .txt file from the Biacore evaluation software, the first window that shows all the traces.
Microsoft Excel is used to convert the .txt file to .csv file for easy import into MatLab.
analyzeSCK.m calls SCK.m to create an SCK object and parses the input data files, and performs analysis using function in SCK.m
analyzeSCK.m also requires other MatLab scripts: donlinmultifit.m, nlinmultifit.m, setfig.m.
analyzeSCK.m determines RNA capture, binding kinetics and affinity. 
analyzeSCK.m creates a .mat file containing a MatLab object with the estimated parameters.

