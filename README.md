# Analyze SPR binding affinity for data from Biacore X100

This depository contains the the MatLab scripts that perform binding affinity and kinetics analysis exported from a Biacore X100 instrument.

To perform the analysis for an example Biacore experiment that uses the Single Cycle Kinetics protocol, run the example script analyzeSCK.m. 
analyzeSCK.m uses input files e.g. FC2-1_SCK.csv as the sensorgram traces and an Excel file e.g. SCK_layout.xlsx detailing the run information (e.g. sample names and ligand concentrations etc.).
The .csv file with sensorgram traces is exported as a .txt file from the Biacore evaluation software, the first window that shows all the traces. Microsoft Excel is then used to convert the .txt file to .csv file for easy import into MatLab.
analyzeSCK.m calls SCK.m to create an SCK object and parses the input data files, and performs analysis using function in SCK.m
analyzeSCK.m also requires other MatLab scripts: donlinmultifit.m, nlinmultifit.m, setfig.m.
analyzeSCK.m determines RNA capture, binding kinetics and affinity. 
analyzeSCK.m creates a .mat file containing a MatLab SCK object with the estimated parameters.

To perform the analysis for an example Biacore experiment that uses the Multi-Cycle Kinetics protocol, run the example script analyzeMCK.m. 
analyzeMCK.m uses input files e.g. FC2-1_MCK.csv as the sensorgram traces and an Excel file e.g. MCK_layout.xlsx detailing the run information (e.g. sample names and ligand concentrations etc.).
The .csv file with sensorgram traces is exported as a .txt file from the Biacore evaluation software, the first window that shows all the traces. Microsoft Excel is then used to convert the .txt file to .csv file for easy import into MatLab.
analyzeMCK.m parses the input data files, and performs analysis using function in SCK.m
analyzeMCK.m also requires other MatLab scripts: donlinmultifit.m, nlinmultifit.m, setfig.m.
analyzeMCK.m determines RNA capture, binding kinetics and affinity. 
