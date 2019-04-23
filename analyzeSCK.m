
sck=SCK(); % initialize object

numstartups=0; % note how many start up cycles were run (skipped in analysis)
timeintervals=[0 362 726 1086 1457]+400; % note when each association/dissocation starts

sck.load('FC2-1_SCK.csv','SCK_layout.xlsx',1,'A4:G24');
% % plot rawdata
% figure(1);clf
% plot(sck.rawdata(:,1:2:end),sck.rawdata(:,2:2:end),'o')

% determine and plot amt of RNA captured
sck.findCapture(numstartups,385:390,1:65)

% trim and normalize data
sck.prepBiadata(numstartups,timeintervals);

% perform kinetic fit to langmuir equation
beta0=[1e4 3e-3 3e1]; % initial guesses
datarange=1:5; % choose the number of binding events that will be useful for fit
showfit=true;
sck.fitKineticBiadata(beta0,datarange,showfit)

% perform dose response fit 
ic50guesses=[.000001 0.1  0.1 10];
datarange=1:5; % choose the number of binding events that will be useful for fit
showfit=true;
sck.fitEquilibriumBiadata(ic50guesses,datarange,showfit);

% display results of estimated parameters
sck.dispFit();

% save object
save('SRet_sck.mat','sck');



