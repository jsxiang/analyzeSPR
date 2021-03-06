clear
warning('OFF')
% specify options
numstartupcycles=1;

capturebaseinterval=56:61; % time interval for determining baseline RU
captureinterval=270:275; % time interval for measuring RNA capture

timestart=338;
contacta=180; % association time, in seconds

startcapturesteadystate=260;

doskipspikes=true; % trim injection spikes
doflankingbeforeandafter=true; % normalize signal to regions flanking binding event
dozerobase=true; % zero baseline
dorefsubtract=true; % subtract traces by a negative control

%% load sensorgrams
data=csvread('FC2-1_MCK.csv',1,0);

%% load experimental information
samplenames={};
ligandname={};
ligandconc=[];
cycnum=[];
cycparam={};
xind=(numstartupcycles*2+1):2:length(data(1,:));
yind=(numstartupcycles*2+2):2:length(data(1,:));
ligandconcunit='\muM';


allligandconcunit=ligandconcunit;

[n,t,r]=xlsread('MCK_layout.xlsx',1,'A1:D26');
allnames=r(:,1);
ligconc=r(:,4);
allligandconc=[];
allligandname={};
for i=1:length(ligconc)
   allligandconc(i)= ligconc{i};
    if ~strcmp(r(i,2),'Buffer')
        ligandname=r(i,2);
    allligandname{end+1}=ligandname{1};
    else
            allligandname{end+1}='Buffer';
    end
end


%% find capture/baseline
setfig('capture');clf
capture=mean(data(captureinterval,yind))-mean(data(capturebaseinterval,yind));
plot(1:length(yind),capture,'x','linewidth',2,'Markersize',12);
% ylim([0,1.2*max(capture)])
xlabel('cycle no.')
ylabel('Capture (RU)')
grid on
set(gca,'fontsize',16)
set(gca,'linewidth',1.5)
saveas(gcf,'capture','png')

%% baseline subtracted data
basesubdata=zeros(size(data(:,yind)));
basesubdata(:,yind)=data(:,yind)-repmat(capture,length(data(:,1)),1);

%% time zeroed data
basesubdata(:,xind)=data(:,xind)-timestart;

allxdata=basesubdata(:,xind);
allydata=basesubdata(:,yind);
%% Prep data
uniquesamplenames=unique(allnames);
biadata=struct;

for i=1:length(uniquesamplenames)
    bufferi=[];
    buffername={};
    bufferdata=[];
    for j=1:length(allnames)
        if abs(allligandconc(j)-0)<1e-8 & regexp(allnames{j},uniquesamplenames{i})
            bufferdata=[bufferdata allydata(1:end-1,j)];
        end
    end

    biadata(i).name=uniquesamplenames{i};
    biadata(i).meanbufferdata=mean(bufferdata,2);
    biadata(i).bufferdata=bufferdata;
end



for i=1:length(biadata)
    legendnames={};
    setfig(uniquesamplenames{i});clf
    hold on
    biadata(i).tdata=[];
    biadata(i).ydata=[];
    biadata(i).skipspikes=[];
    biadata(i).ligandconc=[];
    for j=1:length(allnames)
       if abs(allligandconc(j)-0)<1e-8
           continue;
       elseif regexp(allnames{j},strcat('^',uniquesamplenames{i},'$'))
            fprintf('i=%d,j=%d',i,j)
            tc=find(abs(allxdata(:,j)+timestart-startcapturesteadystate)==min(abs(allxdata(:,j)+timestart-startcapturesteadystate)));
            t0=find(abs(allxdata(:,j)-0)==min(abs(allxdata(:,j)-0)));
            td=find(abs(allxdata(:,j)-contacta*2)==min(abs(allxdata(:,j)-contacta*2)));
            
            if doflankingbeforeandafter
                flankingtimesteps=[(tc):t0-1 (td+100):(td+250)];
            else
            flankingtimesteps=[(tc):t0-1];
            end
            
            ta=allxdata(flankingtimesteps,j); % flanking
%             ya=allydata((260):t0-1,j); % association
            ya=allydata(flankingtimesteps,j)-biadata(i).meanbufferdata(flankingtimesteps);
            ynorm=(ya-min(ya))/(max(ya)-min(ya));

            t=allxdata((tc):td,j); % all time
            ydata=allydata((tc):td,j)-biadata(i).meanbufferdata((tc):td);
            
            if dozerobase
                % first fit exponential decay curve

                F= @(x,xdata) x(1)+x(2).*(exp(-x(3).*xdata/60)); % single exponential decay
                x0=[ynorm(1) ynorm(end) 0.5];
                [x, resnorm,~,exitflag,output]=lsqcurvefit(F,x0,ta,ynorm);

                if exitflag==1
                    ynormfit=F(x,t).*(max(ya)-min(ya))+(min(ya));
                    y=ydata-ynormfit;
                else
                    p=polyfit(ta,ynorm,1);
                    yfit=polyval(p,t);
                    ynormfit=yfit.*(max(ya)-min(ya))+(min(ya));
                    y=ydata-ynormfit;
                end
            else
                y=ydata;
            end
            skipspikea=abs(t-0)<2;
            skipspiked=abs(t-contacta)<3;
            skipspikes=skipspikea|skipspiked;
            if doskipspikes
                plot(t(~skipspikes), y(~skipspikes),'-','linewidth',2)
            else
                plot(t,y,'-','linewidth',2)
            end
            biadata(i).tdata=[biadata(i).tdata t];
            biadata(i).ydata=[biadata(i).ydata y];
            biadata(i).ligandconc(end+1)=allligandconc(j);
            biadata(i).skipspikes=[biadata(i).skipspikes skipspikes];
            legendnames{end+1}=sprintf('%0.3f %s %s ',allligandconc(j),allligandconcunit,allligandname{j});       
                
       end
    end
    
    legend(legendnames,'location','northeast')
    xlabel('time (s)')
    ylabel('Blank subtracted SPR Response Unit')
    title(uniquesamplenames{i},'interpreter','none')
    set(gca,'fontsize',16)
    set(gca,'linewidth',1.5)
    saveas(gcf,sprintf('%s_SPR',uniquesamplenames{i}),'png')
    
end


%% subtract N*N* from sample curves and perform binding affinity and kinetics analyses
refi=[];
for i=1:length(biadata)
    if ~isempty(regexp(biadata(i).name,'N[0-9]+N[0-9]+')) && isempty(regexp(biadata(i).name,'bad'))
    refi(end+1)=i;
    end
    
end

ref=refi(1);
showfit=true;
ligandname='trans-Zeatin';
refdatarange=1:5;
datarange=1:5;
ic50guesses=[.1 0.1  0.1 10 ];
earlytimepoint=true;
tssplacement=5;


runi=1:length(biadata);
runi(runi==ref)=[];


fid=fopen('SPRanalyzed.csv','w');
fprintf('\tName, \tka, \tkd, \tKD_kinetic, \tKD_equilibrium\n');

for i=runi

datalen=min(length(biadata(i).tdata),length(biadata(ref).tdata));
% biadata(i).ydata_refsubtract=biadata(i).ydata(1:datalen,:)-biadata(ref).ydata(1:datalen,:);
biadata(i).tdata_refsubtract=biadata(i).tdata(1:datalen,datarange);
biadata(i).ydata_refsubtract=biadata(i).ydata(1:datalen,datarange)-biadata(ref).ydata(1:datalen,refdatarange);
C=[];
ta=[];
ya=[];
td=[];
yd=[];

legendnames={};
for j=datarange-min(datarange)+1
tss0=find(abs(biadata(i).tdata_refsubtract(:,j)-(1))==min(abs(biadata(i).tdata_refsubtract(:,j)-(1))));
tssend=find(abs(biadata(i).tdata_refsubtract(:,j)-(contacta-1))==min(abs(biadata(i).tdata_refsubtract(:,j)-(contacta-1))));


% ta=[ta biadata(i).tdata_refsubtract([tss0-3 tss0:tssend],j)];
% ya=[ya biadata(i).ydata_refsubtract([tss0-4 tss0:tssend],j)];

ta=[ta biadata(i).tdata_refsubtract([tss0-1 tss0:tssend],j)];
ya=[ya biadata(i).ydata_refsubtract([tss0-1 tss0:tssend],j)];

td0=find(abs(biadata(i).tdata_refsubtract(:,j)-(contacta+6))==min(abs(biadata(i).tdata_refsubtract(:,j)-(contacta+6))));
tdend=find(abs(biadata(i).tdata_refsubtract(:,j)-(contacta*2))==min(abs(biadata(i).tdata_refsubtract(:,j)-(contacta*2))));
tdspan=td0:tdend;

td=[td biadata(i).tdata_refsubtract([td0-1 tdspan],j)];
yd=[yd biadata(i).ydata_refsubtract([tssend tdspan],j)];



C=[C biadata(i).ligandconc(j)/1e6];
beta0=[1e3 3e-3 30]; % initial guesses for global fits only

legendnames{end+1}=sprintf('%0.3f %s %s',biadata(i).ligandconc(j),allligandconcunit,ligandname);     

end


name=sprintf('%s - %s',biadata(i).name,biadata(ref).name);
setfig(name);clf
title(name,'interpreter','none')

plot([ta;td],[ya;yd],'.','MarkerSize',12)
hold on

[beta,r,J,Sigma,mse,errorparam,robustw]=donlinmultifit(C,ta,ya,td,yd,beta0,showfit);
legend(legendnames,'location','best')
ylim([-7.5,40])


% determine steady state affinity
tss=find(abs(allxdata(:,j)-(contacta-3))==min(abs(allxdata(:,j)-(contacta-3))));

tssinterval=(tss-5):tss;
tss=find(abs(biadata(i).tdata_refsubtract(:,1)-(contacta-3))==min(abs(biadata(i).tdata_refsubtract(:,1)-(contacta-3))));
if earlytimepoint
    tss=find(abs(biadata(i).tdata_refsubtract(:,1)-(contacta-tssplacement))==min(abs(biadata(i).tdata_refsubtract(:,1)-(contacta-tssplacement))));
end
    tssinterval=(tss-5):tss;


yss=mean(biadata(i).ydata_refsubtract(tssinterval,:));
yssnormalized=(yss-min(yss))/(max(yss)-min(yss));


ic50timepoint=contacta-tssplacement;
yL=get(gca,'YLim');
hold on
if earlytimepoint
plot(ones(1,100)*ic50timepoint,linspace(yL(1),yL(2)),':','linewidth',2)
else
plot(ones(1,100)*(contacta-3),linspace(yL(1),yL(2)),':','linewidth',2)
end

setfig(strcat(biadata(i).name,' IC50'));clf
hold on

[hillCoeff ec50]=doseResponse(biadata(i).ligandconc(datarange),yssnormalized,ic50guesses);
ylabel('normalized RU')
xlabel(sprintf('ligand concentration (%s)',allligandconcunit))
title(strcat(biadata(i).name,' IC50'),'interpreter','none')
set(gca,'linewidth',1.5)
set(gca,'XScale','log')

hold on
plot(biadata(i).ligandconc(datarange),yssnormalized,'o','linewidth',2)
text(0.9*(mean(log10(biadata(i).ligandconc))),0.55,sprintf('KD_e_q = %0.2f \\muM',ec50),'fontsize',14)
saveas(gcf,strcat(biadata(i).name,'_IC50'),'png')

biadata(i).KDss=real(ec50);
biadata(i).eqRU=yssnormalized;

betamat(i).beta=beta;
betamat(i).ec50=ec50;
fprintf('%s, ',biadata(i).name);
fprintf('%2.4e, %2.4e, %2.4e, ',betamat(i).beta(1),betamat(i).beta(2),betamat(i).beta(2)/betamat(i).beta(1));
fprintf('%2.4e',betamat(i).ec50*1e-6);
fprintf('\n');


fprintf(fid,'%s, ',biadata(i).name);
fprintf(fid,'%2.4e,%2.4e,%2.4e, ',betamat(i).beta(1),betamat(i).beta(2),betamat(i).beta(2)/betamat(i).beta(1));
fprintf(fid,'%2.4e \n',betamat(i).ec50*1e-6);


end
fclose(fid);














%% END







