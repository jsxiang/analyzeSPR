classdef SCK < handle

    properties
        input;
        rawdata;
        capture;
        biadata;
    end

    methods
        function obj=SCK()
        end

        function load(data,file,cyclesfile,sheetnum,colrow)
            if isempty(file)
                file='FC2-1.csv';
            end
            d=csvread(file,1,0);
            fprintf('Loading data from %s...\n',file);

            samplenames={};
            ligandname={};
            ligandconc=[];
            cycnum=[];
            cycparam={};
            ligandconcunit='\muM';
            

            [num,txt,raw] =xlsread(cyclesfile,sheetnum,colrow);

            samplecycnum=[];
            
            for i=1:length(num(:,1))
                samplenames{end+1}=txt{i,1};
                ligandconc=[ligandconc; num(i,:)];
                ligandname{end+1}=txt{i,2};
                
            end
            
            data.input.samplenames=samplenames;
            data.input.ligandconc=ligandconc;
            data.input.ligandname=ligandname;
            data.input.ligandconcunit=ligandconcunit;
            
            data.rawdata=d;
            
        end
        

        function findCapture(data,numstartupcycles,captureinterval,capturebaseinterval)
            % find capture/baseline
            d=data.rawdata;
            xind=(numstartupcycles*2+1):2:length(d(1,:));
            yind=(numstartupcycles*2+2):2:length(d(1,:));
            data.capture=mean(d(captureinterval,yind))-mean(d(capturebaseinterval,yind));
            
            setfig('capture');clf
            plot(1:length(yind),data.capture,'x','linewidth',2,'Markersize',12);
            % ylim([0,1.2*max(capture)])
            xlabel('cycle no.')
            ylabel('Capture (RU)')
            grid on
            set(gca,'fontsize',16)
            set(gca,'linewidth',1.5)
            saveas(gcf,'capture','png')      

        end
        
        
        function prepBiadata(data,numstartupcycles,timeintervals)
            % parse out unique RNA species run
            sample_ligand_cycles=[];
            sample_buffer_cycles=[];
            sample_ligand_names={};
            
            data.biadata=struct;
            
            d=data.rawdata;
            xind=(numstartupcycles*2+1):2:length(d(1,:));
            yind=(numstartupcycles*2+2):2:length(d(1,:));
            
            unique_samplenames=unique(data.input.samplenames);
            
            for i=1:length(unique_samplenames)
                data.biadata(i).ligand=[];
                data.biadata(i).buffer=[];
                data.biadata(i).name=unique_samplenames{i};
                for j=1:length(data.input.samplenames)
                    if regexp(data.input.samplenames{j},unique_samplenames{i})
                    if regexp(data.input.ligandname{j},'buffer')
                        data.biadata(i).buffer(end+1).cyc_num=j;
                        data.biadata(i).buffer(end).ligandconc=data.input.ligandconc(j,:);
                        data.biadata(i).buffer(end).traw=data.rawdata(:,xind(j));
                        data.biadata(i).buffer(end).yraw=data.rawdata(:,yind(j));
                    else
                        data.biadata(i).ligand(end+1).cyc_num=j;
                        data.biadata(i).ligand(end).ligandname=data.input.ligandname{j};
                        data.biadata(i).ligand(end).ligandconc=data.input.ligandconc(j,:);
                        data.biadata(i).ligand(end).traw=data.rawdata(:,xind(j));
                        data.biadata(i).ligand(end).yraw=data.rawdata(:,yind(j));
                    end
                    end
                end
            end
            
            time=timeintervals;
            %
            for k=1:length(data.biadata)
                for j=1:length(data.biadata(k).buffer)
                    intv={};
                    longer=360;
                    ta=zeros(180,5);
                    td=zeros(longer+1,5);
                    ya=zeros(180,5);
                    yd=zeros(longer+1,5);

                    for i=2:length(time)
                        x1=find(abs(data.biadata(k).buffer(j).traw-time(i-1))==min(abs(data.biadata(k).buffer(j).traw-time(i-1))));
                        x2=find(abs(data.biadata(k).buffer(j).traw-time(i))==min(abs(data.biadata(k).buffer(j).traw-time(i))));
                        intv{end+1}=x1:x2;
                    end
                    intv{end+1}=x2:length(data.biadata(k).buffer(j).traw);

                    for i=1:(length(intv)-1)
                            ta(:,i)=(1:180)';
                            td(1:180,i)=(181:360)';
                            ya(:,i)=data.biadata(k).buffer(j).yraw(intv{i}(1:180));
                            yd(1:180,i)=data.biadata(k).buffer(j).yraw(intv{i}(181:360));

                    end
                    
                    % account for long dissociation of the final binding
                    % event
                    
                    ta(:,end)=(1:180)';
                    td(:,end)=(181:(181+longer))';
                    ya(:,end)=data.biadata(k).buffer(j).yraw(intv{end}(1:180));
                    yd(:,end)=data.biadata(k).buffer(j).yraw(intv{end}(181:(181+longer)));
                    

                    data.biadata(k).buffer(j).ta=ta;
                    data.biadata(k).buffer(j).td=td;
                    data.biadata(k).buffer(j).ya=ya;
                    data.biadata(k).buffer(j).yd=yd;
                    
                end
                
                for j=1:length(data.biadata(k).ligand)
                    intv={};
                    ta=zeros(180,5);
                    td=zeros(longer+1,5);
                    ya=zeros(180,5);
                    yd=zeros(longer+1,5);

                    for i=2:length(time)
                        x1=find(abs(data.biadata(k).ligand(j).traw-time(i-1))==min(abs(data.biadata(k).ligand(j).traw-time(i-1))));
                        x2=find(abs(data.biadata(k).ligand(j).traw-time(i))==min(abs(data.biadata(k).ligand(j).traw-time(i))));
                        intv{end+1}=x1:x2;
                    end
                    intv{end+1}=x2:length(data.biadata(k).ligand(j).traw);
                    
                    for i=1:(length(intv)-1)
                            ta(:,i)=(1:180)';
                            td(1:180,i)=(181:360)';
                            ya(:,i)=data.biadata(k).ligand(j).yraw(intv{i}(1:180));
                            yd(1:180,i)=data.biadata(k).ligand(j).yraw(intv{i}(181:360));

                    end
                    
                    % account for long dissociation of the final binding
                    % event
                    ta(:,end)=(1:180)';
                    td(:,end)=(181:(181+longer))';
                    ya(:,end)=data.biadata(k).ligand(j).yraw(intv{end}(1:180));
                    yd(:,end)=data.biadata(k).ligand(j).yraw(intv{end}(181:(181+longer)));
                    
                    
                    data.biadata(k).ligand(j).ta=ta;
                    data.biadata(k).ligand(j).td=td;
                    data.biadata(k).ligand(j).ya=ya;
                    data.biadata(k).ligand(j).yd=yd;
                    
                    
                    data.biadata(k).rep(j).tabs=data.biadata(k).ligand(j).ta;
                    data.biadata(k).rep(j).tdbs=data.biadata(k).ligand(j).td;
                    
                    % substract sample traces by buffer traces
                    if (length(data.biadata(k).buffer)-length(data.biadata(k).ligand))==0
                    
                    data.biadata(k).rep(j).yabs=data.biadata(k).ligand(j).ya-data.biadata(k).buffer(j).ya;
                    data.biadata(k).rep(j).ydbs=data.biadata(k).ligand(j).yd-data.biadata(k).buffer(j).yd;
                    elseif (length(data.biadata(k).buffer)-length(data.biadata(k).ligand))==1
                        bufferya(:,:,1)=data.biadata(k).buffer(j).ya;
                        bufferya(:,:,2)=data.biadata(k).buffer(j+1).ya;
                        meanbufferya=mean(bufferya,3);
                        
                        bufferyd(:,:,1)=data.biadata(k).buffer(j).yd;
                        bufferyd(:,:,2)=data.biadata(k).buffer(j+1).yd;
                        meanbufferyd=mean(bufferyd,3);
                        
                    data.biadata(k).rep(j).yabs=data.biadata(k).ligand(j).ya-meanbufferya;
                    data.biadata(k).rep(j).ydbs=data.biadata(k).ligand(j).yd-meanbufferyd;                        
                        
                    end
                    data.biadata(k).rep(j).ligandconc=data.biadata(k).ligand(j).ligandconc;
                    
                    % zero baseline data
                    scalefactor=mean(data.biadata(k).rep(j).yabs(1:5,1));
                    data.biadata(k).rep(j).yabs=data.biadata(k).rep(j).yabs-scalefactor;
                    data.biadata(k).rep(j).ydbs=data.biadata(k).rep(j).ydbs-scalefactor;                        
                    
                end      
            end
        end
        
        
        function fitKineticBiadata(data,beta0,datarange,showfit)
            betamat=struct;
            
            for k=1:length(data.biadata)
                if regexp(data.biadata(k).name,'N[0-9]*N[0-9]*')
                    controlid=k;
                end
            end

%             fprintf('\tName, \tka, \tkd, \tKD_kinetic\n');
            for k=1:length(data.biadata)
                
                for j=1:length(data.biadata(k).rep)
                    ta= data.biadata(k).rep(j).tabs;
                    td= data.biadata(k).rep(j).tdbs;
                    
                if k==controlid
                    ya= data.biadata(k).rep(j).yabs;
                    yd= data.biadata(k).rep(j).ydbs;
                else % normalize to negative control
                    if length(data.biadata(controlid).rep)==length(data.biadata(k).rep)
                    ya= data.biadata(k).rep(j).yabs-data.biadata(controlid).rep(j).yabs;
                    yd= data.biadata(k).rep(j).ydbs-data.biadata(controlid).rep(j).ydbs;
                    else 
                    ya= data.biadata(k).rep(j).yabs-data.biadata(controlid).rep(1).yabs;
                    yd= data.biadata(k).rep(j).ydbs-data.biadata(controlid).rep(1).ydbs;
                    end                        
                end
                    ligandconc=data.biadata(k).rep(j).ligandconc*1e-6;
                    
                    cyclename=sprintf('%s-%d',data.biadata(k).name,j);
                    setfig(strcat('kinetic fig: ',cyclename));clf;hold on
                    
                    legendnames={};
                    for i=1:length(ta(1,:))
                        
%                     % show assoc and dissoc as separate colors
%                     plot([ta(:,i)],[ya(:,i)],'-','linewidth',2)
%                     plot([td(td(:,i)~=0,i)],[yd(td(:,i)~=0,i)],'-','linewidth',2)
%                     % show long tail
%                     plot([ta(:,i);td(td(:,i)~=0,i)],[ya(:,i);yd(td(:,i)~=0,i)],'-','linewidth',2)
                    % clip tail
                    plot([ta(1:180,i);td(1:180,i)],[ya(1:180,i);yd(1:180,i)],'-','linewidth',2)
%                     legendnames{end+1}=sprintf('%1.3f %s %s',data.biadata(k).rep(j).ligandconc(i),data.input.ligandconcunit,data.biadata(k).ligand(j).ligandname);
                    legendnames{end+1}=sprintf('%1.3f %s',data.biadata(k).rep(j).ligandconc(i),data.input.ligandconcunit);
                    end
                    legend(legendnames{:},'location','northeast')
                    
                    try
                    [beta,r,J,Sigma,mse,errorparam,robustw]=donlinmultifit(ligandconc(datarange),ta(:,datarange),ya(:,datarange),td(:,datarange),yd(:,datarange),beta0,showfit);
                    betamat(end+1).beta=beta;
                    data.biadata(k).rep(j).ka=betamat(end).beta(1);
                    data.biadata(k).rep(j).kd=betamat(end).beta(2);
                    data.biadata(k).rep(j).kD=betamat(end).beta(2)/betamat(end).beta(1);
                    catch
                        fprintf('bad fit: %s-%d\n',data.biadata(k).name,j)
                    end
                    
%                     xL=get(gca,'XLim');
                    xlim([-100,600])
                    set(gca,'fontsize',24)
                    set(gca,'XTick',[0 200 400 600])
                    ylim([-20,60])
                end
                
            end
        end
        
        
        function fitEquilibriumBiadata(data,ic50guesses,datarange,showfit)
            
            for k=1:length(data.biadata)
                for j=1:length(data.biadata(k).rep)
                    yss=mean(data.biadata(k).rep(j).yabs(end-10:end-3,datarange));
                    yssnormalized=(yss-min(yss))/(max(yss)-min(yss));
%                     setfig('dummy');clf

                    data.biadata(k).rep(j).yssnormalized=yssnormalized;
                
                    ligconc=data.biadata(k).rep(j).ligandconc(datarange)*1e-6;
                    if showfit
                        cyclename=sprintf('%s-%d',data.biadata(k).name,j);
                        setfig(strcat('Eq. IC50 fig: ',cyclename));clf;hold on
                        [hillCoeff, ec50]=doseResponse(ligconc,yssnormalized,ic50guesses);


                        ylabel('Normalized RU')
                        xlabel(sprintf('ligand concentration (%s)',data.input.ligandconcunit))
                        title(strcat(cyclename,' IC50'),'interpreter','none')
                        set(gca,'linewidth',1.5)
                        set(gca,'XScale','log')

                        hold on
                        plot(ligconc,yssnormalized,'o','linewidth',2)
                        text(0.9*round(mean(log10(ligconc))),0.55,sprintf('KD_e_q = %0.2f',ec50),'fontsize',14)

                    else
                        setfig('dummy');clf
                        [hillCoeff, ec50]=doseResponse(ligconc,yssnormalized,ic50guesses);
                    end
                    
                    data.biadata(k).rep(j).KDss=real(ec50);                
                end
            end
        end
        
        function dispFit(data)
            fprintf('\tName, \tka, \tkd, \tKD_kinetic\n');
            for k=1:length(data.biadata)
                for j=1:length(data.biadata(k).rep)
                    cyclename=sprintf('%s-%d',data.biadata(k).name,j);
                    fprintf('%20s, %2.4e, %2.4e, %2.4e, %2.4e\n',cyclename,...
                        data.biadata(k).rep(j).ka,data.biadata(k).rep(j).kd,...
                        data.biadata(k).rep(j).kD,data.biadata(k).rep(j).KDss);

                end
            end

        end
            
        
    end
end