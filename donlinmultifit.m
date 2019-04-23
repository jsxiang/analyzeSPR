function [beta,r,J,Sigma,mse,errorparam,robustw]=donlinmultifit(C,x1,y1,x2,y2,beta0,showfit)
addpath ./nlinmultifit

% Prepare input for NLINMULTIFIT
x_cell={};
y_cell={};
mdl_cell={};

% Define fitting functions and parameters
for i=1:length(C)
% ka=beta(1);
% kd=beta(2);
% Rmax=beta(3);
% Ri=beta(4);
% Rinf=beta(5);

% from O'Shannessy et al, 
mdl_cell{end+1} = @(beta,x) (C(i).*beta(1).*beta(3))./(C(i).*beta(1) + beta(2)) .* (1-exp(-((C(i).*beta(1)+beta(2)).*(x)))) + beta(4); % association phase
mdl_cell{end+1} = @(beta,x) ((C(i).*beta(1).*beta(3))./(C(i).*beta(1) + beta(2)) + beta(4)) .* exp(-beta(2).*(x-x(1))) + beta(5); % dissociation phase


x_cell{end+1}=x1(:,i);
x2(x2(:,i)==0,i)=nan; % uncomment for SCK analysis, removing zeros
x_cell{end+1}=x2(:,i);

y_cell{end+1}=y1(:,i);
y2(x2(:,i)==0,i)=nan; % uncomment for SCK analysis, removing zeros
y_cell{end+1}=y2(:,i);

end

beta0all=[beta0 ones(1,2*i)];


% perform fitting
[beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit(x_cell, y_cell, mdl_cell, beta0all);
                
% 
% Calculate model predictions and confidence intervals
for i=1:length(C)

% Calculate parameter confidence intervals
ci = nlparci(beta,r,'Jacobian',J);

hold all;
box on;

    [ypred1,delta1] = nlpredci(mdl_cell{i*2-1},x1(:,i),beta,r,'covar',Sigma);
    [ypred2,delta2] = nlpredci(mdl_cell{i*2},x2(~isnan(x2(:,i)),i),beta,r,'covar',Sigma);
    if showfit
    plot([x1(:,i); x2(~isnan(x2(:,i)),i)],[ypred1; ypred2],'Color',[0.2 0.2 0.5],'linewidth',1.5);
    plot(x1(:,i),ypred1+delta1,'Color',[0.5 0.5 0.5],'LineStyle',':','linewidth',1);
    plot(x1(:,i),ypred1-delta1,'Color',[0.5 0.5 0.5],'LineStyle',':','linewidth',1);
    plot(x2(~isnan(x2(:,i)),i),ypred2+delta2,'Color',[0.5 0.5 0.5],'LineStyle',':','linewidth',1);
    plot(x2(~isnan(x2(:,i)),i),ypred2-delta2,'Color',[0.5 0.5 0.5],'LineStyle',':','linewidth',1);
    end
    
set(gca,'linewidth',1.5)
set(gca,'fontsize',16)
ylabel('SPR Response Unit')
xlabel('Time (seconds)')

end
t=sprintf('KD_k_i_n_e_t_i_c = %0.3e M',beta(2)/beta(1));
text(x2(1),-3.5,t,'fontsize',14)

beta;
end

                
                
                
                
                
                
                
                
                
                