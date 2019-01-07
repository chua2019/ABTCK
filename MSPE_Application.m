clear all; close all;
% load covariance function
addpath ./cov_functions
% add functions
addpath ./function_fold
addpath ./function_fold/conditionalLHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load functions to generate random variables 
addpath ./Random_var; 
% load tree operations 
addpath ./Multi_GP_operations;
addpath ./Multi_tree_NON_nested;
addpath ./help_tree_operations; 
addpath ./Multi_prediction 
addpath ./G_example/wrfdata;

addpath ./G_example/Sent_to_Alex;
%%


warning('off','all')



for jj=1:60

load sgp25km_daily_rain_June_mean
%sgp50km_daily_rain_obs_mean 
%load ./timeaveraged_data/sgp50km_lat_lon2d.mat 
load sgp25km_lat_lon2d 
load sgp25km_lat_lon2d 
load sgp25km_calibr

load sgp50km_daily_rain_June_mean
%sgp50km_daily_rain_obs_mean 
%load ./timeaveraged_data/sgp50km_lat_lon2d.mat 
load sgp50km_lat_lon2d 
load sgp50km_lat_lon2d 
load sgp50km_calibr

xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

n_P50=size(sgp50km_calibr,1);

% DESCRIPTIVE TEXT
p1_50km =10*sgp50km_calibr(1:n_P50,1)/(max(sgp50km_calibr(:,1))-min(sgp50km_calibr(1:n_P50,1)));
%exp((sgp50km_calibr(1:n_P50,1)+(-min(sgp50km_calibr(:,1)))));
p2_50km =10*sgp50km_calibr(1:n_P50,2)/(max(sgp50km_calibr(:,2))-min(sgp50km_calibr(1:n_P50,2)));
%log((sgp50km_calibr(1:n_P50,2)+1));

p3_50km = 10*sgp50km_calibr(1:n_P50,3)/(max(sgp50km_calibr(:,3))-min(sgp50km_calibr(1:n_P50,3)));
p4_50km = 10*sgp50km_calibr(1:n_P50,4)/(max(sgp50km_calibr(:,4))-min(sgp50km_calibr(1:n_P50,4)));
p5_50km = 10*sgp50km_calibr(1:n_P50,5)/(max(sgp50km_calibr(:,5))-min(sgp50km_calibr(1:n_P50,5)));
% p6_12km = sgp12km_calibr(1:n_P25,6) ;
pAll_50kmAll =[p1_50km p2_50km p3_50km p4_50km p5_50km];
pAll_50km=pAll_50kmAll;%(1:49,:);
n_P50=size(pAll_50km,1);

ind=40; % till 59 seem to be ok 
Ally_50km = (squeeze(sgp50km_daily_rain_June_mean(1:n_P50,:,:))) ;%log(squeeze(sgp50km_daily_rain_June_mean(1:n_P50,:,:))+0.5) ;
y_50km = (squeeze(sgp50km_daily_rain_June_mean(ind,:,:))) ;
x1_50km = squeeze(sgp50km_lat_lon2d(1,:,:)) ;
x2_50km = squeeze(sgp50km_lat_lon2d(2,:,:)) ;

All_50km = [x1_50km(:) x2_50km(:) Ally_50km(1:n_P50,:)'];
All_50km(any(isnan(All_50km), 2), :) = [];
Anew_50km = All_50km;
space_in50km=pAll_50km(1:n_P50,:);
space_s50km=Anew_50km(:,1:2);
MeanZt_50ALL=mean((Anew_50km(:,3:(n_P50+2))))';
%figure; scatter(space_in50km(:,1),space_in50km(:,4),20,(MeanZt_50ALL)); colorbar;
n_P25=size(sgp25km_calibr,1);

% DESCRIPTIVE TEXT
p1_25km =10*sgp25km_calibr(1:n_P25,1)/(max(sgp50km_calibr(:,1))-min(sgp50km_calibr(1:n_P50,1)));
%p1_25km = exp((sgp25km_calibr(1:n_P25,1)+(-min(sgp25km_calibr(:,1)))));
p2_25km =10*sgp25km_calibr(1:n_P25,2)/(max(sgp50km_calibr(:,2))-min(sgp50km_calibr(1:n_P50,2)));
%p2_25km =log((sgp25km_calibr(1:n_P25,2)+1));
p3_25km = 10*sgp25km_calibr(1:n_P25,3)/(max(sgp50km_calibr(1:n_P50,3))-min(sgp50km_calibr(1:n_P50,3)));
p4_25km = 10*sgp25km_calibr(1:n_P25,4)/(max(sgp50km_calibr(1:n_P50,4))-min(sgp50km_calibr(1:n_P50,4)));
p5_25km = 10*sgp25km_calibr(1:n_P25,5)/(max(sgp50km_calibr(1:n_P50,5))-min(sgp50km_calibr(1:n_P50,5)));
pAll_25kmAll =[p1_25km p2_25km p3_25km p4_25km p5_25km];
pAll_25km=pAll_25kmAll;%(1:49,:);
n_P25=size(pAll_25km,1);


ind=40; % till 59 seem to be ok 
Ally_25km = (squeeze(sgp25km_daily_rain_June_mean(1:n_P25,:,:))) ;%log(squeeze(sgp25km_daily_rain_June_mean(1:n_P25,:,:))+0.5) ;
%y_25km = (squeeze(sgp25km_daily_rain_June_mean(ind,:,:))) ;
x1_25km = squeeze(sgp25km_lat_lon2d(1,:,:)) ;
x2_25km = squeeze(sgp25km_lat_lon2d(2,:,:)) ;

 
 All_25km = [x1_25km(:) x2_25km(:) Ally_25km(1:n_P25,:)'];
All_25km(any(isnan(All_25km), 2), :) = [];
Anew_25km = All_25km;
space_in25km=pAll_25km(1:n_P25,:);
space_s25km=Anew_25km(:,1:2);
MeanZt_25ALL=mean((Anew_25km(:,3:(n_P25+2))))';
% figure; scatter(space_in25km(:,1),space_in25km(:,4),20,(MeanZt_25ALL)); colorbar;
% figure; scatter(space_in50km(:,1),space_in50km(:,4),20,(MeanZt_50ALL)); colorbar;
%%

space_in25kmALL=pAll_25km(1:n_P25,:);
space_in50km=pAll_50km(1:n_P50,:);
% space_s25km=Anew_25km(:,1:2);
% space_s50km=Anew_50km(:,1:2);
MeanZt_25ALL=mean((Anew_25km(:,3:(n_P25+2))))';
n_inn=40;

rng(jj+35)
%ipick_Pred=randsample(nn_inpt,8);

All_25km = [x1_25km(:) x2_25km(:) Ally_25km(1:n_P25,:)'];
All_25km(any(isnan(All_25km), 2), :) = [];
site_l1=All_25km(:,1:2);
z=ones(size(site_l1,1),1);
xobj=ones(length(z),1);% z can be the output or one of the locations%data(:,6);     % class/categorigal variable
vobj=1;
niter=1000;    % no iterations
[ipick_Pred,xsam,xobs,oL,isam]=cLHS(35,space_in25kmALL,xobj,vobj,niter);

MeanZt_25PredictAA=MeanZt_25ALL(ipick_Pred,:);
space_in25kmPred=space_in25kmALL(ipick_Pred,:);
MeanZt_25=MeanZt_25ALL;
space_in25km=space_in25kmALL;
MeanZt_25(ipick_Pred,:) =[];
space_in25km(ipick_Pred,:) =[];

OLA_pred=sortrows([space_in25kmPred MeanZt_25PredictAA],1);
MeanZt_50=mean(Anew_50km(:,3:(n_P50+2)))'; 

%%

LBlim_aa(1)=min(min(pAll_25km(:,1)),min(pAll_50km(:,1)))-0.1;
LBlim_aa(2)=min(min(pAll_25km(:,2)),min(pAll_50km(:,2)))-0.1;
LBlim_aa(3)=min(min(pAll_25km(:,3)),min(pAll_50km(:,3)))-0.1;
LBlim_aa(4)=min(min(pAll_25km(:,4)),min(pAll_50km(:,4)))-0.1;
LBlim_aa(5)=min(min(pAll_25km(:,5)),min(pAll_50km(:,5)))-0.1;


UBlim_aa(1)=max(max(pAll_25km(:,1)),max(pAll_50km(:,1)))+0.1;
UBlim_aa(2)=max(max(pAll_25km(:,2)),max(pAll_50km(:,2)))+0.1;
UBlim_aa(3)=max(max(pAll_25km(:,3)),max(pAll_50km(:,3)))+0.1;
UBlim_aa(4)=max(max(pAll_25km(:,4)),max(pAll_50km(:,4)))+0.1;
UBlim_aa(5)=max(max(pAll_25km(:,5)),max(pAll_50km(:,5)))+0.1;




site_e=OLA_pred(:,1:5); site_ee=site_e; n_ee=size(OLA_pred,1) ;
MeanZt_25Predict=OLA_pred(:,end);
% mm=49; x1=[LBlim_aa(1):(UBlim_aa(1)-LBlim_aa(1))/mm:UBlim_aa(1)]'; x2=[LBlim_aa(2):(UBlim_aa(2)-LBlim_aa(2))/mm:UBlim_aa(2)]'; 
% [X1,X2]=meshgrid(x1,x2); site_e=[X1(:), X2(:)];  

%%

SITE{1}{1}= space_in50km; SITE{1}{2}=space_in25km;
ZT{1}{1}=MeanZt_50; ZT{1}{2}=MeanZt_25; 
p=5; a_old=1; BB=300000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%        MCMC algorith    %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
funname='gauss'; %funname2='anisotropic';
funname2 = 'constant';%  funname2 ='linear'; % funname2='quadratic'; %
Delta=0.4; 
Dim_x=2;
a_NODES=1;
 th_bound.min=0;
 th_bound.max=1;
 AA_E{1}{1}=2*eye(p); AA_E{1}{2}=.5*eye(p);

 LIM_ALL{1}=[LBlim_aa(1) UBlim_aa(1); LBlim_aa(2) UBlim_aa(2);LBlim_aa(3) UBlim_aa(3); LBlim_aa(4) UBlim_aa(4);LBlim_aa(5) UBlim_aa(5)]; 
 
 
 

 NODES{1}.number =1; %  
 NODES{1}.inter =0; % exterior 
 NODES{1}.parent_child=0; % Left child 
 NODES{1}.direction=1;
 NODES{1}.rul=LIM_ALL{1}(1,2);
 
aalpha=0.2; bbeta=5; %aalpha=0.2;  bbeta=10; %

 PR_CURR{1}{1}=10^(-15); PR_CURR{1}{2}=10^(-15);        
prior{1}.a=2.5; prior{1}.b=2.5;
prior{2}.a=2.5; prior{2}.b=2.5;
%a1=1; b1=10; a2=10; b2=10; probm =0.5;% a1=2; b1=2; a2=20; b2=2; probm =0.5;

% pprop.a1=1;
% pprop.b1=10;%pprop.b1=0.5;%
% pprop.a2=5;
% pprop.b2=2;
% pprop.probm=0.5;

pprop.a1=2;
pprop.b1=1;%pprop.b1=0.5;%
pprop.a2=20;
pprop.b2=2;
pprop.probm=0.5;

% for ii=1:1000
% xx(ii)=rand_mix_gamma(pprop.a1,pprop.b1,pprop.a2,pprop.b2,pprop.probm);
% end
% figure;hist(xx)
% tau1=0.00005; tau2=0.005;
tau1=0.0001; tau2=0.001;%tau2=0.0005;

S=2; low_l=5; a_new=1;

nn{1}{1}=size(ZT{1}{1},1);
nn{1}{2}=size(ZT{1}{2},1);


[Zt_realGP,site_realGP,prCurrGP,AAGP,a_newGP,lim_allGP, KrigSGP,Var_krigSGP,Krig_siteGP, Krig_W_etaMATGP, Krig_site_etaGP]=...
    BCoKrtigNON_nested(PR_CURR,SITE,ZT,nn,AA_E,LIM_ALL,a_new,a_NODES,NODES,pprop, aalpha, bbeta,S,low_l,site_ee,n_ee,BB,Delta,tau1,tau2,funname);

Pred_MeanGP=mean(Krig_W_etaMATGP(1000:BB,:))';

%Pred_Mean-MeanZt_12Predict
%AA_E=AAGP;

MSPE_GP(jj)=mean((Pred_MeanGP-MeanZt_25Predict).^2);
%tau1=0.00005; tau2=0.004;

[Zt_real,site_real,prCurr,AA,a_new,lim_all, KrigS,Var_krigS,Krig_site, Krig_W_etaMAT, Krig_site_eta]=...
    BTMultiNON_nested(PR_CURR,SITE,ZT,nn,AA_E,LIM_ALL,a_new,a_NODES,NODES,pprop, aalpha, bbeta,S,low_l,site_ee,n_ee,BB,Delta,tau1,tau2,funname);


Pred_Mean=mean(Krig_W_etaMAT(1000:BB,:))';

%Pred_Mean-MeanZt_12Predict

MSPE(jj)=mean((Pred_Mean-MeanZt_25Predict).^2)

%ERROR(jj)=mean(mean(abs(Krig_site_eta-site_ee))); MSPE-MSPE_GP;
close all;

colormap(jet)
site_eeR=(max(sgp50km_calibr(:,1))-min(sgp50km_calibr(1:n_P50,1))).*site_ee/10;
  fighdl1_R=figure;
set(fighdl1_R,'DefaultAxesFontSize', 30)
set(findall(fighdl1_R,'type','text'),'fontSize',30)
scatter(site_eeR(:,1),site_eeR(:,2),30,(MeanZt_25Predict),'fill'); colorbar; caxis([min(MeanZt_25Predict) max(MeanZt_25Predict)])
%colormap('bone');
xlabel(['$P_d$'], 'interpreter','latex');
ylabel(['$P_e$'], 'interpreter','latex');
zlabel(['Prec. Mean'], 'interpreter','latex');
set(fighdl1_R, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
 saveas(fighdl1_R, ['Results\Application_predictions\AARPrediction_25_nstart_',num2str(jj),'PdPe.pdf']); 
 saveas(fighdl1_R, ['Results\Application_predictions\AARPrediction_25_nstart_',num2str(jj),'PdPe.fig']); 
 
   fighdl1_R=figure;
set(fighdl1_R,'DefaultAxesFontSize', 30)
set(findall(fighdl1_R,'type','text'),'fontSize',30)
scatter(site_eeR(:,1),site_eeR(:,3),30,(MeanZt_25Predict),'fill'); colorbar; caxis([min(MeanZt_25Predict) max(MeanZt_25Predict)])
%colormap('bone');
xlabel(['$P_d$'], 'interpreter','latex');
ylabel(['$P_t$'], 'interpreter','latex');
zlabel(['Prec. Mean'], 'interpreter','latex');
set(fighdl1_R, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
 saveas(fighdl1_R, ['Results\Application_predictions\AARPrediction_25_nstart_',num2str(jj),'PdPt.pdf']); 
 saveas(fighdl1_R, ['Results\Application_predictions\AARPrediction_25_nstart_',num2str(jj),'PdPt.fig']); 
   fighdl2_R=figure;
set(fighdl2_R,'DefaultAxesFontSize', 30)
Krig_site_etaGPR=(max(sgp50km_calibr(:,1))-min(sgp50km_calibr(1:n_P50,1))).*Krig_site_etaGP/10;
set(findall(fighdl2_R,'type','text'),'fontSize',30)
scatter(Krig_site_etaGPR(:,1),Krig_site_etaGPR(:,2),30,(Pred_MeanGP),'fill'); colorbar; caxis([min(MeanZt_25Predict) max(MeanZt_25Predict)])
%colormap('bone');
xlabel(['$P_d$'], 'interpreter','latex');
ylabel(['$P_e$'], 'interpreter','latex');
zlabel(['Prec. Mean'], 'interpreter','latex');
set(fighdl2_R, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
 saveas(fighdl2_R, ['Results\Application_predictions\AAGPPrediction_25_nstart_',num2str(jj),'PdPe.pdf']); 
 saveas(fighdl2_R, ['Results\Application_predictions\AAGPPrediction_25_nstart_',num2str(jj),'PdPe.fig']); 

   fighdl2_R=figure;
set(fighdl2_R,'DefaultAxesFontSize', 30)
Krig_site_etaGPR=(max(sgp50km_calibr(:,1))-min(sgp50km_calibr(1:n_P50,1))).*Krig_site_etaGP/10;
set(findall(fighdl2_R,'type','text'),'fontSize',30)
scatter(Krig_site_etaGPR(:,1),Krig_site_etaGPR(:,3),30,(Pred_MeanGP),'fill'); colorbar; caxis([min(MeanZt_25Predict) max(MeanZt_25Predict)])
%colormap('bone');
xlabel(['$P_d$'], 'interpreter','latex');
ylabel(['$P_t$'], 'interpreter','latex');
zlabel(['Prec. Mean'], 'interpreter','latex');
set(fighdl2_R, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
 saveas(fighdl2_R, ['Results\Application_predictions\AAGPPrediction_25_nstart_',num2str(jj),'PdPt.pdf']); 
 saveas(fighdl2_R, ['Results\Application_predictions\AAGPPrediction_25_nstart_',num2str(jj),'PdPt.fig']); 

   fighdl3_R=figure;
set(fighdl3_R,'DefaultAxesFontSize', 30)
Krig_site_etaR=(max(sgp50km_calibr(:,1))-min(sgp50km_calibr(1:n_P50,1))).*Krig_site_eta/10;
set(findall(fighdl3_R,'type','text'),'fontSize',30)
scatter(Krig_site_etaR(:,1),Krig_site_etaR(:,2),30,(Pred_Mean),'fill'); colorbar; caxis([min(MeanZt_25Predict) max(MeanZt_25Predict)])
%colormap('bone');
xlabel(['$P_d$'], 'interpreter','latex');
ylabel(['$P_e$'], 'interpreter','latex');
zlabel(['Prec. Mean'], 'interpreter','latex');
set(fighdl3_R, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
 saveas(fighdl3_R, ['Results\Application_predictions\AATGPPrediction_25_nstart_',num2str(jj),'PdPe.pdf']); 
 saveas(fighdl3_R, ['Results\Application_predictions\AATGPPrediction_25_nstart_',num2str(jj),'PdPe.fig']); 

    fighdl3_R=figure;
set(fighdl3_R,'DefaultAxesFontSize', 30)
Krig_site_etaR=(max(sgp50km_calibr(:,1))-min(sgp50km_calibr(1:n_P50,1))).*Krig_site_eta/10;
set(findall(fighdl3_R,'type','text'),'fontSize',30)
scatter(Krig_site_etaR(:,1),Krig_site_etaR(:,3),30,(Pred_Mean),'fill'); colorbar; caxis([min(MeanZt_25Predict) max(MeanZt_25Predict)])
%colormap('bone');
xlabel(['$P_d$'], 'interpreter','latex');
ylabel(['$P_t$'], 'interpreter','latex');
zlabel(['Prec. Mean'], 'interpreter','latex');
set(fighdl3_R, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
 saveas(fighdl3_R, ['Results\Application_predictions\AATGPPrediction_25_nstart_',num2str(jj),'PdPt.pdf']); 
 saveas(fighdl3_R, ['Results\Application_predictions\AATGPPrediction_25_nstart_',num2str(jj),'PdPt.fig']); 

pause(.0001);
end


% save('MSPE_full35')
% 
% save('MSPE_GP_full35')



