


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
addpath ./Multi_tree_nested;
addpath ./help_tree_operations; 
addpath ./Multi_prediction
% Load PDE
addpath ./G_example/Ex_2_c; 
%%
%warning('off', 'MSGID')
warning('off')

% % Generate data from a two dimentional input and two dimentional output.
% Low 
 R_low =heat_transfer_simulT(0.07, 0, 1); 
%Intermidiate 
  R_int =heat_transfer_simulT(0.025, 1, 1); 
%High
 R_high = heat_transfer_simulT(0.015,2, 1);
% % First example
 u0 = interpolateSolution(R_low,[0.55 1]',[2 3]') ;
 uu0 = interpolateSolution(R_int,[0.55 1]',[2 3]') ;
 uuu0 = interpolateSolution(R_high,[0.55 1]',[2 3]');

  
 %%
 Dim_x=2;
lim_all=[0 1; 0 3];
rng(111); 
n1=220; n2=130; n3=65;% no of samples

LH_S1=lhsdesign(n1,Dim_x)';
site_l1=[diag(lim_all(:,1))*ones(Dim_x, n1)+diag((lim_all(:,2)-lim_all(:,1)))*LH_S1]';
z=ones(n1,1);
xobj=ones(length(z),1);% z can be the output or one of the locations%data(:,6);     % class/categorigal variable
vobj=1;%[1 2 3 4];     % value of the class
    
niter=10000;    % no iterations
[ipick,xsam,xobs,oL,isam]=cLHS(n2,site_l1,xobj,vobj,niter);
site_l2=site_l1(ipick,:);
%plot coordinates of sampled data
figure;
plot(site_l1(:,1), site_l1(:,2),'.');

hold on;
plot(site_l2(:,1), site_l2(:,2), 'o');
hold off;
xobj2=ones(length(site_l2),1);
[ipick2,xsam2,xobs2,oL2,isam2]=cLHS(n3,site_l2,xobj2,vobj,niter);
site_l3=site_l2(ipick2,:);
%plot coordinates of sampled data
fighdlA=figure;
     set(findall(gcf,'type','text'),'fontSize',20)
     set(gcf,'DefaultAxesFontSize', 20)

plot(site_l1(:,1), site_l1(:,2),'.','MarkerSize', 15);

hold on;
plot(site_l2(:,1), site_l2(:,2), 'o','MarkerSize', 10);
hold on;
plot(site_l3(:,1), site_l3(:,2), '+','MarkerSize', 10);
hold off;
set(fighdlA, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);

saveas(fighdlA, ['Results\Simulation1\Design_plot',num2str(n3),'.pdf']);

% Sites 
SITE{1}{1}= site_l1; SITE{1}{2}=site_l2; SITE{1}{3}=site_l3;
% Low
ZT{1}{1} = interpolateSolution(R_low,SITE{1}{1}(:,1),SITE{1}{1}(:,2))+norm_r(0,0.4,n1);%+norm_r(0,0.5,n1);%+norm_r(0,0.01,m)
% Interm
ZT{1}{2}= interpolateSolution(R_int,SITE{1}{2}(:,1),SITE{1}{2}(:,2))+norm_r(0,0.3,n2);%+norm_r(0,0.1,n2);%norm_r(0,0.01,n)
% High
ZT{1}{3}= interpolateSolution(R_high,SITE{1}{3}(:,1),SITE{1}{3}(:,2))+norm_r(0,0.2,n3);%+norm_r(0,0.5,n3);%norm_r(0,0.01,n)


fighdl6_R=figure;
scatter3(SITE{1}{3}(:,1),SITE{1}{3}(:,2), ZT{1}{3});
xlabel(['x'], 'interpreter','latex');
ylabel(['$t_1$'], 'interpreter','latex');
zlabel(['z'], 'interpreter','latex');
hold off; set(fighdl6_R, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);

figure; 


plot(SITE{1}{3}(:,1),SITE{1}{3}(:,2),'.');

% tasd=lhsdesign(n3,Dim_x);
% figure; plot(tasd(:,1),tasd(:,2),'.');
%%
mm=89; x1=[0:1/mm:1]'; x2=[0:3/mm:3]'; 
[X1,X2]=meshgrid(x1,x2); site_e=[X1(:), X2(:)];  n_e =(mm+1)^2;
p=2; a_old=1; BB=20000;
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
AA_E{1}{1}=1*eye(p); AA_E{1}{2}=1*eye(p); AA_E{1}{3}=1*eye(p);
LIM_ALL{1}=[-0.001 1.001; -0.001 3.001]; 

 NODES{1}.number =1; %  
 NODES{1}.inter =0; % exterior 
 NODES{1}.parent_child=0; % Left child 
 NODES{1}.direction=1;
 NODES{1}.rul=LIM_ALL{1}(1,2);
 
 aalpha=0.2; bbeta=5; % aalpha=0.2;  bbeta=10; % aalpha=0.8; bbeta=2; %
mm_e=59; x1=[0:1/mm_e:1]'; x2 =[0:3/mm_e:3]'; [X1,X2]=meshgrid(x1,x2);
% [X1,X2]=meshgrid(x1,x2);
site_e_eta=[X1(:), X2(:)];
 PR_CURR{1}{1}=10^(-15); PR_CURR{1}{2}=10^(-15); PR_CURR{1}{3}=10^(-15); 

 
 tau1=0.1;
 tau2=0.05;
 prior{1}.a=2.5; prior{1}.b=2.5;
 prior{2}.a=2.5; prior{2}.b=2.5;


%a1=1; b1=10; a2=10; b2=10; probm =0.5;% a1=2; b1=2; a2=20; b2=2; probm =0.5;
% pprop.a1=2;
% pprop.b1=1;%pprop.b1=0.5;%
% pprop.a2=10;
% pprop.b2=2;
% pprop.probm=0.5;
pprop.a1=1;
pprop.b1=10;%pprop.b1=0.5;%
pprop.a2=5;
pprop.b2=2;
pprop.probm=0.5;

S=3; low_l=7; a_new=1;

nn{1}{1}=size(ZT{1}{1},1);
nn{1}{2}=size(ZT{1}{2},1);
nn{1}{3}=size(ZT{1}{3},1);



[Zt,site,prCurr,AA,a_new,lim_all, KrigS,Var_krigS,Krig_rho, Krig_W_etaMAT, Krig_site_eta]=...
    BTMulti(PR_CURR,SITE,ZT,nn,AA_E,LIM_ALL,a_new,a_NODES,NODES,pprop, aalpha, bbeta,S,low_l,site_e,n_e,BB,tau1,tau2,Delta,funname);

hold on 
for ii=1:a_new
     hold on; 
    l1=lim_all{ii}(1,1); l2=lim_all{ii}(1,2);        
    %l3=lim_all{ii}((Dim_x-p),1); l4=lim_all{ii}((Dim_x-p),2);
    l3=lim_all{ii}((Dim_x),1); l4=lim_all{ii}((Dim_x),2);
%text(xcurr(1,i), xcurr(2,i),num2str(i), )
%text((l1+l2)/2,(l3+l4)/2,num2str(ii), 'Color','r');
    plot([l1 l1],[l3 l4]);plot([l2 l2],[l3 l4]);    
    plot([l1 l2], [l3 l3]); plot([l1 l2],[l4 l4]);
end   

hold off
% 
% figure; plot(S_AA(1,:,1,1))
% figure; plot(S_AA(1,:,2,2))
% figure; plot(S_AA(2,:,1,1))
% figure; plot(S_AA(2,:,2,2))
% figure; plot(BETA_v(1,:,1))
% figure; plot(BETA_v(2,:,1))
% figure; plot(BETA_v(2,:,2))
% figure; plot(SIG_MAT(2,:))
% figure; plot(prCurr_MAT(2,:))

size(Krig_W_etaMAT)
mm=89;
Kriging_W_eta =mean(Krig_W_etaMAT(1000:BB,:))';%Krig_W_etaMAT(10000,:)'; % 
fighdl32=figure; 

     set(findall(gcf,'type','text'),'fontSize',20)
     set(gcf,'DefaultAxesFontSize', 20)
XX1=[]; XX2=[]; KW_U=[];
for jj=1:(mm+1)
     KW_U=[KW_U Kriging_W_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)))];
     XX1=[XX1 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),1)];
     XX2=[XX2 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),2)];    
end
 meshc(XX2,XX1,KW_U);% mesh(XX1,XX2,KW_U);
 
     %title 'Bayesian treed multiple co-kriging prediction'
    orient landscape
    xlabel('$x_2$','interpreter','latex')
    ylabel('$x_1$','interpreter','latex')
    zlabel('$u(x_1,x_2)$','interpreter','latex')
%  xlabel(['x'], 'interpreter','latex')
% ylabel(['t_1'], 'interpreter','latex')
% zlabel(['$\eta$'], 'interpreter','latex')
set(fighdl32, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
%     set(gcf, 'PaperPosition', [0 0 25 25]); 
%     set(gcf, 'PaperSize', [25 25]); 
%     set(gca,'FontSize',20)
%  
%       set(findall(gcf,'type','text'),'fontSize',20)
%      set(gcf,'DefaultAxesFontSize', 20)
view(130, 30);
saveas(fighdl32, ['Results\Simulation1\3Level3sec1GorgiosEx1_n2',num2str(n2),'_n3_',num2str(n3),'.pdf']); 
saveas(fighdl32, ['Results\Simulation1\3Level3sec1GorgiosEx1_n2',num2str(n2),'_n3_',num2str(n3),'.fig']); 
%saveas(fighdl32, ['Results\Simulation1\GP1GorgiosEx1_n_',num2str(n),'.pdf']); 





fighdl1=figure;      
set(findall(gcf,'type','text'),'fontSize',20)
set(gcf,'DefaultAxesFontSize', 20);
     imagesc(XX1(:),XX2(:),KW_U); colorbar; set(gca,'YDir','normal'); %caxis([-0.5 2.5]); 
set(fighdl1, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
 %view(130, 30);
% saveas(fighdl1, ['Results\Simulation1\Level3IMGorgiosEx1_n2_',num2str(n2),'.pdf']); 
% saveas(fighdl1, ['Results\Simulation1\Level3IMGorgiosEx1_n2_',num2str(n2),'.fig']); 

%%


Real_val =  interpolateSolution(R_high,Krig_site_eta(:,1),Krig_site_eta(:,2)) ;
MSPE=nanmean(((Kriging_W_eta-Real_val).^2));
%%

size(Krig_rho{2})
Kriging_W_rho =mean(Krig_rho{2})';%Krig_W_etaMAT(10000,:)'; % 
fighdl32=figure; 
     set(findall(gcf,'type','text'),'fontSize',20)
     set(gcf,'DefaultAxesFontSize', 20)
XX1=[]; XX2=[]; KW_Rho=[];
for jj=1:(mm+1)
     KW_Rho=[KW_Rho Kriging_W_rho((1+(jj-1)*(mm+1)):(jj*(mm+1)))];
     XX1=[XX1 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),1)];
     XX2=[XX2 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),2)];    
end
%  meshc(XX2,XX1,KW_Rho);% mesh(XX1,XX2,KW_U);
%  
%      %title 'Bayesian treed multiple co-kriging prediction'
%     orient landscape
%     xlabel('$x_2$','interpreter','latex')
%     ylabel('$x_1$','interpreter','latex')
%     zlabel('$u(x_1,x_2)$','interpreter','latex')
% %  xlabel(['x'], 'interpreter','latex')
% % ylabel(['t_1'], 'interpreter','latex')
% % zlabel(['$\eta$'], 'interpreter','latex')
% %set(fighdl32, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
%     set(gcf, 'PaperPosition', [0 0 25 25]); 
%     set(gcf, 'PaperSize', [25 25]); 
%     set(gca,'FontSize',20)
%  view(130, 30);

fighdl1=figure;      
set(findall(gcf,'type','text'),'fontSize',20)
set(gcf,'DefaultAxesFontSize', 20);
     imagesc(XX1(:),XX2(:),KW_Rho); colorbar; set(gca,'YDir','normal'); %caxis([-0.5 2.5]); 
set(fighdl1, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
% saveas(fighdl1, ['Results\Simulation1\Level3Rhos12ec2GorgiosEx1_n_',num2str(n3),'.fig']);
% saveas(fighdl1, ['Results\Simulation1\Level3Rhos12ec2GorgiosEx1_n_',num2str(n3),'.pdf']);

%%



size(Krig_rho{3})
Kriging_W_rho =mean(Krig_rho{3})';%Krig_W_etaMAT(10000,:)'; % 
fighdl32=figure; 
     set(findall(gcf,'type','text'),'fontSize',20)
     set(gcf,'DefaultAxesFontSize', 20)
XX1=[]; XX2=[]; KW_Rho=[];
for jj=1:(mm+1)
     KW_Rho=[KW_Rho Kriging_W_rho((1+(jj-1)*(mm+1)):(jj*(mm+1)))];
     XX1=[XX1 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),1)];
     XX2=[XX2 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),2)];    
end
%  meshc(XX2,XX1,KW_Rho);% mesh(XX1,XX2,KW_U);
%  
%      %title 'Bayesian treed multiple co-kriging prediction'
%     orient landscape
%     xlabel('$x_2$','interpreter','latex')
%     ylabel('$x_1$','interpreter','latex')
%     zlabel('$u(x_1,x_2)$','interpreter','latex')
% %  xlabel(['x'], 'interpreter','latex')
% % ylabel(['t_1'], 'interpreter','latex')
% % zlabel(['$\eta$'], 'interpreter','latex')
% %set(fighdl32, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
%     set(gcf, 'PaperPosition', [0 0 25 25]); 
%     set(gcf, 'PaperSize', [25 25]); 
%     set(gca,'FontSize',20)
%  view(130, 30);

fighdl1=figure;      
set(findall(gcf,'type','text'),'fontSize',20)
set(gcf,'DefaultAxesFontSize', 20);
     imagesc(XX1(:),XX2(:),KW_Rho); colorbar; set(gca,'YDir','normal'); %caxis([-0.5 2.5]); 
set(fighdl1, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
% saveas(fighdl1, ['Results\Simulation1\Level3Rhos23ec2GorgiosEx1_n_',num2str(n3),'.fig']);
% saveas(fighdl1, ['Results\Simulation1\Level3Rhos23ec2GorgiosEx1_n_',num2str(n3),'.pdf']);





