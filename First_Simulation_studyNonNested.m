

clear all; close all;
% load covariance function
addpath ./cov_functions
% add functions
addpath ./function_fold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load functions to generate random variables 
addpath ./Random_var; 
% load tree operations 
addpath ./Multi_GP_operations;
addpath ./Multi_tree_NON_nested;
addpath ./help_tree_operations; 
addpath ./Multi_prediction
% Load separable model MCMC_operations
%addpath ./GP_operations; 
%%
% warning('off', 'MSGID')
warning('off')
 funname0='first';
% % Generate data from a two dimentional input and two dimentional output.
 x1=[-2:(1/20):6]; x2=[-2:(1/20):6];
 %x=[x1' x2']; sort('descend');  [X1] = meshgrid(x);
 [X1,X2] = meshgrid(x1,x2); 
Xxx= [reshape(X1,161^2,1),reshape(X2, 161^2,1)];
% % First example
 y1=2*(funn_ALL(Xxx,funname0)+1)+0.5*(exp(sin((0.9.*((Xxx(:,1)+2)/8+0.48)).^(10))));%+0.5*(Xxx(:,1)/3).*sin(10*cos(2*(Xxx(:,1)/12)*pi))+1+norm_r(0,0.001,size(Xxx,1));%+0.1*sin(10*cos(2*(Xxx(:,1)/12)*pi));%+0.1*cos(2*(Xxx(:,2))*pi)+0.5+0.1*Xxx(:,2);
  %figure; scaterplot(x1,x2,y1); set(gca,'YDir','normal') 
  
  Dim_x=2;
mm=sqrt(size(Xxx,1))-1;
  XY_eta= [Xxx, y1]; 
XY_11_eta=sortrows(XY_eta,2);
XY_12_eta=sortrows(XY_11_eta,1);
Kriging_W_etaS1=XY_12_eta(:,(Dim_x+1):(Dim_x+1));
Krig_site_eta1=XY_12_eta(:,1:Dim_x);
    
%figure
XX11=[]; XX12=[]; KW_U1=[];
for jj=1:(mm+1)
     KW_U1=[KW_U1 Kriging_W_etaS1((1+(jj-1)*(mm+1)):(jj*(mm+1)))];
     XX11=[XX11 Krig_site_eta1((1+(jj-1)*(mm+1)):(jj*(mm+1)),1)];
     XX12=[XX12 Krig_site_eta1((1+(jj-1)*(mm+1)):(jj*(mm+1)),Dim_x)];    
end

 figure;  imagesc(XX11(:),XX12(:),KW_U1); colorbar; set(gca,'YDir','normal'); %caxis([-0.5 2.5]);   
  
  
  %%

   y2=4*(funn_ALL(Xxx,funname0))+0.1*exp(sin((0.9.*(Xxx(:,1)/6+0.48)).^(10)))+0.5; %+0.5*(Xxx(:,1)/3).*sin(10*cos(2*(Xxx(:,1)/12)*pi))+1+norm_r(0,0.001,size(Xxx,1));%+0.1*sin(10*cos(2*(Xxx(:,1)/12)*pi));%+0.1*cos(2*(Xxx(:,2))*pi)+0.5+0.1*Xxx(:,2);
  %figure; scaterplot(x1,x2,y1); set(gca,'YDir','normal') 
  
  Dim_x=2;
mm=sqrt(size(Xxx,1))-1;
  XY_eta= [Xxx, y2]; 
XY_11_eta=sortrows(XY_eta,2);
XY_12_eta=sortrows(XY_11_eta,1);
Kriging_W_etaS1=XY_12_eta(:,(Dim_x+1):(Dim_x+1));
Krig_site_eta1=XY_12_eta(:,1:Dim_x);
    
%figure
XX11=[]; XX12=[]; KW_U2=[];
for jj=1:(mm+1)
     KW_U2=[KW_U2 Kriging_W_etaS1((1+(jj-1)*(mm+1)):(jj*(mm+1)))];
     XX11=[XX11 Krig_site_eta1((1+(jj-1)*(mm+1)):(jj*(mm+1)),1)];
     XX12=[XX12 Krig_site_eta1((1+(jj-1)*(mm+1)):(jj*(mm+1)),Dim_x)];    
end

 figure;  imagesc(XX11(:),XX12(:),KW_U2); colorbar; set(gca,'YDir','normal'); %caxis([-0.5 2.5]); 
  
 %%
funname0='first';
lim_all=[-2 6; -2 6]; 
rng(112); 
n=40;
LH_S1=lhsdesign(n,2)';
Dim_x=2; Dim_y=2; M=Dim_y; 
m=120;
rng(111)
LH_S=lhsdesign(m,Dim_x)';
rng(112)
SITE1=[diag(lim_all(:,1))*ones(Dim_x, m)+diag((lim_all(:,2)-lim_all(:,1)))*LH_S]';
x=[diag(lim_all(:,1))*ones(Dim_x, n)+diag((lim_all(:,2)-lim_all(:,1)))*LH_S1]';   
delta_x=0.5*exp(sin((0.9.*(x(:,1)/6+0.48)).^(10)))+1;%(x(:,1)/3).*sin(10*cos(2*(x(:,1)/12))*pi)+1;%sin(2*(x(:,1)/12)*pi)+0.2*cos(2*(x(:,2)/2)*pi)+0.5*x(:,2); %discrepancy %ones(size(x,1),1);%0.7;% 0.7;%0.5*sin(2*(x(:,1)/3)*pi)+0.1*cos(2*(x(:,2)/6)*pi)+0.5;%
% First example
 data_y=2*funn_ALL(x,funname0)+delta_x+0.2+norm_r(0,0.001,n);%funnWilCal2d_Dis(x,Theta_v)+delta_x+norm_r(0,0.1,n);+norm_r(0,0.01,n)
% A second equation
var_dat=var(data_y);
% First example
data_eta = 4*funn_ALL(SITE1,funname0)+0.2*exp(sin((0.9.*(SITE1(:,1)/6+0.48)).^(10)))+0.5+norm_r(0,0.001,m); %funnWilCal2d_Dis(SITE(:,1),SITE(:,2)) +norm_r(0,0.02,m); % +norm_r(0,0.2,m)%+norm_r(0,0.01,m)
% A second equation
kk=6;
fighdl6_R=figure;
set(fighdl6_R,'DefaultAxesFontSize', 20)
set(findall(fighdl6_R,'type','text'),'fontSize',20)
%mesh(x1,x2,y1);%colorbar; caxis([-0.001 0.001]);  
hold on; scatter3(SITE1(:,1),SITE1(:,2),data_eta);
scatter3(x(:,1),x(:,2), data_y,200,'r.');
xlabel(['x'], 'interpreter','latex');
ylabel(['$t_1$'], 'interpreter','latex');
zlabel(['z'], 'interpreter','latex');
hold off; set(fighdl6_R, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);

%%
mm=49; x1=[-2:8/mm:6]'; x2=[-2:8/mm:6]'; 
[X1,X2]=meshgrid(x1,x2); site_e=[X1(:), X2(:)];  n_e =(mm+1)^2;
SITE{1}{1}= SITE1; SITE{1}{2}=x;
ZT{1}{1}=data_eta; ZT{1}{2}=data_y; 
p=2; a_old=1; BB=10000;
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
  mu_Z1=mean(data_eta); 
 Mu_Z1= kron(mu_Z1,ones(m,1));
 mu_Z2=mean(data_y); 
 Mu_Z2= kron(mu_Z2,ones(n,1)); 
 AA_E{1}{1}=1*eye(p); AA_E{1}{2}=1*eye(p);
LIM_ALL{1}=[-2.001 6.001; -2.001 6.001]; 
 mu_Z1=mean(data_eta); 
 NODES{1}.number =1; %  
 NODES{1}.inter =0; % exterior 
 NODES{1}.parent_child=0; % Left child 
 NODES{1}.direction=1;
 NODES{1}.rul=LIM_ALL{1}(1,2);
 
aalpha=0.05;  bbeta=20; %aalpha=0.8; bbeta=2; 
mm_e=59; x1=[-2:8/mm_e:6]'; x2 =[-2:8/mm_e:6]'; [X1,X2]=meshgrid(x1,x2);
% [X1,X2]=meshgrid(x1,x2);
site_e_eta=[X1(:), X2(:)];
 PR_CURR{1}{1}=10^(-15); PR_CURR{1}{2}=10^(-15);        
prior{1}.a=2.5; prior{1}.b=2.5;
prior{2}.a=2.5; prior{2}.b=2.5;
tau1=0.05;
tau2=0.0005;

%pprop.a1=1; pprop.b1=10; pprop.a2=10; pprop.b2=10; pprop.probm =0.5;% a1=2; b1=2; a2=20; b2=2; probm =0.5;
pprop.a1=2;
pprop.b1=1;%pprop.b1=0.5;%
pprop.a2=20;
pprop.b2=2;
pprop.probm=0.5;


S=2; low_l=3; a_new=1;

nn{1}{1}=size(ZT{1}{1},1);
nn{1}{2}=size(ZT{1}{2},1);


[Zt_real,site_real,prCurr,AA,a_new,lim_all, KrigS,Var_krigS,Krig_site, Krig_W_etaMAT, Krig_site_eta]=...
    BTMultiNON_nested(PR_CURR,SITE,ZT,nn,AA_E,LIM_ALL,a_new,a_NODES,NODES,pprop, aalpha, bbeta,S,low_l,site_e,n_e,BB,Delta,tau1,tau2,funname);

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

hold off;
mm=49;
Kriging_W_eta =mean(Krig_W_etaMAT)';%Krig_W_etaMAT(10000,:)'; % 
fighdl32=figure; 
     set(findall(gcf,'type','text'),'fontSize',20)
     set(gcf,'DefaultAxesFontSize', 20)
XX1=[]; XX2=[]; KW_U=[];
for jj=1:(mm+1)
     KW_U=[KW_U Kriging_W_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)))];
     XX1=[XX1 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),1)];
     XX2=[XX2 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),2)];    
end
 mesh(XX1,XX2,KW_U);
 xlabel(['x'], 'interpreter','latex')
ylabel(['t_1'], 'interpreter','latex')
zlabel(['$\eta$'], 'interpreter','latex')
set(fighdl32, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);


fighdl1=figure;      
set(findall(gcf,'type','text'),'fontSize',20)
set(gcf,'DefaultAxesFontSize', 20);
     imagesc(XX1(:),XX2(:),KW_U); colorbar; set(gca,'YDir','normal'); %caxis([-0.5 2.5]); 
set(fighdl1, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
%saveas(fighdl1, ['Results\Simulation1\NSBTsample_n_',num2str(m),'_nstart_',num2str(n),'.pdf']); 
%saveas(fighdl1, ['Results\Simulation1\Sttsample_n_',num2str(m),'_nstart_',num2str(n),'.pdf']); 

fighdl6= figure;   
plot(x(:,1),x(:,2),'r*');
 hold on; 
plot(SITE1(:,1),SITE1(:,2),'bo');
xlabel('x_1');
ylabel('x_2');
hold off;
set(gcf,'NextPlot','add');axes; h = title('Sampling'); set(gca,'Visible','off'); set(h,'Visible','on'); 
set(fighdl6, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
%saveas(fighdl6, ['Results\Simulation1\sample_n_',num2str(size(SITE1_NEW(:,1))),'_nstart_',num2str(n),'.pdf']); 

%%


Real_val = 2*funn_ALL(Krig_site_eta,funname0)+0.5*exp(sin((0.9.*(Krig_site_eta(:,1)/6+0.48)).^(10)))+1.2;
MSPE=mean(((Kriging_W_eta-Real_val).^2));
MSPE 


