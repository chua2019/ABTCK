function [S_AA,prCurr_MAT,BETA_v,SIG_MAT]= Multi_GIBBS_pred(Zt,site,AA,S,Delta,BB,site_e,tau1,tau2)
    Dim_x=size(site{1},2); 
S_AA=zeros(S,BB,Dim_x,Dim_x);
SIG_MAT=  zeros(S,BB);
prCurr_MAT=zeros(S,BB,1);

%    case 'constant'
BETA_v= zeros(S,BB,2);
for i=1:BB
[AA, SIG, prCurr, mu, beta_X]= Multi_MH_s(Zt,site,AA,S,Delta,tau1,tau2);
%beta_X{2}
beta_X{1}=[beta_X{1}; 1];
for t=1:S;
     S_AA(t,i,:,:) = AA{t};
     SIG_MAT(t,i) =SIG{t}; 
     BETA_v(t,i,:)= beta_X{t};            
     prCurr_MAT(t,i)=prCurr{t};
     
     
     beta{t}=beta_X{t}(1);
     rho{t}=beta_X{t}(2);    
end

[Kriging_W, Var_krig]=Krig_Multi(Zt,beta,rho,SIG,AA,site,site_e,S,tau2);

M=1;mm=sqrt(size(site_e,1))-1;
XY_eta= [site_e, Kriging_W]; 
XY_11_eta=sortrows(XY_eta,2);
XY_12_eta=sortrows(XY_11_eta,1);
Kriging_W_etaS=XY_12_eta(:,(Dim_x+1):(Dim_x+M));
%Krig_W_etaMAT(i,:)=Kriging_W_etaS;%
Krig_site_eta=XY_12_eta(:,1:Dim_x);

if mod(i,50)<=0
%     AA{1}
%     AA{2}
    figure
XX1=[]; XX2=[]; KW_U=[];
for jj=1:(mm+1)
     KW_U=[KW_U Kriging_W_etaS((1+(jj-1)*(mm+1)):(jj*(mm+1)))];
     XX1=[XX1 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),1)];
     XX2=[XX2 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),Dim_x)];    
end

imagesc(XX1(:),XX2(:),KW_U); colorbar; set(gca,'YDir','normal'); %caxis([0.95 1.65]); 

pause(.0001); % <- a time consuming OP  
fprintf(1,'iter %d\n',i);
end
end