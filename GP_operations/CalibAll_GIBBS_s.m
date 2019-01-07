function [S_MAT1, S_MAT2,sigma_2, TAU_2,prCurr_MAT,mean_Z1_MAT,beta_X,SIGMA_MAT,SIGm_eta,SIGm_delta, Theta_v, Delta_v]= ...
    CalibAll_GIBBS_s(Zt,site1,site2,SIG_eta,A_eta,SIG_delta, A_delta,tau2_Z,tau2_Y,Delta1,lim_all,BB,p,funname,funname2)
    Dim_x=size(site2,2); 
    n1=size(site1,1);
    q=2; % change this
%BB=4000;
S_MAT1=zeros(BB,Dim_x,Dim_x);
S_MAT2=zeros(BB,p,p);
Theta_v=zeros(BB,(Dim_x-p));

SIGm_eta= zeros(BB,1);
SIGm_delta=  zeros(BB,1);
switch funname2
    case 'constant'
Delta_v= zeros(BB,2);
    case 'linear'
Delta_v= zeros(BB,4);
end
sigma_2=zeros(BB,1);TAU_2=zeros(BB,1); theta_beta=zeros(BB,1);
mean_Z1_MAT=zeros(BB,q); SIGMA_MAT=zeros(BB,q,q); prCurr_MAT=zeros(BB,1);
for i=1:BB
    [A_eta,A_delta,theta_v, site1, SIG_eta, SIG_delta, tau2_Z, tau2_Y, prCurr, mu,SIGMA, beta_X]= ...
        CalibAll_MH_s(Zt,site1,site2,SIG_eta,A_eta,SIG_delta, A_delta,tau2_Z,tau2_Y,Delta1,lim_all,p,funname,funname2);
                 S_MAT1(i,:,:) = A_eta;
                 S_MAT2(i,:,:) = A_delta;
                 Theta_v(i,:) = theta_v;
                 SIGMA_MAT(i,:,:) =SIGMA;
                 mean_Z1_MAT(i,:)=mu(1,:); 
                 SIGm_eta(i,:)=SIG_eta;
                 SIGm_delta(i,:)=SIG_delta;                
                 Delta_v(i,:)= beta_X;            
                 prCurr_MAT(i)=prCurr;
                 TAU_2(i,1)=tau2_Y;
                 TAU_2(i,2)=tau2_Z;
         %fprintf(1,'iter %d\n',i);
end
end
