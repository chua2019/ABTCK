function [A_curr,A_curr1, SIG_2, SIG1_2,tau2_Z, tau2_Y, prCurr, mu,SIGMA, beta_X]= CalibAll_MH_GP(Zt,site1,site2,SIG_2,A_curr,SIG1_2, A_curr1,tau2_Z,tau2_Y,Delta1,p,funname,funname2)
% Separate the 2D area into smaller areas in the y direc.   
Dim_x=size(site1,2); Dim_y=size(Zt,2); 
nn1=size(site1,1); nn2= size(site2,1);
nn =size(Zt,1);
site=[site1; site2];
% Mean variance.

X=mean_basis(nn,nn1,nn2,site, funname2);
[mu, SIGMA, beta_X, prCurr]=CalibMV_mu_s(Zt,SIG_2,A_curr,SIG1_2, A_curr1,tau2_Z,tau2_Y, site1,site2,p,X,funname);

% spatial correlation eta
for k1=1:Dim_x              
    [A_curr, prCurr]=Calib_Diag_s(Zt, prCurr, SIG_2,A_curr,SIG1_2, A_curr1,tau2_Z, tau2_Y, site1,site2,p,X,k1,Delta1.diag, funname);               
end
% varianc eta
[SIG_2, prCurr]=Calib_Var_s(Zt, prCurr, SIG_2,A_curr,SIG1_2, A_curr1,tau2_Z, tau2_Y, site1,site2,p,X,Delta1.var,funname);
% spatial correlation delta
for k2=1:p
     [A_curr1, prCurr]=Calib_Diag_sDelta(Zt, prCurr, SIG_2,A_curr,SIG1_2, A_curr1,tau2_Z, tau2_Y, site1,site2,p,X,k2,Delta1.diag,funname);            
end
% variance delta 
[SIG1_2, prCurr]=Calib_Var_sDelta(Zt, prCurr, SIG_2,A_curr,SIG1_2, A_curr1,tau2_Z, tau2_Y, site1,site2,p,X,Delta1.var,funname);
% calibration parameter
% theta_v = zeros((Dim_x-p),1);
%for i=(p+1):Dim_x
%[site1, prCurr]=Calib_Theta_s(Zt, prCurr, SIG_2,A_curr,SIG1_2, A_curr1,tau2_Z, tau2_Y, site1,site2,p,X,th_bound{(i-p)},i,funname);
%theta_v(i-p) =  site1(1,i); %theta_v
%end
% Nugget
[tau2_Z, prCurr]=Calib_Nugget_obs(Zt, prCurr, SIG_2,A_curr,SIG1_2, A_curr1,tau2_Z,tau2_Y, site1,site2,p,X,Delta1.nugget,funname); 
[tau2_Y, prCurr]=Calib_Nugget_cc(Zt, prCurr, SIG_2,A_curr,SIG1_2, A_curr1,tau2_Z,tau2_Y, site1,site2,p,X,Delta1.nugget,funname); 
end