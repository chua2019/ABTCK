function [Cmat]=CrosCov_Calib(SIG_2,AA,SIG1_2, AA1, site1,site2,site_e,p,funname)
% Following Higdon 2005 + Conti 2010 --- Separable model 
% site1 represent the n observations
% site2 represent the m computer experiments (models)
% p represents the  input parametes
% l represents the calibration parametes
% tau2_Y represent the observation error

site=[site1; site2];
n_m=size(site,1); n = size(site1,1); m = size(site2,1);

theta.beta=1; theta.betat=0; theta.nu=1.5; theta.sigma_2=1; theta1=theta;

Sigma_eta2= SIG_2.*(Cov_func_ALLD(theta, AA,site2,site_e,funname));% +tau2_Y.*eye(n_m)%

SYd_11=[];
if (size(site1,1)>0)
SYd_11=SIG1_2.*(Cov_func_ALLD(theta1, AA1, site1(:,1:p),site_e(:,1:p),funname));%+tau2_Y.*eye(n);
end
SYd_21=SIG_2.*(Cov_func_ALLD(theta, AA,site1,site_e, funname));
Cmat=[(SYd_11+SYd_21); Sigma_eta2];
end





