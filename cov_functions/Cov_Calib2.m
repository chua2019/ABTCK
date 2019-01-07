
function [Cmat]=Cov_Calib2(SIG_eta,A_eta,SIG_delta, A_delta,tau2_Z,tau2_Y, site1,site2,p,funname)
% Following Kaennedy & O'hagan 
% site1 represent the n observations
% site2 represent the m computer experiments (models)
% p represents the  input parametes
% l represents the calibration parametes
% tau2_Y represent the observation error

site=[site1; site2];
n_m=size(site,1); n = size(site1,1); m = size(site2,1);

theta.beta=1; theta.betat=0; theta.nu=1.5; theta.sigma_2=1; theta1=theta; %this is not important
Sigma_eta=SIG_eta.*(Cov_func_ALLD(theta, A_eta, site,site,funname)+blkdiag(0.*eye(n),tau2_Y.*eye(m))); %
if (n==0)
    SYd_11=[];
else
    SYd_11=SIG_delta.*(Cov_func_ALLD(theta1, A_delta, site1(:,1:p),site1(:,1:p),funname))+tau2_Z.*eye(n);
end
SYd_2=zeros(m); SYd_12=zeros(n,m);

Sigma_Y=[SYd_11 SYd_12;SYd_12' SYd_2];
Cmat=Sigma_eta+Sigma_Y;
end
