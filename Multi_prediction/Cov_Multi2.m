
function [Cmat]=Cov_Multi2(sigma2,AA,rho, site)
% Following Kaennedy & O'hagan 
% site1 represent the n observations
% site2 represent the m computer experiments (models)
% p represents the  input parametes
% l represents the calibration parametes
% tau2_Y represent the observation error

SYd_11=V_tt(sigma2,AA,rho,site{1},site{1},1);
SYd_2=V_tt(sigma2,AA,rho,site{2},site{2},2);
SYd_12=(rho{2}^2).*V_tt(sigma2,AA,rho,site{1},site{2},2);

Sigma_Y=[SYd_11 SYd_12;SYd_12' SYd_2];
Cmat=Sigma_eta+Sigma_Y;
end
