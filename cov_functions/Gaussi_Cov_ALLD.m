function [Cmat]=Gaussi_Cov_ALLD(theta, Lambda, site1,site2)
n = size(site1,1); m = size(site2,1);
Dim=size(site1,2);
SIGMA=Lambda;%SIGMA=Lambda'*Lambda;%
inv_S=inv(SIGMA);
DD=0;
for k=1:Dim
SITE_M1n{k}=kron(site1(:,k),ones(1,m));
SITE_M1m{k}=kron(site2(:,k),ones(1,n));
DSITE_1{k}=SITE_M1n{k}-SITE_M1m{k}';
D2{k} = inv_S(k,k)*DSITE_1{k}.^2;
DD=(D2{k}+DD);
end
sigma_2=theta.sigma_2;
beta=theta.beta;
Q=(DD)./(beta^2);
m3=exp(-Q); 
  result =sigma_2.*m3; 
  Cmat=result;  
end
