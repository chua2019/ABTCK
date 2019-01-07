function [R_phi_prop]=R_phi_ALLD(theta1, Lambda1, site1, site2,t1, t2)
n = size(site1,1);m = size(site2,1);
Dim=size(site1,2);
inv_S=inv(Lambda1); 
DD=0;
for k=1:Dim
SITE_M1n{k}=kron(site1(:,k),ones(1,m));
SITE_M1m{k}=kron(site2(:,k),ones(1,n));
DSITE_1{k}=SITE_M1n{k}-SITE_M1m{k}';
D2{k} = inv_S(k,k)*DSITE_1{k}.^2;
DD=(D2{k}+DD);
end
Aniso_Dist=sqrt(DD);
% Aniso_Dist1 = pdist2(X,Y,'mahalanobis',SIGMA);  % mahalanobis distance
% %Aniso_Dist = squareform(Aniso_Dist1);
% Aniso_Dist=sqrt(Aniso_Dist1);
alld.Matd = Aniso_Dist; 
alld.Matdt = sqrt((t2-t1)^2);
R_phi_prop=Matern_ST3(theta1,alld);
end
