function [COV]=Pred_covMulti(sigma2,AA,rho,site,S)
% Compute the covariance function for the whole region 
% where l_1=1:a_new
% Use this to estimate papameters

% theta: parameter of the GP with variance
% Lambda: Parameter of the anisotropoic matrix for mahalanobis distance
% site_new: observation data  sites
% funname: 'gauss' or 'matern'

a=S;

for i=1:S
    for j=1:S 
        if (i==j)
%             COMMM=sigma2{i}.*Gauss_Cov(1, AA{i}, site{i},site{i});
%             rho_com=1;
%             if i>1
%             for k=1:(i-1)
%             rho_com=rho_com*rho{i-k}^2;
%              COMMM=COMMM+(rho_com).*sigma2{i-k}.*Gauss_Cov(1, AA{i-k}, site{i},site{i});
%             end
%             end
            CCC{i,j}=V_tt(sigma2,AA,rho,site{i},site{j},i);
        else 
            if i>j
            for k=j:(i-1)
             rho_com=rho_com.*rho{k}^2;
            end
            CCC{i,j}= rho_com.*V_tt(sigma2,AA,rho,site{i},site{j},j);
            else
            for k=i:(j-1)
             rho_com=rho_com.*rho{k}^2;
            end
            CCC{i,j}= rho_com.*V_tt(sigma2,AA,rho,site{i},site{j},i);
            end
        end          
    end
end

for j=1:a
    Cov_M1{1,j}=CCC{j,1};
    if a>1
for i=2:a
Cov_M1{i,j}=[Cov_M1{i-1,j} CCC{j,i}];
end
    end
end

Cov_M11{1}=Cov_M1{a,1};
if a>1
for j=2:a
Cov_M11{j}=[Cov_M11{j-1}; Cov_M1{a,j}];
end
end
COV=Cov_M11{a};
end
