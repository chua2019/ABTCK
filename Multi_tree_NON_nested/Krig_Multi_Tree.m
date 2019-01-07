function [Kriging_W,site_W]=Krig_Multi_Tree(Zt,beta,rho,sigma2,AA,site,site_e,a_new,S,tau2)

% Kriging with nuget effect.
% site_new = the acctual sites where we have data.
% site_new_e is the grided locations where we do prediction.
% ***** Usualy l_1 and l_2 equal with a*****%
% nn is a cell with the number of points in the Element (region)

%Zt_prev{1}=[]; Zt_prevE{1}=[];
% Covariances for site;
for k=1:a_new
    %rho{k}
    %AA{k}{1}
    %AA{k}{2}
    % beta{k}
%     Zt{k}{1}
[Kriging_Wi{k}, Var_krig{k}]=Krig_Multi(Zt{k},beta{k},rho{k},sigma2{k},AA{k},site{k},site_e{k},S,tau2);
% Zt_prev{k+1}=[Zt_prev{k}  ZT(:,k)];
% Zt_prevE{k+1}=[Zt_prevE{k}  Kriging_W{k}];
end

         [Kriging_W]=[Kriging_Wi{1}];
         [site_W]=[site_e{1}];
         for i=2:a_new
             [Kriging_W]=[Kriging_W; Kriging_Wi{i}];
             [site_W]=[site_W ; site_e{i}];
         end   

end
