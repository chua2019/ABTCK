function [Var_Zs1]=Var_Zs(sigma2,AA,rho,x,S,tau2)
       Var_Zs1=sigma2{S}.*Gauss_Cov2(1, AA{S}, x, x,tau2);
            for k=1:(S-1)           
%                 rho_com{k}=1;
%                 for j=k:(S-1)
%                     rho_com{k}=rho_com{k}*rho{j}^2;
%                 end
                rho_com{k}=1;
                for j=(k):(S-1)
                    rho_com{k}=rho_com{k}*rho{j+1}^2;
                end
             Var_Zs1=Var_Zs1+(rho_com{k}).*sigma2{k}.*Gauss_Cov2(1, AA{k}, x, x,tau2);
            end
end