            
function [COMMM]=V_tt(sigma2,AA,rho,site1,site2,i,tau2)
COMMM=sigma2{i}.*Gauss_Cov2(1, AA{i}, site1,site2,tau2);
            %rho_com=1;
            if i>1
            for k=1:(i-1)
                rho_com{k}=1;
                for kk1=k:(i-1)
                rho_com{k}=rho_com{k}*(rho{kk1+1}^2);
                end
             COMMM=COMMM+(rho_com{k}).*sigma2{k}.*Gauss_Cov2(1, AA{k}, site1,site2,tau2);
            end
            end
end