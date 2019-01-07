function [f_basis]=F_primeX(site,rho, S)
       if S==1 
           rho_com=1;
           nn1=size(site,1);
           
          X1 = [ones(nn1,1)];  f_basis=[X1];
       end
       if S==2 
           rho_com=1;
           nn1=size(site,1);

%         for j=1:(S-1)  % carefull here
%          rho_com=rho_com.*rho{j};
%         end
        for j=1:(S-1)
         rho_com=rho_com.*rho{j+1};
        end
           
          X1 = [rho_com*ones(nn1,1)]; X2=[ones(nn1,1)]; f_basis=[X1  X2];
       end
       if S==3  
           nn1=size(site,1);
           rho_com1=1;
%         for j=2:(S-1)
%          rho_com1=rho_com1.*rho{j};
%         end
%         rho_com2=1;
%         for j=2:(S-1)
%          rho_com2=rho_com2.*rho{j};
%         end
        for j=1:(S-1)
         rho_com1=rho_com1.*rho{j+1};
        end
        rho_com2=1;
        for j=2:(S-1)
         rho_com2=rho_com2.*rho{j+1};
        end
          X1 = rho_com1*ones(nn1,1); 
          X2=rho_com2*ones(nn1,1); 
          X3=ones(nn1,1);
          f_basis=[X1  X2 X3];
       end
           
end