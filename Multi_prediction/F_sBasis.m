function [F_basis]=F_sBasis(site,rho, S)
       if S==1 
           nn1=size(site{1},1);
          X1 = [ones(nn1,1)];  F_basis=[X1];
       end

       if S==2 
           nn1=size(site{1},1);
           nn2=size(site{2},1);
          X1 = [ones(nn1,1); rho{2}*ones(nn2,1)]; X2=[zeros(nn1,1); ones(nn2,1)]; F_basis=[X1  X2];
       end
       if S==3  
           nn1=size(site{1},1);
           nn2=size(site{2},1);
           nn3=size(site{3},1);
          X1 = [ones(nn1,1); rho{2}*ones(nn2,1); rho{2}*rho{3}*ones(nn3,1)]; 
          X2=[zeros(nn1,1); ones(nn2,1); rho{3}*ones(nn3,1)]; 
          X3=[zeros(nn1,1); zeros(nn2,1); rho{3}*ones(nn3,1)];
          F_basis=[X1  X2 X3];
       end
           
end