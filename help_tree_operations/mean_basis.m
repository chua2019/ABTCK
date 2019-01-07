

function [X_s]=mean_basis(nn,nn1,nn2,site_new, funname2)
switch funname2
case 'constant'    
    if (nn1==0)
    X1 = [ones(nn,1)];  X_s=[X1];
    else
    X1 = [ones(nn,1)]; X2=[ones(nn1,1); zeros(nn2,1)]; X_s=[X1  X2];
    %Dim_x=size(site_new,2); X1 = [ones(nn,1)]; X2=[[ones(nn1,1); zeros(nn2,1)], [site_new(1:nn1,:); zeros(nn2,Dim_x)]]; X_s=[X1  X2];
    end  
case 'linear'
    if (nn1 ==0)
    X1 = [ones(nn,1) site_new];  X_s=[X1];
    else
    X1 = [ones(nn,1) site_new]; X2=[ones(nn1,1); zeros(nn2,1)]; X_s=[X1  X2];
    %Dim_x=size(site_new,2); X1 = [ones(nn,1) site_new]; X2=[[ones(nn1,1); zeros(nn2,1)], [site_new(1:nn1,:); zeros(nn2,Dim_x)]]; X_s=[X1  X2];
    end
case 'quadratic'
    if (nn1 ==0)
    X1 = [ones(nn,1) site_new site_new.^2];  X_s=[X1];
    else
    X1 = [ones(nn,1) site_new site_new.^2]; X2=[ones(nn1,1); zeros(nn2,1)]; X_s=[X1  X2];
    end

end 