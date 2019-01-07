function [X1]=mean_basis_eta(site_new2, funname2)
nn2=size(site_new2,1);
switch funname2
case 'constant'    
    X1 = [ones(nn2,1)]; 
case 'linear'
    X1 = [ones(nn2,1) site_new2]; 
case 'quadratic'
    X1 = [ones(nn2,1) site_new2 site_new2.^2];  
end 