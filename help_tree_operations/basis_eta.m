function [L_eta]=basis_eta(site2, funname2)
q=size(site2,2);
switch funname2
case 'constant' 
L_eta=1;
case 'linear'
L_eta=1+q;  
case 'quadratic'
L_eta=1+2*q;
end
end