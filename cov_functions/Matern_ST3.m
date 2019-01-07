function [Cmat]=Matern_ST3(theta,alld)
sigma_2= theta.sigma_2; nu=theta.nu;
beta=theta.beta; %betat=theta.betat; 
  sh2 = (alld.Matd)./beta;%sqrt((alld.Matd./beta).^2+(alld.Matdt./betat).^2);    
  m3=(sh2).^nu.*besselk(nu,sh2)./(2.^(nu-1).*gamma(nu));%((sh2).^nu).*besselk(nu,sh2)./(2.^(nu-1).*gamma(nu)); %
  m3(sh2==0) =1;  
  result =sigma_2.*m3; 
  Cmat=result;  
end
