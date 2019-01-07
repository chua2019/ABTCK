function [Sigma_eta_hat] =sigma_eta_hat(Zt,A_eta, tau2_Y, site,X,funname)
n=size(Zt,1); q=size(Zt,2); 
% Covariances
theta.sigma_2=1; theta.nu=1; theta.beta=1;
R_phi_prop=Cov_func_ALLD(theta, A_eta, site, site,funname);
C_prop_S = R_phi_prop+tau2_Y.*eye(n);
% Empirical value.
  Q_S = chol(C_prop_S); 
% Regression mean
 q_X=size(X,2); m = size(X,2);
  X_2=zeros(n,q_X);
  for j=1:q_X
      XX=X(:,j);%careful with the kron-product.      
      X22 =(Q_S'\(XX)); X2 =(Q_S\X22);      
      X_2(:,j)=X2'; %reshape(X_1,nn,1);
  end 
 beta_X=inv(X'*X_2)*X_2'*Zt; 
 Mu_Zt=X*beta_X;
% Regression Variance 
  Y_2=zeros(n,q); Y=(Zt-Mu_Zt);  
  for i=1:q
      YY=Y(:,i); %vec2mat(Y(:,i),n);
      Y22 =(Q_S'\YY);  Y2 = (Q_S\Y22);       
      Y_2(:,i) = Y2';
  end 
  
Sigma_eta_hat =(1/(n-m))*Y'*Y_2; % 

end










