function [Mu_Zt, Sigma, prCurr]=MV_ALLD_2_s(Zt,theta, AA,sigma2S, site,X_s,funname)
% compute the posterior of Alpha: the coef of the mean equation.
n=size(Zt,1); q=size(Zt,2);  nn=n;
% Covariances
R_phi_prop=Cov_func_ALLD(theta, AA, site, site,funname);
C_prop_S = R_phi_prop+sigma2S.*eye(nn);
% Empirical value.
  Q_S = chol(C_prop_S); 
% Regression mean
  X=X_s; q_X=size(X,2); m = size(X_s,2);
  X_2=zeros(n,q_X);
  for j=1:q_X
      XX=X(:,j);%careful with the kron-product.      
      X22 =(Q_S'\(XX)); X2 =(Q_S\X22);      
      X_1= X2'; X_2(:,j)=X_1; %reshape(X_1,nn,1);
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
Sigma =(1/(n-m))*Y'*Y_2;
C_prop_Y=Sigma; %C_prop_Y=Sigma+0.00005*eye(q);
Q_Y = chol(C_prop_Y);
Q_HAH=chol(X'*X_2);
% log_det at the proposed value.
  logdet_S = 2*sum(log(diag(Q_S)));
  logdet_HAH = 2*sum(log(diag(Q_HAH)));
  logdet_Y = 2*sum(log(diag(Q_Y)));
  % Likelihood 
  prCurr=-(q/2)*logdet_S-(q/2)*logdet_HAH-(n-m)/2*logdet_Y;  
  if (isreal(prCurr)==0)
      MALAKA1=0;
  end

end
