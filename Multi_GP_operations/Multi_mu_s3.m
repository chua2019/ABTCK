function [Mu_Zt, Sigma, beta_X, prCurr]=Multi_mu_s3(Zt,AA, site,X,tau1,tau2)
% z(x,t)=eta(x,t)+delta(x)+e(x)
% p is the dimention of x in Kennedy-O'Hagan
% site1 are the observed sites with the unknown thetas
% site2 are the simulated sites of (x,t).
% total site=[site1; site2]
% careful how you define X = [H_1(y1) 0; rho H_1(D_2) H_2(D_2)] -- Usualy 
% A_eta parameters of eta 
% A_delta parameters of delta
% SIG_eta the variance of eta
% SIG_delta the variance of delta
n=size(Zt,1); q=size(Zt,2);  q_X=size(X,2); 
% X=X_un;
% for kk=1:q_X
% X(:,kk)=(X_un(:,kk)-mean(X_un(:,kk)))./std(X_un(:,kk));
% end
m=size(X,2);

C_prop_S =Gauss_Cov2(1,AA, site,site,tau2);

Q_S = chol(C_prop_S); 
% Regression mean

  X_2=zeros(n,q_X);
%       XX=X;%careful with the kron-product.      
%       X22 =(Q_S'\(XX)); X2 =(Q_S\X22);      
%       X_1= X2'; X_2=X_1; %reshape(X_1,nn,1);
  for j=1:q_X
      XX=X(:,j);%careful with the kron-product.      
      X22 =(Q_S'\(XX)); X2 =(Q_S\X22);      
      X_1= X2'; X_2(:,j)=X_1; %reshape(X_1,nn,1);
  end 
% when the functions are the same in some subregions we obser some
% numerical instabilities -- this can be solved better separately
% targetting the subregion
 beta_X=(chol(X'*X_2+eye(q_X)*tau1)\((chol(X'*X_2+eye(q_X)*tau1)'\X_2')))*Zt; %+eye(q_X)*0.1%beta_X=inv(X'*X_2)*X_2'*Zt;  %
 %beta_X=inv(X'*X_2+eye(q_X)*0.1)*X_2'*Zt;
 Mu_Zt=X*beta_X;
% Empirical value.

  Z=zeros(n,q); Y=(Zt-Mu_Zt); 
  for i=1:q
      YY=Y(:,i); %vec2mat(Y(:,i),n);
      Y22 =(Q_S'\YY);  Y2 = (Q_S\Y22);       
      Z(:,i) = Y2';
  end 
Sigma =(1/(n-m))*Y'*Z;
  C_prop_Y=Sigma; %C_prop_Y=Sigma+0.00005*eye(q);

  Q_Y = chol(C_prop_Y);
  logdet_S = 2 * sum(log(diag(Q_S)));
  logdet_Y = 2 * sum(log(diag(Q_Y)));
  bs_Y1=Q_Y'\(Y'); bs_Y =Q_Y\bs_Y1;
  %bs_Y=Q_Y'\(Y'); 
  % Likelihood 
  prCurr=-(q/2)*logdet_S-(n/2)*logdet_Y-1/2*trace(bs_Y*Z);

% Covariances
% C_prop_S=Gauss_Cov(SIG,AA, site,site);
% Q_S = chol(C_prop_S); 
% % Regression mean 
% q_X=size(X,2); m = size(X,2);
% X_2=zeros(n,q_X);
% for j=1:q_X
%   XX=X(:,j);%careful with the kron-product.      
%   X22 =(Q_S'\(XX)); X2 =(Q_S\X22);      
%   X_2(:,j)=X2'; %reshape(X_1,nn,1);
% end   
% HAH=(X_2'*X); 
% Q_HAH=chol(HAH);
% 
% beta_X=inv(X'*X_2)*X_2'*Zt; 
% Mu_Zt=X*beta_X;
%  
% % log of Exp part
% Y_2=zeros(n,q); Y=(Zt-Mu_Zt);  
% for i=1:q
%   YY=Y(:,i); %vec2mat(Y(:,i),n);
%   Y22 =(Q_S'\YY);  Y2 = (Q_S\Y22);       
%   Y_2(:,i) = Y2';
% end 
% log_exp_part =(1/2)*Y'*Y_2; %   
% Sigma =(1/(n-m))*Y'*Y_2; % 
% 
%   logdet_S = 2*sum(log(diag(Q_S)));
%   logdet_HAH = 2*sum(log(diag(Q_HAH)));  
% % Likelihood
% prCurr=-(q/2)*logdet_S-(q/2)*logdet_HAH-log_exp_part;  
if (isreal(prCurr)==0)
  MALAKA1=0;
end

end
