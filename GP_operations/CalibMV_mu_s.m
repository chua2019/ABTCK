function [Mu_Zt, Sigma, beta_X, prCurr]=CalibMV_mu_s(Zt,SIG_eta,A_eta,SIG_delta, A_delta,tau2_Z, tau2_Y, site1,site2,p,X,funname)
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
n=size(Zt,1); q=size(Zt,2); 
% Covariances
C_prop_S=Cov_Calib2(SIG_eta,A_eta,SIG_delta, A_delta,tau2_Z, tau2_Y, site1,site2,p,funname);
Q_S = chol(C_prop_S); 
% Regression mean 
q_X=size(X,2); m = size(X,2);
X_2=zeros(n,q_X);
for j=1:q_X
  XX=X(:,j);%careful with the kron-product.      
  X22 =(Q_S'\(XX)); X2 =(Q_S\X22);      
  X_2(:,j)=X2'; %reshape(X_1,nn,1);
end   
HAH=(X_2'*X); 
Q_HAH=chol(HAH);

beta_X=inv(X'*X_2)*X_2'*Zt; 
Mu_Zt=X*beta_X;
 
% log of Exp part
Y_2=zeros(n,q); Y=(Zt-Mu_Zt);  
for i=1:q
  YY=Y(:,i); %vec2mat(Y(:,i),n);
  Y22 =(Q_S'\YY);  Y2 = (Q_S\Y22);       
  Y_2(:,i) = Y2';
end 
log_exp_part =(1/2)*Y'*Y_2; %   
Sigma =(1/(n-m))*Y'*Y_2; % 
%C_prop_Y=Sigma; 
%Q_Y = chol(C_prop_Y);
% logdet_Y = 2*sum(log(diag(Q_Y)));
% log_det at the proposed value.
  logdet_S = 2*sum(log(diag(Q_S)));
  logdet_HAH = 2*sum(log(diag(Q_HAH)));  
% Likelihood
prCurr=-(q/2)*logdet_S-(q/2)*logdet_HAH-log_exp_part;  
if (isreal(prCurr)==0)
  MALAKA1=0;
end

end

