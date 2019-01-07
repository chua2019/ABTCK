
function [site1, prCurr]=Calib_Theta_s(Zt, prCurr, SIG_eta,A_eta,SIG_delta, A_delta,tau2_Z, tau2_Y, site1,site2,p,X,lim_all,i,funname)
% M-H algorithm for Diag of theta (Kennedy and O'Hagan 2001).
% A_eta refere to eta -- dimention (p+l)x(p+1)
Dim_x =size(site1,2); n_1 =size(site1,1);
site1_prop=site1;
th_bound.max=lim_all(i,2);
th_bound.min=lim_all(i,1);
site1_prop(:,i)=ones(n_1,1)*tnormrnd(site1(1 ,i),(th_bound.max-th_bound.min)/6,th_bound.min,th_bound.max);

q=size(Zt,2); n=size(Zt,1);
q_X=size(X,2); m = size(X,2);
C_prop_S=Cov_Calib2(SIG_eta,A_eta,SIG_delta, A_delta,tau2_Z,tau2_Y, site1_prop,site2,p,funname);

Q_S = chol(C_prop_S);
X_2=zeros(n,q_X);
  for j=1:q_X
      XX=X(:,j);%careful with the kron-product.      
      X22 =(Q_S'\(XX)); X2 =(Q_S\X22);      
      X_2(:,j)=X2';  %reshape(X_1,nn,1);
  end   
 beta_X=inv(X'*X_2)*X_2'*Zt; Mu_Z=X*beta_X;
 HAH=(X'*X_2); 
 Q_HAH=chol(HAH);  
Y_2=zeros(n,q); Y=(Zt-Mu_Z);
  for j=1:q
      YY=Y(:,j); %vec2mat(Y(:,i),n);
      Y22 =(Q_S'\YY);  Y2 = (Q_S\Y22);       
      %ZZ= Y2'; Z(:,i)=reshape(ZZ,n,1);
      Y_2(:,j) = Y2';
  end
    log_exp_part =(1/2)*Y'*Y_2; % 
% C_prop_Y = 1/(n-m)*Y'*Y_2;
% Q_Y = chol(C_prop_Y);  
%  logdet_Y = 2*sum(log(diag(Q_Y)));
% log_det at the proposed value.
  logdet_S = 2*sum(log(diag(Q_S)));
  logdet_HAH = 2*sum(log(diag(Q_HAH)));
% Likelihood 
prProp=-(q/2)*logdet_S-(q/2)*logdet_HAH-log_exp_part; 
if (isreal(prProp)==1)
% from phi_curr --> phi_prop.
log_pr_tran = -tnormlike(site1_prop(1,i),site1(1,i),(th_bound.max-th_bound.min)/6,th_bound.min,th_bound.max);
% from phi_prop --> phi_curr.
log_pr_rtran = -tnormlike(site1(1,i),site1_prop(1,i),(th_bound.max-th_bound.min)/6,th_bound.min,th_bound.max);
% prior Beta
% prior parameters a and b betapdf(x,a,b)

% log_prior_Curr= log(betapdf((site1(1,i)-th_bound.min)/(th_bound.max-th_bound.min),a,b));
% log_prior_Prop=log(betapdf((site1_prop(1,i)-th_bound.min)/((th_bound.max-th_bound.min)),a,b));

% MH ratio
MH_ratio =(prProp-prCurr)+(log_pr_rtran-log_pr_tran);%+(log_prior_Prop-log_prior_Curr);
u = log(rand);
    if (u<MH_ratio)
    % accept
        site1 = site1_prop; prCurr = prProp;
    end
end
end

