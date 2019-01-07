
function [tau2_Y, prCurr]=Calib_Nugget_cc(Zt, prCurr, SIG_eta,A_eta,SIG_delta, A_delta,tau2_Z,tau2_Y, site1,site2,p,X,Delta,funname)
% M-H algorithm for thenugget.
% A_eta refere to eta -- dimention (p+l)x(p+1)
% A_delta refere to delta -- dimention pxp
Dim_x =size(site1,2);
up_bound = 6*10^(-1); low_bound=10^(-20);
    phi_curr = tau2_Y; 
    phi_prop1=logn_r(log(phi_curr), Delta,1); %lognrnd(log(phi_curr), Delta); 
    if (phi_prop1>low_bound && phi_prop1<up_bound)
            tau2_Y_prop= phi_prop1;  
        elseif phi_prop1>up_bound
            tau2_Y_prop= up_bound; 
        elseif phi_prop1<low_bound
            tau2_Y_prop= low_bound; 
    end             
q=size(Zt,2); n=size(Zt,1);
q_X=size(X,2); m=size(X,2);
C_prop_S=Cov_Calib2(SIG_eta,A_eta,SIG_delta, A_delta,tau2_Z,tau2_Y_prop, site1,site2,p,funname);
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
  for i=1:q
      YY=Y(:,i); %vec2mat(Y(:,i),n);
      Y22 =(Q_S'\YY);  Y2 = (Q_S\Y22);       
      %ZZ= Y2'; Z(:,i)=reshape(ZZ,n,1);
      Y_2(:,i) = Y2';
  end
    log_exp_part =(1/2)*Y'*Y_2; % 
% log_det at the proposed value.
  logdet_S = 2*sum(log(diag(Q_S)));
  logdet_HAH = 2*sum(log(diag(Q_HAH)));
% Likelihood 
prProp=-(q/2)*logdet_S-(q/2)*logdet_HAH-log_exp_part; 
if (isreal(prProp)==1)
% from phi_curr --> phi_prop.
log_pr_tran = -log(tau2_Y_prop)-(log(tau2_Y_prop)-log(tau2_Y))^2/(2*Delta^2);
% from phi_prop --> phi_curr.
log_pr_rtran = -log(tau2_Y)-(log(tau2_Y)-log(tau2_Y_prop))^2/(2*Delta^2);
% prior
% prior_curr=0;
% for ii=1:Dim_x
%     prior_curr=prior_curr+(-1)*log(1+tau2_Y^2);
% end
% prior_prop=0;
% for ii=1:Dim_x
%     prior_prop=prior_prop+(-1)*log(1+tau2_Y_prop^2);
% end
prior_curr=(-1)*(tau2_Y)/0.3;
prior_prop= (-1)*(tau2_Y_prop)/0.3;
% MH ratio
MH_ratio =(prProp-prCurr)+(log_pr_rtran-log_pr_tran)+(prior_prop-prior_curr);
u = log(rand);
    if (u<MH_ratio)
    % accept
        tau2_Y = tau2_Y_prop; prCurr = prProp;
    end
end
end