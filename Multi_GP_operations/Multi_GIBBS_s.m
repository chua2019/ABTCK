function [S_AA,prCurr_MAT,BETA_v,SIG_MAT]= Multi_GIBBS_s(Zt,site,AA,S,Delta,BB,tau1,tau2)
    Dim_x=size(site{1},2); 
S_AA=zeros(S,BB,Dim_x,Dim_x);
SIG_MAT=  zeros(S,BB);
prCurr_MAT=zeros(S,BB,1);

%    case 'constant'
BETA_v= zeros(S,BB,2);
for i=1:BB
[AA, SIG, prCurr, mu, beta_X]= Multi_MH_s(Zt,site,AA,S,Delta,tau1,tau2);
%beta_X{2}
beta_X{1}=[beta_X{1}; 0];
for t=1:S;
     S_AA(t,i,:,:) = AA{t};
     SIG_MAT(t,i) =SIG{t}; 
     BETA_v(t,i,:)= beta_X{t};            
     prCurr_MAT(t,i)=prCurr{t};
end

%fprintf(1,'iter %d\n',i);
end
end
