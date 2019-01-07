function [Kriging_W, Var_krig]=Krig_Multi(Zt,beta,rho,sigma2,AA,site,site_e,S,tau2)
% Kriging with nuget effect.
% [site1; site2] real and computer code data.  
% site = the acctual sites where we have data.
% site_e is the grided locations where we do prediction.
% ***** Usualy l_1 and l_2 equal with a*****%
% nn is a cell with the number of points in the Element (region)
pseudo_sigma=cell(S,1);
for i=1:S
pseudo_sigma{i}=1;
end
Kriging_W=[];Var_krig=[];
         if (size(site_e,1)>0) % && size(site1,1)>0 && size(site2,1)>0
            % for j=1:Dim_y
         COV_M=V_Multi(sigma2,AA,rho,site,S,tau2);        
         c=t_star(sigma2,AA,rho,site_e,site,S,tau2);                
%          I=ones(nn,1);
%          Ind=(COV_M\I)';
%          m_error=(1-Ind*c)/(Ind*I);
%          cc=(c+I*m_error);

         lambda=(COV_M\c)'; %inv_C=inv(COV_M);lambda=c'*inv_C; %kron(eye(Dim),c')*inv_C; 

         FF_X=F_sBasis(site,rho, S);
         ZT{1}=Zt{1};
         for t=2:S
         ZT{t}=[ZT{t-1}; Zt{t}];
         end
         f_basis=F_primeX(site_e,rho, S);
         alpha{1}=beta{1};
         for t=2:S
         alpha{t}=[alpha{t-1}; beta{t}];
         end
         mu= FF_X*alpha{S};
         mu1=f_basis*alpha{S}; 
         %mu1(1,:)
         Kriging_W= mu1+lambda*(ZT{S}-mu); % 
         Var_Zs1=Var_Zs(sigma2,AA,rho,site_e,S,tau2);
         Var_krig = diag((Var_Zs1-lambda*c));%kron(diag((C_o-lambda*c)),(SIG_eta+SIG_delta))+m_error';        
         end
   
end
