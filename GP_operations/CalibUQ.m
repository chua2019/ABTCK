function [Kriging_W]=CalibUQ(Zt,alpha,SIG_2,A_curr,SIG1_2, A_curr1,tau2_Y, site1, site2,mu_e, var_e,p,funname,funname2)
% True only if the covaraince function is exponential.

site = [site1; site2];
         if size(site_e,1)>0
            % for j=1:Dim_y
         COV_M=Cov_Calib(SIG_2,A_curr,SIG1_2, A_curr1,tau2_Y, site1,site2,p,funname);        
         c=CrosCov_Calib(SIG_2,A_curr,SIG1_2, A_curr1, site1,site2,site_e,p,funname);        
         lambda=(COV_M\c)'; %inv_C=inv(COV_M);lambda=c'*inv_C; %kron(eye(Dim),c')*inv_C;                   
    switch funname2
        case 'constant'
            X = ones(size(site,1),2);
            X_e = ones(size(site_e,1),2);          
        case 'linear'
            X = [ones(size(site,1),1) siteones(size(site,1),1)];
            X_e = [ones(size(site_e,1),1) site_e ones(size(site_e,1),1)];
        case 'quadratic'
            X = [ones(size(site,1),1) site site.^2 ones(size(site,1),1)];
            X_e = [ones(size(site_e,1),1) site_e site_e.^2 ones(size(site_e,1),1)];
    end
         mu1= X_e*alpha;
         mu= X*alpha;
         Kriging_W=mu1+lambda*(Zt-mu);       
         end
   
end
