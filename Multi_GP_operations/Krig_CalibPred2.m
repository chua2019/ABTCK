function [Kriging_W, Var_krig]=Krig_CalibPred2(Zt,alpha,SIG_eta,A_eta,SIG_delta, A_delta,tau2_Z, tau2_Y, site1, site2,site_e,p,funname,funname2)
% Kriging with nuget effect.
% [site1; site2] real and computer code data.  
% site = the acctual sites where we have data.
% site_e is the grided locations where we do prediction.
% ***** Usualy l_1 and l_2 equal with a*****%
% nn is a cell with the number of points in the Element (region)

site = [site1; site2]; nn=size(site,1); nn1=size(site1,1); nn2=size(site2,1);

Kriging_W=[];
         if (size(site_e,1)>0) % && size(site1,1)>0 && size(site2,1)>0
            % for j=1:Dim_y
         COV_M=Cov_Calib2(SIG_eta,A_eta,SIG_delta, A_delta,tau2_Z,tau2_Y, site1,site2,p,funname);        
         c=CrosCov_Calib(SIG_eta,A_eta,SIG_delta, A_delta, site1,site2,site_e,p,funname);        
         site_e_pseudo=[];
         C_o=Cov_Calib2(SIG_eta,A_eta,SIG_delta, A_delta,tau2_Z,tau2_Y, site_e,site_e_pseudo,p,funname); 
         
         I=ones(nn,1);
         Ind=(COV_M\I)';
         m_error=(1-Ind*c)/(Ind*I);
         cc=(c+I*m_error);
         lambda=(COV_M\cc)'; %inv_C=inv(COV_M);lambda=c'*inv_C; %kron(eye(Dim),c')*inv_C;      
         X =mean_basis(nn,nn1,nn2,site, funname2);
    switch funname2
        case 'constant'
            %X = ones(size(site,1),2);
            if (size(site1,1)==0)
                X_e = ones(size(site_e,1),1);  
            else
                X_e = ones(size(site_e,1),2);  
            end
        case 'linear'
            %X = [ones(size(site,1),1) site  ones(size(site,1),1)];
            if (size(site1,1)==0)
                X_e = [ones(size(site_e,1),1) site_e];
            else
                X_e = [ones(size(site_e,1),1) site_e ones(size(site_e,1),1)]; 
            end
        case 'quadratic'
            %X = [ones(size(site,1),1) site site.^2 ones(size(site,1),1)];
            X_e = [ones(size(site_e,1),1) site_e site_e.^2 ones(size(site_e,1),1)];
    end
         mu1= X_e*alpha;
         mu= X*alpha;
         Kriging_W=mu1+lambda*(Zt-mu);       

         Var_krig = diag((C_o-lambda*c))+m_error';%kron(diag((C_o-lambda*c)),(SIG_eta+SIG_delta))+m_error';
         
         end
   
end
