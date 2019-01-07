function Cmat=Cov_func_ALLD(theta, SIGMA, site1,site2,funname)
switch funname
    case 'gauss'
        if (size(site1,1)==0 && size(site2,1)==0)
        Cmat=[];
        else
        Cmat=Gaussi_Cov_ALLD(theta, SIGMA, site1,site2);
        end
    case 'matern'
        t1=1; t2=1;
        Cmat=R_phi_ALLD(theta, SIGMA, site1,site2, t1, t2);
end
end

