function [obj,isam,dif,iobj]=LHSobj(nsample,nvar,xsam,xobs,xedge,ibest,w1,w2,corr,w3,vobj,cobj)
% objective function for LHS

% data in quantiles
if(w1==0),       
    dif=0;
else
    for j=1:nvar,
        hsam(:,j) = histc(xsam(:,j),xedge(:,j));    
        isam(:,j)= hsam(1:nsample,j);
    end
    dif=sum(abs(isam-ibest)');
end

% correlation
if(w2==0),
    dc=0;
else,
    csam=corrcoef(xsam);
    dc=sum(sum(abs(corr-csam)));
end

% object data
if(w3==0),
    do=0;
    iobj=[];
else,
    nobj=length(vobj);
    for j=1:nobj,
        iobj(j)=sum(xobs==vobj(j));
    end
    do=sum(abs(iobj-cobj));
end 

obj=w1.*sum(dif)+ w2.*dc + w3.*do;