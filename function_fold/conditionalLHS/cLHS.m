function [ipick,xsam,xobs,oL,isam]=cLHS(nsample,x,xobj,vobj,niter,w1,w2,w3)
% [ipick,xsam,xobs,oL,isam]=cLHS(nsample,x,xobj,vobj,niter,w1,w2,w3)
% conditioned Latin hypercube sampling 
% Input:
% nsample   = no of sample
% x = continuous data
% xobj = class/object data
% vobj = values of the object
% optional:
%   niter = no iterations
%   w1 = weight for continuous data
%   w2 = weight for corelation among data
%   w3 = weight for object data
%
%Output:
% ipick = index of the selected samples (from data matrix x)
% xsam = sampled continuous data
% xobs = sampled class data
% oL = objective function at each iteration
% isam = matrix containing the count for each strata
%
% (c) 2005 Australian Center for Precision Agriculture

% parameters
if(nargin==4),
    niter=10000;
end
if(nargin<6),
    w1=1;
    w2=1;
    w3=1;
end;


% annealing schedule
temp=1;
tdecrease=0.95;
metrop=1;

[ndat,nvar]=size(x);

% the edge of the strata
step=100/nsample;
Pl=[0:step:100]';
pc=([0.5:ndat]-0.5)./ndat;
for j=1:nvar,
    xedge(:,j)=prctile(x(:,j),[0:step:100]');
    [Y,I] = sort(x(:,j));
    xquant(I,j)=pc';
end
% target value
ibest=ones(nsample,nvar);

% data correlation
corr=corrcoef(x);

% for object/class variable
if(length(xobj)==0),
    nobj=0;
    cobj=0;
    w3=0;
    xobs=[];
else
    nobj=length(vobj);
    for j=1:nobj
        cobj(j)=sum(xobj==vobj(j));
    end
    cobj=cobj./ndat*nsample;    % target value
end

% initialise, pick randomly
nunsam=ndat-nsample;
idat=randperm(ndat);
ipick=idat(1:nsample);
iunsam=idat(nsample+1:ndat);
xsam=x(ipick,:);
if(nobj>0), xobs=xobj(ipick,:); end;

% objective function
[obj,isam,dif,iobj]=LHSobj(nsample,nvar,xsam,xobs,xedge,ibest,w1,w2,corr,w3,vobj,cobj);
    
tic
for il=1:niter,
   % il
    idx=randperm(nunsam);
    iunsam=iunsam(idx);

    objsave=obj;
    isave=ipick;
    iusave=iunsam;
    dsave=dif;
    
    if (rand<0.5), 
         % pick a random sample & swap with reservoir
        iw=randperm(nsample);        
        ich=iunsam(1);
        jch=ipick(iw(1));
        ipick(iw(1))=ich;
        iunsam(1)=jch;
        xsam(iw(1),:)=x(ich,:);
        if(nobj>0), xobs(iw(1),:)=xobj(ich,:); end;
    
    else
        % remove the worse sampled & resample
        worse=max(dif);
        iworse=find(dif==worse);
        nworse=length(iworse);

        % swap with reservoir
        jch=ipick(iworse);
        kch=iunsam(:,1:nworse);
        ipick(iworse)=kch;
        iunsam(:,1:nworse)=jch;    
        xsam(iworse,:)=x(kch,:);
        if(nobj>0), xobs(iworse,:)=xobj(kch,:); end;

    end

%   calc obj
    [obj,isam,dif,iobj]=LHSobj(nsample,nvar,xsam,xobs,xedge,ibest,w1,w2,corr,w3,vobj,cobj);
    
%   compare with previous iterations    
    de=obj-objsave;
    metrop=exp(-de/temp)+rand*temp;

    if(obj==0), 
        oL(il)=obj;
        break; 
    end;
        
    if or(de<=0, rand<metrop), % accept change
    else
        % revert, no changes
        ipick=isave;
        iunsam=iusave;
        xsam=x(ipick,:);
        if(nobj>0), xobs=xobj(ipick,:); end
        obj=objsave;
        dif=dsave;
        
    end
    
    oL(il)=obj;
    if mod(il,10)==0,
        temp=temp*tdecrease;
    end;
    
end

% calc the final obj function
[obj,isam,dif,iobj]=LHSobj(nsample,nvar,xsam,xobs,xedge,ibest,w1,w2,corr,w3,vobj,cobj);
ipick=ipick';

toc