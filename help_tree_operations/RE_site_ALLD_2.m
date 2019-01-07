function [site_new, Zt_new, nn]= RE_site_ALLD_2(site1,Zt,Lim_ALL,a)

% Separate the 2D area into smaller areas in the y direc. 
% b is the chosen dimention to split.
n=size(site1(:,1),1);
for i=1:a
    L1=Lim_ALL{i}(:,1)'; L2=Lim_ALL{i}(:,2)';  
    k=1;
    if size(site1,2)==2
        for j=1:n
           if ((site1(j,1)>=L1(1))  && (site1(j,1)<L2(1)) && (site1(j,2)>=L1(2))&& (site1(j,2)<L2(2)))  
            site_new{i}(k,:) = site1(j,:);
            Zt_new{i}(k,:)=Zt(j,:);
            k=k+1;
           end
        end
    elseif size(site1,2)==1
        for j=1:n
           if ((site1(j,1)>=L1(1))  && (site1(j,1)<L2(1)))  
            site_new{i}(k,:) = site1(j,:);
            Zt_new{i}(k,:)=Zt(j,:);
            k=k+1;
           end
        end
    end
    nn{i}=k-1;
end
end
