function [site_new, Zt_new, nn]= RE_site_ALLD(site1,Zt,Lim_ALL,a)

% Separate the 2D area into smaller areas in the y direc. 
% b is the chosen dimention to split.
n=size(site1(:,1),1);
for i=1:a
    L1=Lim_ALL{i}(:,1)'; L2=Lim_ALL{i}(:,2)';  
    k=1;site_new{i}=[]; Zt_new{i}=[];
    if size(site1,2)==2
        for j=1:n
           if ((site1(j,1)>=L1(1))  && (site1(j,1)<L2(1)) && (site1(j,2)>=L1(2))&& (site1(j,2)<L2(2)))  
            site_new{i}(k,:) = site1(j,:);
            Zt_new{i}(k,:)=Zt(j,:);
            k=k+1;
           end
        end
    end
    if size(site1,2)==1
        for j=1:n
           if ((site1(j,1)>=L1(1))  && (site1(j,1)<L2(1)))  
            site_new{i}(k,:) = site1(j,:);
            Zt_new{i}(k,:)=Zt(j,:);
            k=k+1;
           end
        end
    end
    if size(site1,2)==3
       for j=1:n
           if ((site1(j,1)>=L1(1))  && (site1(j,1)<L2(1)) && (site1(j,2)>=L1(2))&& (site1(j,2)<L2(2)) && (site1(j,3)>=L1(3))&& (site1(j,3)<L2(3))) 
            site_new{i}(k,:) = site1(j,:);
            Zt_new{i}(k,:)=Zt(j,:);
            k=k+1;
           end
        end
    end
    if size(site1,2)==4
       for j=1:n
           if ((site1(j,1)>=L1(1))  && (site1(j,1)<L2(1)) && (site1(j,2)>=L1(2))&& (site1(j,2)<L2(2)) && (site1(j,3)>=L1(3))&& (site1(j,3)<L2(3)) && (site1(j,4)>=L1(4))&& (site1(j,4)<L2(4))) 
            site_new{i}(k,:) = site1(j,:);
            Zt_new{i}(k,:)=Zt(j,:);
            k=k+1;
           end
        end
    end
    if size(site1,2)==5
       for j=1:n
           if ((site1(j,1)>=L1(1))  && (site1(j,1)<L2(1)) && (site1(j,2)>=L1(2))&& (site1(j,2)<L2(2)) && (site1(j,3)>=L1(3))&& (site1(j,3)<L2(3)) && (site1(j,4)>=L1(4))&& (site1(j,4)<L2(4)) && (site1(j,5)>=L1(5))&& (site1(j,5)<L2(5))) 
            site_new{i}(k,:) = site1(j,:);
            Zt_new{i}(k,:)=Zt(j,:);
            k=k+1;
           end
        end
    end
    if size(site1,2)==6
       for j=1:n
           if ((site1(j,1)>=L1(1))  && (site1(j,1)<L2(1)) && (site1(j,2)>=L1(2))&& (site1(j,2)<L2(2)) && (site1(j,3)>=L1(3))&& (site1(j,3)<L2(3)) && (site1(j,4)>=L1(4))&& (site1(j,4)<L2(4)) && (site1(j,5)>=L1(5))&& (site1(j,5)<L2(5)) ...
                   && (site1(j,6)>=L1(6))&& (site1(j,6)<L2(6))) 
            site_new{i}(k,:) = site1(j,:);
            Zt_new{i}(k,:)=Zt(j,:);
            k=k+1;
           end
        end
    end
    if size(site1,2)==7
       for j=1:n
           if ((site1(j,1)>=L1(1))  && (site1(j,1)<L2(1)) && (site1(j,2)>=L1(2))&& (site1(j,2)<L2(2)) && (site1(j,3)>=L1(3))&& (site1(j,3)<L2(3)) && (site1(j,4)>=L1(4))&& (site1(j,4)<L2(4)) && (site1(j,5)>=L1(5))&& (site1(j,5)<L2(5)) ...
                   && (site1(j,6)>=L1(6))&& (site1(j,6)<L2(6)) && (site1(j,7)>=L1(7))&& (site1(j,7)<L2(7))) 
            site_new{i}(k,:) = site1(j,:);
            Zt_new{i}(k,:)=Zt(j,:);
            k=k+1;
           end
        end
    end
    if size(site1,2)==8
       for j=1:n
           if ((site1(j,1)>=L1(1))  && (site1(j,1)<L2(1)) && (site1(j,2)>=L1(2))&& (site1(j,2)<L2(2)) && (site1(j,3)>=L1(3))&& (site1(j,3)<L2(3)) && (site1(j,4)>=L1(4))&& (site1(j,4)<L2(4)) && (site1(j,5)>=L1(5))&& (site1(j,5)<L2(5)) ...
                   && (site1(j,6)>=L1(6))&& (site1(j,6)<L2(6)) && (site1(j,7)>=L1(7))&& (site1(j,7)<L2(7))  && (site1(j,8)>=L1(8))&& (site1(j,8)<L2(8))) 
            site_new{i}(k,:) = site1(j,:);
            Zt_new{i}(k,:)=Zt(j,:);
            k=k+1;
           end
        end
    end
    nn{i}=k-1;
end
end
