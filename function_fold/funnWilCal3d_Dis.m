function [y]=funnWilCal3d_Dis(x1,x2)

discon=0.45;
if (size(x2,1) ==1)
for i=1:size(x2,1) 
    for j=1:size(x2,2) 
if x2(i,j)<discon
    y(i,j)=(1-x1(i,j)).*cos(pi.*x2);
end
if x2(i,j)>=discon
    y(i,j)=(1-x1(i,j)).*cos(pi.*x2)+0.5;
end
    end  
end
else
for i=1:size(x2,1) 
    for j=1:size(x2,2) 
if x2(i,j)<discon
    y(i,j)=(1-x1(i,j)).*cos(pi.*x2(i,j));
end
if x2(i,j)>=discon
    y(i,j)=(1-x1(i,j)).*cos(pi.*x2(i,j))+0.5;
end
    end
end
end