function [y]=funn4WilCald_Dis(x1,x2)

discon1=0.30;
discon2=0.40;

if (size(x2,1) ==1)
for i=1:size(x1,1) 
    for j=1:size(x1,2)         
        if x1(i,j)<discon1
        y(i,j)=(1-x1(i,j)).*cos(pi.*x2);
        end              
            if (x1(i,j)>=discon1 && x2<discon2)
                y(i,j)=(1-x1(i,j)).*cos(pi.*x2)-0.5;
            end
            if (x1(i,j)>=discon1 && x2>=discon2)
                y(i,j)=(1-x1(i,j)).*cos(pi.*x2)+0.5;
            end
    end
end
else
for i=1:size(x2,1) 
    for j=1:size(x2,2) 
        if x1(i,j)<discon1
            y(i,j)=(1-x1(i,j)).*cos(pi.*x2(i,j));
        end
        if (x1(i,j)>=discon1 && x2(i,j)<discon2)   
            y(i,j)=(1-x1(i,j)).*cos(pi.*x2(i,j))-0.5;
        end
        if (x1(i,j)>=discon1 && x2(i,j)>=discon2) 
            y(i,j)=(1-x1(i,j)).*cos(pi.*x2(i,j))+0.5;
        end
    end
end

end