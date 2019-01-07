function [y]=funnLoic(x1,x2,x3,funname)
a=7;
b=2;
dis=0;
c=2;
switch funname
    case 'level1'
        y=(sin(pi.*x1)+(sin(pi.*x2).^2));
    case 'level2'

            y=sin(pi.*x1)+a.*(cos(pi.*x2).^2)+(b.*(pi.*x3).^2);

    case 'level3'
        for i=1:size(x3,1)
       if x3(i)>dis 
        y(i,:)=2*sin(pi.*x1(i))+(a.*(sin(pi.*x2(i)).^2))+(b.*(pi.*x3(i)).^2).*(sin(pi.*x1(i)));
       else
        y(i,:)=((sin(pi.*x1(i))+a.*(sin(pi.*x2(i)).^2)))+exp(sin((0.9.*(x1(i)+0.48)).^(10)));%+b.*(pi.*x3).^4
       end
        end
end
