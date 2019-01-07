function [y]=funn_ALL(x,funname)

switch funname
    case 'first'
y=x(:,1).*exp(-x(:,1).^2-x(:,2).^2);
    case 'second'
y=exp(sin((0.9.*(x(:,1)+0.48)).^(10)))+x(:,2).*0.99+0.1;
    case 'third'
y=exp(sin((0.9.*(x(:,1)+0.48)).^(10)))+x(:,2).*x(:,3)+x(:,4);
    case 'funn2out2'
        y=funn2out2(x(:,1),x(:,2));
    case 'ko-2'
        n=size(x,1);
        for i=1:n
            [Y(i,:)] = ko_solver(x(i,:),10);
        end
        y=Y(:,3);
    case 'Mult_ko-2'
        n=size(x,1);
        for i=1:n
            [Y(i,:)] = ko_solver(x(i,:),10);
        end
        y=Y;
        
    case 'Mult_time_ko-2'
n=size(x,1); y=[]; 
        for t=1:T
            for i=1:n
                [Y(i,t,:)] = ko_solver(x(i,:),t);
            end
           %y(((t-1)*n+1):(n*t),:)=Y(:,t,1:3);
        end
%         i=1:3 aaa(i)=Y(:,:,3)';
%          y=reshape(aaa,T*n,1);
        for ii=1:n
         y(((ii-1)*T+1):(ii*T),:)=Y(ii,:,1:3); % depending on the way we present y we determine the cov_matrix. 
        end
    case 'Mult_time_ko-1'
    n=size(x,1); y=[]; 
        for t=1:T
            for i=1:n
                [Y(i,t,:)] = ko_solver(x(i,:),t);
            end
        end
        for ii=1:n
         y(((ii-1)*T+1):(ii*T),:)=Y(ii,:,1:3); % depending on the way we present y we determine the cov_matrix. 
        end   
    case 'calib1'
        y=funnWilCal2d_Dis(x(:,1),x(:,2));
    case 'calib2'
        y=funnWilCal2d(x(:,1),x(:,2));
end
end