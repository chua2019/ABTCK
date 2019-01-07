            
function [t_val]=t_star(sigma2,AA,rho,x,D,S,tau2)
 % Revisit this one      
    for t=1:S
    if t==1
%         for j=1:(S-1)
%          rho_com=rho_com.*rho{j};
%         end
        rho_com=1;
        for j=1:(S-1)
         rho_com=rho_com.*rho{j+1};
        end
        
      t_st{1}=rho_com*sigma2{1}.*Gauss_Cov2(1, AA{1}, x,D{1},tau2); 
    end
    if t>1
        rho_com=1;
        for j=t:(S-1)
         rho_com=rho_com.*rho{j+1};
        end
        
        
        rho_com1{1}=1;
         for jj=1:(S-1)
         rho_com1{1}=rho_com1{1}.*rho{jj+1};
         end          
        t_stPrim{t}{1}=(rho_com1{1}).*sigma2{1}.*Gauss_Cov2(1, AA{1}, x,D{t},tau2);    
        for j=2:t
        
        rho_com1{j}=1;
         for jj=j:(S-1)
         rho_com1{j}=rho_com1{j}.*rho{jj+1};
         end    
         
        t_stPrim{t}{j}=rho{j-1+1}.*t_stPrim{t}{j-1}+(rho_com1{j}).*sigma2{j}.*Gauss_Cov2(1, AA{j}, x,D{t},tau2);
        end
        t_st{t}=rho{t-1+1}.*t_stPrim{t}{t-1}+(rho_com).*sigma2{t}.*Gauss_Cov2(1, AA{t}, x,D{t},tau2);
    end
    end
    t_val=[transpose(t_st{1})];
    
    for j=2:S
        t_val=[t_val; transpose(t_st{j})];
    end
    
end