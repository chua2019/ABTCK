function [quartile1]=quantile_f(All_Var_kr,a)
n=size(All_Var_kr,1);
All_Var_krS=sortrows(All_Var_kr);
val=int16(n*a);
quartile1=All_Var_krS(val);
end




