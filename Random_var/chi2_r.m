
function drow_chi2=chi2_r(df,n)
% drow from chi-squar
drow_chi2=zeros(n,1);
for i=1:n
Z=norm_r(0,1, df);
drow_chi2(i)=sum(Z^2);
end

end



