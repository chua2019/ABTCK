function like=Wishart_like(w,SIGMA,n,d)

% this is the proportional distribution - one constant is missing.
like =(det(SIGMA))^(-n/2)*(det(w))^((n-d-1)/2)*exp(-trace((SIGMA^(-1)*w))/2);
end