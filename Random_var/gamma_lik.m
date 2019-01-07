function like=gamma_lik(x,a,b)
% likelihood
like =(1/gamma(a))*x^(a-1)*exp(-b*x);
end