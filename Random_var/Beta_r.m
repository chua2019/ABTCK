%
% RANDOM BETA NUMBER GENERATOR 
%
% IT GENERATES RANDOM NUMBERS FROM A BETA(a,b)DISTRIBUTION 
% =========================================================================
%
% CALL :
% Z=Beta_r(alpha,beta)
% Y = rngGamma( A, B ) , A > 0
% n = sample size
% =========================================================================
%
% CODED BY :
%
% B. Konomi
% @PNNL
% 28/09/2012
% =========================================================================

function z = Beta_r(a,b,n) 
z=zeros(n,1);
for i=1:n
x = rngGamma(a,1); y = rngGamma(b,1);
z(i)=x/(x+y);
end
end

% =========================================================================
