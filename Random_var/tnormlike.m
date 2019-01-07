function result = tnormlike(x,mu,sigma,left,right)
% PURPOSE: random draws from a normal truncated to (left,right) interval
% ------------------------------------------------------
% USAGE: y = rand_nort(mu,sigma2,left,right)
% where:   mu = mean (nobs x 1)
%      sigma2 = variance (nobs x 1)
%        left = left truncation points (nobs x 1)
%       right = right truncation points (nobs x 1)
% ------------------------------------------------------
% RETURNS: y = (nobs x 1) vector
% ------------------------------------------------------
% NOTES: use y = rand_nort(mu,sigma2,left,mu+5*sigma2)
%        to produce a left-truncated draw
%        use y = rand_nort(mu,sigma2,mu-5*sigma2,right)
%        to produce a right-truncated draw
% ------------------------------------------------------
% 
% adopted from 
% James P. LeSage, Dept of Economics
% and
% Ordinal Data Modeling by Valen Johnson and James Albert
% Springer-Verlag, New York, 1999.

 
  std =  sigma;
% Calculate bounds on probabilities
  lowerProb = normcdf(left,mu,sigma);
  upperProb = normcdf(right,mu,sigma);
% Draw uniform from within (lowerProb,upperProb)
  result=normlike([mu,sigma],x)+log(upperProb-lowerProb); 
 