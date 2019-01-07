function log_like=LogWishart_like(w,SIGMA,n,d)
% log likelihood
P_gama=1;
for i=1:d
    P_gama =P_gama*gamma((n+1-i)/2);
end
log_like =(-n/2)*log(det(SIGMA))+((n-d-1)/2)*(det(w))+(-trace((SIGMA^(-1)*w))/2)-n*d/2*log(2)-d*(d-1)/4*log(pi)-log(P_gama);
end