
function log_prob = logBeta_pprob(a,b,x) 
log_prob=(a-1)*log(x)+(b-1)*log(1-x);
end