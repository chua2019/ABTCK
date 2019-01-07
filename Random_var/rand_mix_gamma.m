function rand_mixg=rand_mix_gamma(a1,b1,a2,b2,probm)
% gramacy and lee a1=1, b1=20; a2=10; b2=10; probm =0.5;
if rand < probm
    rand_mixg= rngGamma(a1,b1); 
else 
    rand_mixg= rngGamma(a2,b2);
end
end