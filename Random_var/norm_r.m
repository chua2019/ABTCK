function draws=norm_r(mu,sd, n)
if mod(n,2)==0
    for i=1:(n/2)
        U1=rand;
        U2=rand;
        T=2*pi*U2;
        z1=((-2*log(U1))^0.5)*cos(T);
        z2=((-2*log(U1))^0.5)*sin(T);
        N1(i)=mu+z1*sd;
        N2(i)=mu+z2*sd;
    end
    draws=[N1,N2]';
else
    for i=1:((n+1)/2)
        U1=rand;
        U2=rand;
        T=2*pi*U2;
        z1=((-2*log(U1))^0.5)*cos(T);
        z2=((-2*log(U1))^0.5)*sin(T);
        N1(i)=mu+z1*sd;
        N2(i)=mu+z2*sd;
    end
    DroWs=[N1,N2]';
    draws=DroWs(1:n);
end

end
