function draw=Wishart_r(SIGMA,df)
d=size(SIGMA,1);
L=chol(SIGMA);
n=df;
X=zeros(d,n);
for i=1:n
    Z=norm_r(0,1,d);
    X(:,i)=L'*Z;
end
W=X*X';
draw=W;
end

