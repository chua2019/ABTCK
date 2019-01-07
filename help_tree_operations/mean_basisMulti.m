

function [Gt]=mean_basisMulti(site,Zt,S)
  nn{1}=size(site{1},1);
    Gt{1}=[ones(nn{1},1)];
for t=2:S
    %Constant mean only 
    nn{t}=size(site{t},1);
     X1 = [ones(nn{t},1)]; 
    ISM=ismember(site{t-1},site{t});
     X2=[Zt{t-1}(ISM(:,1))];
     Gt{t}=[X1  X2];

end 
end