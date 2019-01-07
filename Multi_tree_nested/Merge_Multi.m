function [site,Zt, n, prCurr1,AA, a, lim_ALL,NODES,a_NODES]=...
    Merge_Multi(prCurr,site,Zt,n,AA,S,lim_ALL,a,a_NODES,NODES, aalpha, bbeta,tau1,tau2,funname)
% Separate  overall region into smaller subregions - start here with discontinuity on eta(x,t); and then generalize to d(x) later;
% Delete Separate the Dim_x-Dimentiona area into smaller areas. 
% Possible merging subsets % TREE will record the tree in every step... The TREE will consist of
%  a cells with it own depth. e.g. [0 1] is the left and then right node.
% We need to have that many trees as final nodes. This means we can follow
% the same as for the other parameters.
% If we want to record the depth of the tree we can measure the maximum of
% the length of each components of the tree.

% a : Nuber of external Nodes.
% a_NODES: Number of Nodes.
% NODES: all nodes of a tree - if we know this we can find the tree.
% lim_ALL: for every node we have a rectangular region in the input space.

% NODES{i}.number; % a node numbe k with have (2*k) and (2*k+1) node childs
% Nodes{i}.inter = 0 or 1 : if interior 1 if exterior 0.
% NODES{i}.parent_child % here whether or not is a left or right child (0 if left 1 if right).
% NODES{i}.direction % here the direction we cut it ([0] if it is child).
% NODES{i}.rul  % here the rule we use --> a value between 0 and 1.

Dim_x = size(site{1}{1},2); 
Dim_y = size(Zt{1}{1},2);
a_NODES_new =a_NODES-2;

%NODES_new=cell(a_NODES_new);
% find all the prounable nodes.

for jj=1:a_NODES
    NO_number{jj} = NODES{jj}.number;
end
bb=0;
for  ii=1:a_NODES
    if (NODES{ii}.inter == 1)
        num_1=NODES{ii}.number;
        CH_L=find([NO_number{:}]==2*num_1); % Left child
        CH_R=find([NO_number{:}]==(2*num_1+1)); % Right child
        if (NODES{CH_L}.inter==0 && NODES{CH_R}.inter==0); % if both of the childerens are exterior.
            bb=bb+1; PP(bb)=ii;  depth_each(bb)=floor(log(NODES{ii}.number)/log(2));% or better asta(bb)=ii;
        end
    end
end
    kk1=randi(bb); k = PP(kk1); K = NODES{k}.number;
    depth = depth_each(kk1);
% Proune the left and right child(of parent NODE k).
% find the left and right child to delite
        k_LP=find([NO_number{:}]==2*K); % Left child
        k_RP=find([NO_number{:}]==(2*K+1)); % Right child
a_new=a-1; 
% delete the two NODES and name new NODES
if K==1    
    [nn1{1}]=n{1};
    [AA_new1{1}]=AA{k_LP}; 
    [prProp{1}]= prCurr{1};  
    [lim_ALL_new{1}]=lim_ALL{1};      
    [site_new1{1}]=site{1};       
    [Zt_new1{1}]=Zt{1};        
    
    [NODES_new{1}]=NODES{1};
else     
    [AA_new1{1:a_NODES_new}]=AA{[1:(k_LP-1) (k_LP+2):a_NODES]};   
    [prProp{1:a_NODES_new}]= prCurr{[1:(k_LP-1) (k_LP+2):a_NODES]};
    [lim_ALL_new{1:a_NODES_new}]=lim_ALL{[1:(k_LP-1) (k_LP+2):a_NODES]};     
%  for t=1:S   
%     [site_new1{1:a_NODES_new}{t}]=site1{[1:(k_LP-1) (k_LP+2):a_NODES]}{t};
%  end
    [site_new1{1:a_NODES_new}]=site{[1:(k_LP-1) (k_LP+2):a_NODES]};
    [Zt_new1{1:a_NODES_new}]=Zt{[1:(k_LP-1) (k_LP+2):a_NODES]};        
    [nn1{1:a_NODES_new}] = n{[1:(k_LP-1) (k_LP+2):a_NODES]};
    [NODES_new{1:a_NODES_new}]=NODES{[1:(k_LP-1) (k_LP+2):a_NODES]};           
end

% In the new setting find the parent which we merge the childrens and
% change the parameter to be one of the child parameters.
for jj=1:a_NODES_new
NO_number_new{jj}=NODES_new{jj}.number;
end
k_new=find([NO_number_new{:}]==K); 
NODES_new{k_new}.inter = 0; % make this exterior node.

% GP parameters which change
% randomly chose the child which with give its parameters.

randi(2);
if (randi(2)==1)
    AA_new1{k_new}= AA{k_LP}; 
    %Proposal_AA = LogWishart_like(AA_new1{k_new},AA{k_RP},Dim_y,Dim_y); %Wishart_like(AA_new{k_new},AA{k_RP},Dim_y,Dim_y);
else 
     AA_new1{k_new}= AA{k_RP};
     %Proposal_AA = LogWishart_like(AA_new1{k_new},AA{k_LP},Dim_y,Dim_y); %Wishart_like(AA_new{k_new},AA{k_LP},Dim_y,Dim_y);
end
% X=mean_basis_eta(site_new2{k_new}, funname2);
% Sigma_eta_hat = sigma_eta_hat(Zt_new2{k_new},AA_new1{k_new}, tau2_Y{k}, site_new2{k_new},X,funname);  

for t=S
nn1{k_new}{t}=size(site_new1{k_new}{t},1);
end
% Posterior  value of at the proposed value.

X_s{k_new}=mean_basisMulti(site_new1{k_new},Zt_new1{k_new},S);
for t=1:S
[mu_Zt11, SIGMA11, beta_X11, prProp11{t}]=Multi_mu_s3(Zt_new1{k_new}{t},AA_new1{k_new}{t}, site_new1{k_new}{t},X_s{k_new}{t},tau1,tau2);
 prProp{k_new}{t}=prProp11{t};  
end
% calculate the depth of the existing tree
hh2=0;
for i=1:a_NODES
    if (NODES{i}.inter==0)
        hh2=hh2+1; depth_each(hh2)=floor(NODES{i}.number/2);
    end
end
%depth=max(depth_each); 
% Calculate prior.
% so we can use the exact opposite from the one used in split.
% G is the set of gronable nodes
G=a_new; %Grownable nodes of T' tree.
% 2*size(PP) is the set of prunable nodes
P=2*size(PP,2); % Prounable node of T tree.

% some logical values for aalpha=0.6; bbeta=3; 
Prior_prop=log(1-aalpha*(1+depth)^(-bbeta));
Prior_curr=log(aalpha)+(-bbeta)*log((1+depth))+2*log(1-aalpha*(2+depth)^(-bbeta));
  
%% MH ratio
 prod_prCurr1=0;
 prod_prCurr2=0;
 prod_prProp=0;
for t=1:S
 prod_prCurr1=prod_prCurr1+prCurr{k_LP}{t};
 prod_prCurr2=prod_prCurr2+prCurr{k_RP}{t};
 prod_prProp=prod_prProp+prProp{k_new}{t};
end

  %2*(Prior_prop-Prior_curr) % Proposal_AA
  MH_ratio = (prod_prProp+0-prod_prCurr1-prod_prCurr2)+1*(Prior_prop-Prior_curr)+(log(P)-log(G));
u = log(rand);%0;%
%MH_ratio-u
if (u<MH_ratio)
% accept
    site=site_new1; 
    Zt=Zt_new1;
    n=nn1; 
    prCurr1 = prProp; 
    a=a_new;
    AA=AA_new1;    
    lim_ALL=lim_ALL_new;
    NODES=NODES_new;
    a_NODES = a_NODES_new;
else 
     prCurr1 = prCurr; 
end    
end

