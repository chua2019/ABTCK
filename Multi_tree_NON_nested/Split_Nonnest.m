function [site,Zt, site_real, Zt_real, n, prCurr1, AA,lim_ALL, a,Nodes,a_Nodes]= ...
    Split_Nonnest(prCurr,site,Zt,site_real,Zt_real,n,AA,S,low_l,lim_ALL,a,a_NODES,NODES, pprop, aalpha, bbeta,tau1,tau2,funname)
% Separate  overall region into smaller subregions - start here with discontinuity on eta(x,t); and then generalize to d(x) later;
% L_R: will give the a value to the left or the right child. So we
% will be able to use it in the merge step.
% Split k cell - ellement. n=size(site1(:,1),1);
% a : Nuber of external Nodes.
% a_NODES: Number of Nodes.
% NODES: all nodes of a tree - if we know this we can find the tree.
% lim_ALL: for every node we have a rectangular region in the input space.
% Zt1 is the observation output
% Zt2 is the computer code run output 
% site1 is the observation input (including theta -- D_1(theta))
% site2 is the computer code run output D_2.

% NODES{i}.number; % a node numbe k with have (2*k) and (2*k+1) node childs
% Nodes{i}.inter = 0 or 1 : if interior 1 if exterior 0.
% NODES{i}.parent_child % here whether or not is a left or right child (0 if left 1 if right).
% NODES{i}.direction % here the direction we cut it ([0] if it is child).
% NODES{i}.rul  % here the rule we use --> a value between 0 and 1.

Dim_x = size(site{1}{1},2); Dim_y = size(Zt{1}{1},2);
a_NODES_new =a_NODES+2; hh2=0;
for i=1:a_NODES
    if (NODES{i}.inter==0)
        hh2=hh2+1; fin_nodes(hh2)=i;
        depth_each(hh2)=floor(log(NODES{i}.number)/log(2));
    end
end
a_new=a+1;
kk1=randi(a);
k=fin_nodes(kk1);
% Change the tree and its parameters only for the Node{k} and Node{k+1}
depth=depth_each(kk1);%max(depth_each);
% GP parameters
AA_new1=cell(a_NODES_new,1);
% tau2_Y_new=cell(a_NODES_new,1);
% tau2_Z_new=cell(a_NODES_new,1);
prProp = cell(a_NODES_new,1); 
     % Tree 
     [NODES_new{1:(a_NODES)}]=NODES{1:a_NODES};
     % GP parameters
     %for t=1:S
     for i=1:a_NODES
     [AA_new1{i}]=AA{i};         
     %[AA_new1{1:(a_NODES),1}]=AA{1:a_NODES};          
     %[SIGMA_new1{1:(a_NODES),1}]=SIGMA{1:a_NODES}; 
     end
     %end 
     %Posterior
     [prProp{1:(a_NODES)}]=prCurr{1:a_NODES};
     % Limits 
     [lim_ALL_new{1:(a_NODES)}]=lim_ALL{1:a_NODES};
     % New site and outputs for the new partition
     [site_new1{1:(a_NODES)}]=site{1:a_NODES};     
     %site_new1
    [site_real_new{1:(a_NODES)}]=site_real{1:a_NODES};
    
     [Zt_new1{1:(a_NODES)}]=Zt{1:a_NODES};
     [Zt_real_new{1:(a_NODES)}]=Zt_real{1:a_NODES};
     [nn1{1:(a_NODES)}]=n{1:a_NODES};  
% Split in the Dir-th direction.
Dir=randi(Dim_x);
ran=0.2+0.6*rand;%
        L_new{1}=lim_ALL{k};        
        L_new{1}(Dir,1)=lim_ALL{k}(Dir,1);
        L_new{1}(Dir,2)=lim_ALL{k}(Dir,1)+((ran)*(lim_ALL{k}(Dir,2)-lim_ALL{k}(Dir,1)));      
        
        L_new{2}=lim_ALL{k};            
        L_new{2}(Dir,1)=lim_ALL{k}(Dir,1)+((ran)*(lim_ALL{k}(Dir,2)-lim_ALL{k}(Dir,1)));
        L_new{2}(Dir,2)=lim_ALL{k}(Dir,1)+((lim_ALL{k}(Dir,2)-lim_ALL{k}(Dir,1)));          
        lim_ALL_new{a_NODES+1}=L_new{1}; lim_ALL_new{a_NODES+2}=L_new{2};
                Ssite_new1{1}=cell(1,S); Ssite_new1{2}=cell(1,S); 
        ZZt_new1{1}=cell(1,S); ZZt_new1{2}=cell(1,S); 
       
        for t=1:S
        Nnn1{1}{t}=0; Nnn1{2}{t}=0;
        [Ssite_new1{1}{t}, ZZt_new1{1}{t}, Nnn1{1}{t}]= RE_siteLim(site{k}{t},Zt{k}{t},L_new{1});
        [Ssite_new1{2}{t}, ZZt_new1{2}{t}, Nnn1{2}{t}]= RE_siteLim(site{k}{t},Zt{k}{t},L_new{2});
        
        [site_real_new{(a_NODES+1)}{t}, Zt_real_new{(a_NODES+1)}{t}, nn1{1}{t}]= RE_siteLim(site_real{k}{t},Zt_real{k}{t},L_new{1});
        [site_real_new{(a_NODES+2)}{t}, Zt_real_new{(a_NODES+2)}{t}, nn1{2}{t}]= RE_siteLim(site_real{k}{t},Zt_real{k}{t},L_new{2});
        end

        %Ssite_new2{1}=[]; Ssite_new2{2}=[]; ZZt_new2{1}=[]; ZZt_new2{2}=[];Nnn2{1}=0; Nnn2{2}=0;
%         [Ssite_new2{1}, ZZt_new2{1}, Nnn2{1}]= RE_siteLim(site2{k},Zt2{k},L_new{1});
%         [Ssite_new2{2}, ZZt_new2{2}, Nnn2{2}]= RE_siteLim(site2{k},Zt2{k},L_new{2});
        

% Change the parent node to be interior.
 NODES_new{k}.inter =1; % interior         
% Left child 
 NODES_new{a_NODES+1}.number =2*NODES{k}.number; % as the left chiled
 NODES_new{a_NODES+1}.inter =0; % exterior 
 NODES_new{a_NODES+1}.parent_child=0; % Left child 
 NODES_new{a_NODES+1}.direction=Dir;
 NODES_new{a_NODES+1}.rul=lim_ALL_new{(a_NODES+1)}(Dir,2); % uper bound
% Right child 
 NODES_new{a_NODES+2}.number =2*NODES{k}.number+1; % as the right child
 NODES_new{a_NODES+2}.inter =0; % exterior 
 NODES_new{a_NODES+2}.parent_child=1;% Right child 
 NODES_new{a_NODES+2}.direction=Dir;
 NODES_new{a_NODES+2}.rul=lim_ALL_new{(a_NODES+2)}(Dir,1); % lower bound

 
 Min_n1=min(nn1{1}{:});
 Min_n2=min(nn1{2}{:});
 Min_n=min(Min_n1,Min_n2);
if (Min_n>low_l)% (Nnn1{1}>low_limit1 && Nnn1{2}>low_limit1 && Nnn2{1}>low_limit2 && Nnn2{2}>low_limit2)
    % location and observation of the new node           
    site_new1{1,(a_NODES+1)}=Ssite_new1{1}; site_new1{1,(a_NODES+2)}=Ssite_new1{2};    
    Zt_new1{1,(a_NODES+1)}=ZZt_new1{1}; Zt_new1{1,(a_NODES+2)}=ZZt_new1{2};  
    %for t=1:S
    nn1{1,(a_NODES+1)}=Nnn1{1}; nn1{1,(a_NODES+2)}=Nnn1{2};   
    %end  
    % parameter in the new nodes  
    a1=pprop.a1; b1=pprop.b1; a2=pprop.a2; b2=pprop.a2; probm =pprop.probm;% a1=2; b1=3; a2=8; b2=2; probm =0.5;% % a1=2; b1=1; a2=4; b2=2; probm =0.5;% a1=2; b1=2; a2=20; b2=2; probm =0.5; %
    for t=1:S
    Proposal_AA{t}=0;
    end
    if randi(2)==1
        AA_new1{(a_NODES+1)}=AA{k};
    for h1=1:Dim_x
    %h1=Dir;
    for t=1:S
    AA_new1{(a_NODES+2)}{t}(h1,h1)=rand_mix_gamma(a1,b1,a2,b2,probm);
    Proposal_AA{t}=Proposal_AA{t}+log(probm*gamma_lik(AA_new1{(a_NODES+1)}{t}(h1,h1),a1,b1)+(1-probm)*gamma_lik(AA_new1{(a_NODES+1)}{t}(h1,h1),a2,b2));
    end
    end    
    % this may be changed
    else
        AA_new1{(a_NODES+2)}=AA{k};
    for h1=1:Dim_x
    %h1=Dir;
    for t=1:S
    AA_new1{(a_NODES+1)}{t}(h1,h1)=rand_mix_gamma(a1,b1,a2,b2,probm);
    Proposal_AA{t}=Proposal_AA{t}+log(probm*gamma_lik(AA_new1{(a_NODES+2)}{t}(h1,h1),a1,b1)+(1-probm)*gamma_lik(AA_new1{(a_NODES+2)}{t}(h1,h1),a2,b2));
    end
    end
    AA_new1{(a_NODES+1)}=AA{k}; 
    end
    
    % Compute q(T*,T) and q(T,T*); G is the set of gronable nodes;
    G=a;
    % P is the set of prunable nodes (tree before spliting)
    bb=0; PP=[];
    for  ii=1:a_NODES
        if (NODES{ii}.inter == 1)
            num_1=NODES{ii}.number;
            for jj=1:a_NODES
                NO_number{jj} = NODES{jj}.number;
            end
            CH_L=find([NO_number{:}]==2*num_1); % Left child
            CH_R=find([NO_number{:}]==(2*num_1+1)); % Right child
            if (NODES{CH_L}.inter==0 && NODES{CH_R}.inter==0) % if both of the childerens are exterior.
                bb=bb+1; PP(bb)=ii; % or better asta(bb)=ii;
            end
        end
    end
    P=2*size(PP,1);
    % posterior for the new node n_NODE+1
    

X_s{(a_NODES+1)}=mean_basisMulti(site_new1{(a_NODES+1)},Zt_new1{(a_NODES+1)},S);
X_s{(a_NODES+2)}=mean_basisMulti(site_new1{(a_NODES+2)},Zt_new1{(a_NODES+2)},S);
 

for t=1:S
    [mu_Zt11, SIGMA11, beta_X11, prProp11{t}]=Multi_mu_s3(Zt_new1{(a_NODES+1)}{t},AA_new1{(a_NODES+1)}{t}, site_new1{(a_NODES+1)}{t},X_s{(a_NODES+1)}{t},tau1,tau2);
    [mu_Zt12, SIGMA21, beta_X12, prProp12{t}]=Multi_mu_s3(Zt_new1{(a_NODES+2)}{t},AA_new1{(a_NODES+2)}{t}, site_new1{(a_NODES+2)}{t},X_s{(a_NODES+2)}{t},tau1,tau2);    
end
 prProp{(a_NODES+1)}=prProp11;  
 prProp{(a_NODES+2)}=prProp12;  
 prod_prProp1=0;
 prod_prProp2=0;
 prod_prCurr=0;
 prod_Proposal_AA=0;
for t=1:S
 prod_prProp1=prod_prProp1+prProp{(a_NODES+1)}{t};
 prod_prProp2=prod_prProp2+prProp{(a_NODES+2)}{t};
  prod_prCurr=prod_prCurr+prCurr{k}{t};
  prod_Proposal_AA=prod_Proposal_AA+Proposal_AA{t};
end
%Proposal_AA
    % Calculate prior. 
    Prior_curr=log(1-aalpha*(1+depth)^(-bbeta));
    Prior_prop=log(aalpha)+(-bbeta)*log((1+depth))+5*log(1-aalpha*(2+depth)^(-bbeta));
    % MH
    MH_ratio =(prod_prProp1+prod_prProp2-prod_prCurr-(prod_Proposal_AA))+5*(Prior_prop-Prior_curr)+(log(G)-log(P)); % ftia3e edw to RJ of split.
    u = log(rand);%
    if (u<MH_ratio)
    % accept
        % Updated tree
        Nodes = NODES_new; lim_ALL=lim_ALL_new; a_Nodes=a_NODES_new; 
        % updated locations,observations and number. 
        site=site_new1; Zt=Zt_new1; n=nn1;
        site_real=site_real_new; Zt_real=Zt_real_new;  
        % updated posterior
        prCurr1 = prProp; 
        % updated GP parameters 
        AA=AA_new1; 
        a=a_new; % external nodes
    else 
        prCurr1 = prCurr;  Nodes=NODES; a_Nodes= a_NODES;
    end
else 
    prCurr1 = prCurr; Nodes=NODES; a_Nodes= a_NODES;
end        
end