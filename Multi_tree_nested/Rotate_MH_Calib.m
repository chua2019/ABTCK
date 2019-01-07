function [site1,site2,  Zt1, Zt2, nn, lim_ALL,AA_eta,SIGMA_eta,AA_delta, SIGMA_delta,tau2_Z, tau2_Y, prCurr,NODES]= Rotate_MH_Calib(prCurr,site1,site2,Zt1,Zt2,nn,AA_eta,AA_delta,SIGMA_eta,SIGMA_delta,p, tau2_Z, tau2_Y,lim_ALL,a_NODES,NODES, aalpha,bbeta)
% Add Separate the 2D area into smaller areas in the y direc. 
% BB = the subregion which will be splited.
% Find a parent NODE and change the spliting rull for his childrens.
% we have also to change all the limits of the decendant of this node, site and Zt (observation).
% Do not change the parameters.
% In THIS VERSION DO CHANGE ONLY THE PARENTS WITH EXTERNAL NODES.

% Zt1 is the observation output
% Zt2 is the computer code run output 
% site1 is the observation input (including theta -- D_1(theta))
% site2 is the computer code run output D_2.

% a : Nuber of external Nodes.
% a_NODES: Number of Nodes.
% NODES: all nodes of a tree - if we know this we can find the tree.
% lim_ALL: for every node we have a rectangular region in the input space.

% NODES{i}.number; % a node numbe k with have (2*k) and (2*k+1) node childs
% Nodes{i}.inter = 0 or 1 : if interior 1 if exterior 0.
% NODES{i}.parent_child % here whether or not is a left or right child (0 if left 1 if right).
% NODES{i}.direction % here the direction we cut it .
% NODES{i}.rul  % here the rule we use.

lim_ALL_new =lim_ALL;
prProp=prCurr;
NODES_new = NODES; NODES_NEW = NODES;

SIGMA_eta1 = SIGMA_eta;
AA_eta1  = AA_eta;
AA_delta1=AA_delta;
SIGMA_delta1=SIGMA_delta;

tau2_Y_new =tau2_Y;

% find all the interior nodes and their coresponding parents (exclude node 1) if they spli in one direction.
for jj=1:a_NODES
    NO_number{jj} = NODES{jj}.number;
end

bb=0;
for  ii=2:a_NODES % without the initial node
    if (NODES{ii}.inter == 1)
            num_1=NODES{ii}.number;
            PARENT=find([NO_number{:}]==floor(num_1/2)); % coresponding PARENT NODE
            Children=find([NO_number{:}]==num_1*2);% coresponding left Child NODE
            if (NODES{Children}.direction ==NODES{ii}.direction)
                bb=bb+1; 
                PP(bb)=ii; % the coresponding NODE
                PA_L(bb)=find([NO_number{:}]==floor(num_1/2)); % coresponding PARENT NODE
                ROTA(bb,:) = [PP(bb) PA_L(bb)];
            end
    end
end
if bb>0
kk1=randi(bb); 
k_CH = ROTA(kk1,1); % Child node  
k_P = ROTA(kk1,2); % Parent node
% PARENT AND CHILD NUMBER 
K_CH = NODES{k_CH}.number; % Child number  
K_P =NODES{k_P}.number; % Parent number

if (NODES{k_CH}.parent_child==0)%mod(K_CH,K_P)==0
    K_OCH = K_CH+1;%other_child = right child
else 
    K_OCH = K_CH-1;% other_child = left child 
end
k_OCH=find([NO_number{:}]==floor(K_OCH));

        k_LCH=find([NO_number{:}]==2*K_CH); % Left child of the child inintrior nodes 
        k_RCH=find([NO_number{:}]==(2*K_CH+1)); % Right child  of the child inintrior nodes         
%k_LCH, k_RCH and K_OCH contain the T_1, T_2 and T_3 as in Gramacy and Lee.     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change number, parent_child and  in the NEW
% change the rul
NODES_new{k_CH}.rul = NODES{k_LCH}.rul; 
NODES_new{k_OCH}.rul = NODES{k_LCH}.rul; 
% change the interior directionrot
NODES_new{k_CH}.direction = NODES{k_LCH}.direction; 
NODES_new{k_OCH}.direction = NODES{k_RCH}.direction;  
% change the interior (the oposite to the previouse NODES)
if (NODES{k_CH}.parent_child==0)% if left child 
    NODES_new{k_CH}.inter = NODES_new{k_LCH}.inter;%make node such that 
else 
    NODES_new{k_CH}.inter = NODES_new{k_RCH}.inter;% 
end
NODES_new{k_OCH}.inter = 1; % Always become interior.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NODES_new{k_LCH}.number  = 2*K_OCH; 
NODES_new{k_RCH}.number  = 2*K_OCH+1;

% change the interior direction, rull, number
NODES_new{k_LCH}.direction = NODES{k_CH}.direction; 
NODES_new{k_LCH}.rul = NODES{k_CH}.rul; 

NODES_new{k_RCH}.direction = NODES{k_CH}.direction;  
NODES_new{k_RCH}.rul = NODES{k_CH}.rul; 
% this depends from left or right roation (Do first the left.)
if (NODES{k_CH}.parent_child==0)% if left child ---> right rotation
% C
    prProp{k_CH} = prCurr{k_RCH}; SIGMA_eta1{k_CH} = SIGMA_eta{k_RCH}; AA_eta1{k_CH} = AA_eta{k_RCH};
    SIGMA_delta1{k_CH} = SIGMA_delta{k_RCH}; AA_delta1{k_CH} = AA_delta{k_RCH};
    tau2_Y_new{k_CH} =tau2_Y{k_RCH}; 
% A    
    prProp{k_RCH} = prCurr{k_OCH}; SIGMA_eta1{k_RCH} = SIGMA_eta{k_OCH}; AA_eta1{k_RCH} = AA_eta{k_OCH};
    SIGMA_delta1{k_RCH} = SIGMA_delta{k_OCH}; AA_delta1{k_RCH} = AA_delta{k_OCH};
     tau2_Y_new{k_RCH} =tau2_Y{k_OCH};
    NODES_new{k_RCH}.inter = NODES{k_OCH}.inter;
% B 
    prProp{k_LCH} = prCurr{k_RCH}; SIGMA_eta1{k_LCH} = SIGMA_eta{k_OCH}; AA_eta1{k_LCH} = AA_eta{k_RCH};
    SIGMA_delta1{k_LCH} = SIGMA_delta{k_OCH}; AA_delta1{k_LCH} = AA_delta{k_RCH};
    tau2_Y_new{k_LCH} =tau2_Y{k_RCH};
    NODES_new{k_LCH}.inter = NODES{k_RCH}.inter;
else% if right child --> left rotaion
% C
    prProp{k_CH} = prCurr{k_LCH}; SIGMA_eta1{k_CH} = SIGMA_eta{k_LCH}; AA_eta1{k_CH} = AA_eta{k_LCH};
    SIGMA_delta1{k_CH} = SIGMA_delta{k_LCH}; AA_delta1{k_CH} = AA_delta{k_LCH};
    tau2_Y_new{k_CH} =tau2_Y{k_LCH};    
% A    
    prProp{k_RCH} = prCurr{k_LCH}; SIGMA_eta1{k_RCH} = SIGMA_eta{k_LCH}; AA_eta1{k_RCH} = AA_eta{k_LCH};
    SIGMA_delta1{k_RCH} = SIGMA_delta{k_LCH}; AA_delta1{k_RCH} = AA_delta{k_LCH};
    tau2_Y_new{k_RCH} =tau2_Y{k_LCH};
    NODES_new{k_RCH}.inter = NODES{k_LCH}.inter;
% B 
    prProp{k_LCH} = prCurr{k_OCH}; SIGMA_eta1{k_LCH} = SIGMA_eta{k_OCH}; AA_eta1{k_LCH} = AA_eta{k_OCH};
    SIGMA_delta1{k_LCH} = SIGMA_delta{k_OCH}; AA_delta1{k_LCH} = AA_delta{k_OCH};
    tau2_Y_new{k_LCH} =tau2_Y{k_OCH};
    NODES_new{k_LCH}.inter = NODES{k_OCH}.inter;   
end

%k_LCH, k_RCH and K_OCH contain the T_1, T_2 and T_3 as in Gramacy and Lee.    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:a_NODES
        dd(i)=floor(log(NODES{i}.number)/log(2)); %depth
end
% find the set of nodes which are in the A and B block (pinaka)

if (NODES{k_CH}.parent_child==0) % right rotation
 A_k=k_LCH; B_k=k_RCH; C_k=k_OCH; % left child of the child 
else 
 A_k=k_OCH; B_k=k_LCH; C_k=k_RCH; % right child of the child 
end

% find the set of nodes which are in the A 
l2=0;
for i=1:a_NODES
    if dd(i)>dd(A_k)
    if (2^(dd(i)-dd(A_k))*NO_number{A_k} <=NO_number{i} && 2^(dd(i)-dd(A_k))*(NO_number{A_k}+1) >NO_number{i}) % between [2K,
        l2=l2+1; A_block(l2) = i; A_block_Number(l2)=NO_number{i}; depA(l2)=(dd(i)-dd(A_k)); % the left part set.
        LR_A(l2)=NODES{i}.parent_child;
        for kk=1:depA(l2);
            PARENTa(1)=i;
            num1=PARENTa(kk);
            PARENTa(kk+1)=find([NO_number{:}]==floor(NO_number{num1}/2));
            LR_AALL{l2}(depA(l2)+1-kk)=NODES{PARENTa(kk)}.parent_child;% the right part set.
        end
    end
    end
end
% find the set of nodes which are in the B
r2=0;
for j=1:a_NODES
    if dd(j)>dd(B_k)
    if (2^(dd(j)-dd(B_k))*(NO_number{B_k}) <=NO_number{j} && 2^(dd(j)-dd(B_k))*(NO_number{B_k}+1) >NO_number{j}) 
        r2=r2+1; B_block(r2) = j; B_block_Number(r2)=NO_number{j}; depB(r2)=(dd(j)-dd(B_k));
        LR_B(r2)=NODES{j}.parent_child;% the right part set.
        %[dd(i) B_block]
    end
    end
end  
% find the set of nodes which are in the C
c2=0;
for j=1:a_NODES
    if dd(j)>dd(C_k)
    if (2^(dd(j)-dd(C_k))*(NO_number{C_k}) <=NO_number{j} && 2^(dd(j)-dd(C_k))*(NO_number{C_k}+1) >NO_number{j}) 
        c2=c2+1; C_block(c2) = j; C_block_Number(c2)=NO_number{j}; depC(c2)=(dd(j)-dd(C_k)); LR_C(c2)=NODES{j}.parent_child;
        for kk=1:depC(c2);
            PARENTc(1)=j;
            num1=PARENTc(kk);
            PARENTc(kk+1)=find([NO_number{:}]==floor(NO_number{num1}/2));%floor(num1/2)
            LR_CALL{c2}(depC(c2)+1-kk)=NODES{PARENTc(kk)}.parent_child;% the right part set.
        end
    %C_block
    end
    end
end  

% extend the numbers in the new Nodes  
% mod(k_CH,2)
NO_B=[];
NO_A=[];
NO_C=[];

NO_B_N=[];
NO_A_N=[];
NO_C_N=[];

if (mod(K_CH,2)==0) % right rotation same as (NODES{k_CH}.parent_child==0)%
    if l2>0
        for ii=1:l2
            %NODES_new{A_block(ii)}.number=NO_number{A_block(ii)}-2^(dd(A_k)-1)*2^(depA(ii));
            sum_LR=0;
            for kk=1:depA(ii)
                sum_LR=sum_LR+LR_AALL{ii}(kk)*2^(depA(ii)-kk);
            end
            PALL=find([NO_number{:}]==floor((NO_number{A_k}/2)));
            NODES_new{A_block(ii)}.number=2^(depA(ii))*NO_number{PALL}+sum_LR;

            NO_A_N(ii)=NODES_new{A_block(ii)}.number;
            NO_A(ii)=NODES{A_block(ii)}.number;
        end
    end
    % reduce 
    if r2>0
    for jj=1:r2
        NODES_new{B_block(jj)}.number=NO_number{B_block(jj)}+2^(depB(jj));
        NO_B_N(jj)=NODES_new{B_block(jj)}.number;
        NO_B(jj)=NODES{B_block(jj)}.number;        
    end
    end
    % change the place of the the right and make it left.
    if c2>0
    for jj=1:c2
         %NODES_new{C_block(jj)}.number=NO_number{C_block(jj)}+2^(dd(C_k))*2^(depC(jj));
        sum_LR=0;
        for kk=1:depC(jj)
            sum_LR=sum_LR+LR_CALL{jj}(kk)*2^(depC(jj)-kk);
        end
        PALL=find([NO_number{:}]==floor((NO_number{C_k}*2+1)));
       NODES_new{C_block(jj)}.number=2^(depC(jj))*NO_number{PALL}+sum_LR;
        NO_C_N(jj)=NODES_new{C_block(jj)}.number;
        NO_C(jj)=NODES{C_block(jj)}.number;
    end
    end
else % left rotation
    %increase this
    if l2>0
        for ii=1:l2
            %NODES_new{A_block(ii)}.number=NO_number{A_block(ii)}+2^(dd(A_k))*2^(depA(ii));
        sum_LR=0;
        for kk=1:depA(ii)
            sum_LR=sum_LR+LR_AALL{ii}(kk)*2^(depA(ii)-kk);
        end
        PALL=find([NO_number{:}]==floor((NO_number{A_k}*2)));
        NODES_new{A_block(ii)}.number=2^(depA(ii))*NO_number{PALL}+sum_LR;            
            
            NO_A_N(ii)=NODES_new{A_block(ii)}.number;
            NO_A(ii)=NODES{A_block(ii)}.number;
        end
    end
    % reduce the 
    if c2>0
    for jj=1:c2
        %NODES_new{C_block(jj)}.number=NO_number{C_block(jj)}-2^(dd(C_k)-1)*2^(depC(jj));
        sum_LR=0;
        for kk=1:depC(jj)
            sum_LR=sum_LR+LR_CALL{jj}(kk)*2^(depC(jj)-kk);
        end
        PALL=find([NO_number{:}]==floor(NO_number{C_k}/2));
       NODES_new{C_block(jj)}.number=2^(depC(jj))*NO_number{PALL}+sum_LR;
        NO_C_N(jj)=NODES_new{C_block(jj)}.number;
        NO_C(jj)=NODES{C_block(jj)}.number;
    end
    end
        % change the place of the the right and make it left.
    if r2>0
    for jj=1:r2
        %NODES_new{B_block(jj)}=NODES{B_block(jj)};
        NODES_new{B_block(jj)}.number=NO_number{B_block(jj)}-2^(depB(jj));
        NO_B_N(jj)=NODES_new{B_block(jj)}.number;
        NO_B(jj)=NODES{B_block(jj)}.number;
    end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jj=1:a_NODES
    NEW_NO_number{jj} = NODES_new{jj}.number; NEW_NO(jj)=NODES_new{jj}.number;       
    A_rull(jj)=NODES{jj}.rul;   
    NEW_rull(jj)=NODES_new{jj}.rul; 
    NEW_dir(jj)=NODES_new{jj}.direction; 
    NEW_PA(jj)=NODES_new{jj}.parent_child;   
    NEW_INT(jj)=NODES_new{jj}.inter;   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change lim_ALL, Zt_new, site_new and nn with the new tree 
    lim_ALL_new{1} =lim_ALL{1}; LIM_ALL_NEW{1}=lim_ALL{1};
    Zt_new1_1{1} = Zt1{1}; ZT_NEW1{1}=Zt1{1};
    Zt_new2_1{1} = Zt2{1}; ZT_NEW2{1}=Zt2{1};
    
    site_new1_1{1}= site1{1};  SITE_NEW1{1}=site1{1};
    site_new2_1{1}= site2{1};  SITE_NEW2{1}=site2{1};
    
    site1_1 = site1{1}; site2_1 = site2{1};    
    Zt1_1=Zt1{1}; Zt2_1=Zt2{1};
    nn2{1}=nn{1};NN2{1}=nn2{1}; 

    h=1;
for num_1=sort(NEW_NO(2:a_NODES)) % exclude the first one
    %num_1=NODES_new{jj}.number;
    jj=find([NEW_NO_number{:}]==num_1);
    J=find([NEW_NO_number{:}]==floor(num_1/2)); % coresponding PARENT NODE.
    %J=[NEW_NO_number{:}]==floor(num_1/2); % coresponding PARENT NODE.
    %[jj J]
    lim_ALL_new{jj} = lim_ALL_new{J}; % Let say you have the parent node for the moment.    
    % assign the spliting rul.
    if  (NODES_new{jj}.parent_child == 0);%mod((num_1),2)==0 %
        lim_ALL_new{jj}(NODES_new{jj}.direction,2)= NODES_new{jj}.rul; %number is even (left child ) % change the uper bound
    else
        lim_ALL_new{jj}(NODES_new{jj}.direction,1)= NODES_new{jj}.rul; %number is odd (right child)
    end   
% Change the Zt_new and site_new variable.    
     [site_new1_1{jj}, Zt_new1_1{jj}, nn1{jj}]= RE_siteLim(site1_1,Zt1_1,lim_ALL_new{jj});
     [site_new2_1{jj}, Zt_new2_1{jj}, nn2{jj}]= RE_siteLim(site2_1,Zt2_1,lim_ALL_new{jj});
%%% Create new sorrted tree %%%
    % change tree 
    h=h+1;
    NN2{h}=nn2{jj};
    NODES_NEW{h}= NODES_new{jj};
    LIM_ALL_NEW{h}=lim_ALL_new{jj};
    
    ZT_NEW1{h}=Zt_new1_1{jj}; ZT_NEW2{h}=Zt_new2_1{jj};    
    SITE_NEW1{h}=site_new1_1{jj};  SITE_NEW2{h}=site_new2_1{jj};    
    
    % change parameters
    prPROP{h}=prProp{jj};
    SIGMA_etaNEW{h}=SIGMA_eta1{jj}; 
    AA_etaNEW{h}=AA_eta1{jj};
    
    SIGMA_deltaNEW{h}=SIGMA_delta1{jj}; 
    AA_deltaNEW{h}=AA_delta1{jj};
    tau2_Y_NEW{h}=tau2_Y_new{jj};
end
% Compute priors 
for i=1:a_NODES
        dd(i)=floor(log(NODES{i}.number)/log(2)); %depth
        if (NODES{i}.inter==1)
            P_split(i)=aalpha*(1+dd(i))^(-bbeta);
        else
            P_split(i)=1-aalpha*(1+dd(i))^(-bbeta);
        end
end
P_Tree_curr= prod(P_split); 
  for i=1:a_NODES
        dd(i)=floor(log(NODES_new{i}.number)/log(2)); %depth
        if (NODES_new{i}.inter==1)
            P_split_new(i)=aalpha*(1+dd(i))^(-bbeta);
        else
            P_split_new(i)=(1-aalpha*(1+dd(i))^(-bbeta));
        end
  end
  P_Tree_prop =  prod(P_split_new);
  MH_ratio =log(P_Tree_prop) -log(P_Tree_curr); % MH log ratio
  u = log(rand);
    if (u<MH_ratio)
    % accept
        NODES= NODES_NEW; 
        lim_ALL = LIM_ALL_NEW; 
        Zt1=ZT_NEW1; Zt2=ZT_NEW2;
        site1=SITE_NEW1; site2=SITE_NEW2; 
        nn =NN2; 
        
        prCurr=prPROP;
        SIGMA_eta = SIGMA_etaNEW;
        AA_eta = AA_etaNEW;
        SIGMA_delta = SIGMA_deltaNEW;
        AA_delta = AA_deltaNEW;
        tau2_Y =tau2_Y_NEW;
    end
end
end

