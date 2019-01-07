function [site,Zt, nn, lim_ALL, prCurr,NODES]= ...
    Swap_Multi(prCurr,site,Zt,nn,AA,S,low_l,lim_ALL, a_NODES,NODES,tau1,tau2,funname)

% Add Separate the 2D area into smaller areas in the y direc. 
% BB = the subregion which will be splited.

% Find a parent NODE and change the spliting rull for his childrens.
% we have also to change all the limits of the decendant of this node, site and Zt (observation).
% Do not change the parameters.
% In THIS VERSION DO CHANGE ONLY THE PARENTS WITH EXTERNAL NODES.

% a : Nuber of external Nodes.
% a_NODES: Number of Nodes.
% NODES: all nodes of a tree - if we know this we can find the tree.
% lim_ALL: for every node we have a rectangular region in the input space.

% NODES{i}.number; % a node numbe k with have (2*k) and (2*k+1) node childs
% Nodes{i}.inter = 0 or 1 : if interior 1 if exterior 0.
% NODES{i}.parent_child % here whether or not is a left or right child (0 if left 1 if right).
% NODES{i}.direction % here the direction we cut it ([0] if it is child).
% NODES{i}.rul  % here the rule we use --> a value between 0 and 1.

lim_ALL_new =lim_ALL;
prProp=prCurr;
NODES_new = NODES;


%AA_eta1 = AA;

% find all the interior nodes and their coresponding parents (exclude node 1) if they spli in one direction.
for jj=1:a_NODES
    NO_number{jj} = NODES{jj}.number;
end

bb=0;
for  ii=2:a_NODES % without the initial node
    if (NODES{ii}.inter == 1)
            num_1=NODES{ii}.number;
            PARENT=find([NO_number{:}]==floor(num_1/2)); % coresponding PARENT NODE
            Children=find([NO_number{:}]==num_1*2); % coresponding PARENT NODE
            if (NODES{Children}.direction ~=NODES{ii}.direction)
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
K_P = NODES{k_P}.number; % Parent number
%[k_P, K_CH]
if (NODES{k_CH}.parent_child==0)%mod(K_CH,K_P)==0
    K_OCH = K_CH+1;%other_child = right child
else 
    K_OCH = K_CH-1;% other_child = left child 
end
k_OCH=find([NO_number{:}]==floor(K_OCH));
%[k_CH K_OCH]

% find the left and right child to change the rull
        k_LCH=find([NO_number{:}]==2*K_CH); % Left child of the child inintrior nodes 
        k_RCH=find([NO_number{:}]==(2*K_CH+1)); % Right child  of the child inintrior nodes 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NODES_new{k_CH}.parent_child  = NODES{k_CH}.parent_child; 
NODES_new{k_OCH}.parent_child   = NODES{k_OCH}.parent_child;
% change the rul
NODES_new{k_CH}.rul = NODES{k_LCH}.rul; 
NODES_new{k_OCH}.rul = NODES{k_LCH}.rul; 
% change the interior direction
NODES_new{k_CH}.direction = NODES{k_LCH}.direction; 
NODES_new{k_OCH}.direction = NODES{k_RCH}.direction;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change the interior direction, rull, number
NODES_new{k_LCH}.direction = NODES{k_CH}.direction;  % this depends from left or right roation (Do first the left.)
NODES_new{k_LCH}.rul = NODES{k_CH}.rul; 

NODES_new{k_RCH}.direction = NODES{k_CH}.direction;  % this depends from left or right roation (Do first the left.)
NODES_new{k_RCH}.rul = NODES{k_CH}.rul; 

for jj=1:a_NODES
    NEW_NO_number{jj} = NODES_new{jj}.number;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rebuild the limits. 
    lim_ALL_new{1} =lim_ALL{1};
    nn2{1}=nn{1}; nn1{1}=nn{1};
    Zt_new1_1{1} = Zt{1};   
    site_new1_1{1}= site{1};     
    site1_1 = site{1};
    Zt1_1=Zt{1};
    
    
for jj=2:a_NODES
    num_1=NODES_new{jj}.number;
    J=find([NO_number{:}]==floor(num_1/2)); % coresponding PARENT NODE.
    lim_ALL_new{jj} = lim_ALL_new{J}; % Let say you have the parent node for the moment.
    % assign the spliting rul.
    if  (NODES_new{jj}.parent_child == 0); % mod((num_1/2),2)==1 %
        lim_ALL_new{jj}(NODES_new{jj}.direction,2)= NODES_new{jj}.rul; %number is even (left child ) % change the uper bound
    else
        lim_ALL_new{jj}(NODES_new{jj}.direction,1)= NODES_new{jj}.rul; %number is odd (right child)
    end   
% Change the Zt_new and site_new variable.   
for t=1:S
     [site_new1_1{jj}{t}, Zt_new1_1{jj}{t}, nn1{jj}{t}]= RE_siteLim(site1_1{t},Zt1_1{t},lim_ALL_new{jj});  
end
site_ne{jj}=site_new1_1{jj};
end 


%%% Include both the restrictions?
for jj=1:a_NODES
 Min_n1(jj)=min(cell2mat(nn1{jj}));
end
 Min_n=min(Min_n1);

     if  (Min_n>low_l)%((min([nn2{:}])>4))  %((min([nn1{:}])>8) && (min([nn2{:}])>15))% 
         Dim_y = size(Zt_new1_1{1},2);
    % posterior for the changed nodes ()
        % Find the external NODES so we can make the GP updates.  
    hh2=0;
    for ii=1:a_NODES
        if (NODES{ii}.inter==0)
            hh2=hh2+1; fin_nodes(hh2)=NODES{ii};
            depth_each(hh2)=floor(NODES{ii}.number/2);   
            
            prCurr_ex(hh2)=0; 
            
            for t=1:S
            prCurr_ex(hh2)=prCurr_ex(hh2)+(prCurr{ii}{t});
            end
        end
    end    
   % Find the external NODES_new so we can make the GP updates.  
    hh2=0;
    for ii=1:a_NODES
        if (NODES_new{ii}.inter==0)
            hh2=hh2+1; fin_nodes(hh2)=NODES_new{ii};
            depth_each(hh2)=floor(NODES_new{ii}.number/2);
            
            Zt_new_ex1{hh2}= Zt_new1_1{ii};            
            site_new_ex1{hh2}= site_new1_1{ii};            
            AA_ex1{hh2}= AA{ii};    
            %lim_ALL_ex{hh2}=lim_ALL{ii};                      

        % posterior for the changed node k_LP 
        X_s{hh2}=mean_basisMulti(site_new_ex1{hh2},Zt_new_ex1{hh2},S); 
        prProp_ex(hh2)=0; 
        for t=1:S
        [mu_Zt11, SIGMA11, beta_X11, prProp{ii}{t}]=Multi_mu_s3(Zt_new_ex1{hh2}{t},AA_ex1{hh2}{t}, site_new_ex1{hh2}{t},X_s{hh2}{t},tau1,tau2);
        prProp_ex(hh2)=prProp_ex(hh2)+prProp{ii}{t};%prProp_ex11{t};   %prProp{ii}{t}=(prProp_ex11{t}); 
        end
            
        end
    end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Check
for jj=1:a_NODES
    NEW_NO_number1{jj} = NODES_new{jj}.number;   
end
  % MH
MH_ratio =(sum(prProp_ex)-sum(prCurr_ex));
u = log(rand);
if (u<MH_ratio)
% accept
    site=site_new1_1; 
    Zt=Zt_new1_1; 
    nn=nn1; prCurr = prProp; lim_ALL = lim_ALL_new;
    NODES = NODES_new;
    
end
     end  
end

%end


