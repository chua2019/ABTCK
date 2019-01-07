function [site,Zt,site_real,Zt_real, nn, lim_ALL, prCurr,NODES]= ...
    Change_Non_nest(prCurr,site,Zt,site_real,Zt_real,nn,AA,S,low_l,lim_ALL, a_NODES,NODES,tau1,tau2,funname)
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
% NODES{i}.direction % here the direction we cut it .
% NODES{i}.rul  % here the rule we use --> a value between 0 and 1.

lim_ALL_new =lim_ALL;
prProp=prCurr;
NODES_new=NODES;
% find all the interior nodes if you want the overall.

for jj=1:a_NODES
NO_number{jj} = NODES{jj}.number;
end
bb=0;
for  ii=1:a_NODES
    if (NODES{ii}.inter == 1)
            bb=bb+1; PP(bb)=ii;
    end
end

kk1=randi(bb); k = PP(kk1); 
K = NODES{k}.number;
% find the left and right child to change the rull
        k_LP=find([NO_number{:}]==2*K); % Left child
        k_RP=find([NO_number{:}]==(2*K+1)); % Right child
% Find the the min of the max and the max of the min of the particular direction .

rull_k=NODES{k_LP}.rul;
direk = NODES{k_LP}.direction;

bb1=0; bb2=0;
rull_MIN_ALL(1)= lim_ALL{k}(direk,1);
rull_MAX_ALL(1)= lim_ALL{k}(direk,2);
for  ii=1:a_NODES
    if  (NODES{ii}.direction==direk)
        if (NODES{ii}.rul<rull_k)
            bb1=bb1+1; rull_MIN_ALL(bb1+1)=NODES{ii}.rul; % 
        end
        if (NODES{ii}.rul>rull_k)
            bb2=bb2+1; rull_MAX_ALL(bb2+1)=NODES{ii}.rul; % 
        end        
    end
end

LowerBound = max(rull_MIN_ALL);
UpperBound = min(rull_MAX_ALL);
% New rull
 New_rul=tnormrnd(rull_k,(UpperBound-LowerBound)/6,LowerBound,UpperBound);

 NODES_new{k_LP}.rul = New_rul;
 NODES_new{k_RP}.rul = New_rul;
 
% Rebuild the limits. 
    lim_ALL_new{1} =lim_ALL{1};
    nn1{1}=nn{1};  
    for t=1:S
    n_Axri{1}{t}=size(site_real{1}{t},1);
    end
    Zt_new1_1{1} = Zt{1}; site_new1_1{1}= site{1};  
    Zt_real_new{1} = Zt_real{1}; site_real_new{1}= site_real{1}; 
    site1_1 = site{1}; Zt1_1=Zt{1};
    site_real_1 = site_real{1}; Zt_real_1=Zt_real{1};

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
     
     [site_real_new{jj}{t}, Zt_real_new{jj}{t}, n_Axri{jj}{t}]= RE_siteLim(site_real_1{t},Zt_real_1{t},lim_ALL_new{jj});  
end
site_ne{jj}=site_new1_1{jj};
end


%%% Include both the restrictions?
for jj=1:a_NODES
 Min_n1(jj)=min(n_Axri{jj}{:});
end
 Min_n=min(Min_n1);
%[Min_n, rull_k]

     if  (Min_n>low_l)%((min([nn2{:}])>4))  %((min([nn1{:}])>8) && (min([nn2{:}])>15))% 
         Dim_y = size(Zt_new1_1{1},2);
    % posterior for the changed nodes ()
        % Find the external NODES so we can make the GP updates.  
%     hh2=0;
%     for ii=1:a_NODES
%         if (NODES{ii}.inter==0)
%             hh2=hh2+1; fin_nodes(hh2)=NODES{ii};
%             depth_each(hh2)=floor(NODES{ii}.number/2);   
%             
%             prCurr_ex(hh2)=0; 
%              %XX_s{ii}=mean_basisMulti(site{ii},Zt{ii},S);
%             for t=1:S
%             %[mu_Zt11, SIGMA11, beta_X11, prCurr{ii}{t}]=Multi_mu_s(Zt{ii}{t},AA{ii}{t}, site{ii}{t},XX_s{ii}{t});
%             prCurr_ex(hh2)=prCurr_ex(hh2)+(prCurr{ii}{t});
%             end
%         end
%     end    
   % Find the external NODES_new so we can make the GP updates.  
    hh2=0;
    for ii=1:a_NODES
        if (NODES_new{ii}.inter==0)
            hh2=hh2+1; fin_nodes(hh2)=NODES_new{ii};
            depth_each(hh2)=floor(NODES_new{ii}.number/2);            
            Zt_new_ex1{hh2}= Zt_new1_1{ii};  site_new_ex1{hh2}= site_new1_1{ii};  
            AA_ex1{hh2}= AA{ii};    
            %lim_ALL_ex{hh2}=lim_ALL{ii};                      

             prCurr_ex(hh2)=0; 
             %XX_s{ii}=mean_basisMulti(site{ii},Zt{ii},S);
            for t=1:S
            %[mu_Zt11, SIGMA11, beta_X11, prCurr{ii}{t}]=Multi_mu_s(Zt{ii}{t},AA{ii}{t}, site{ii}{t},XX_s{ii}{t});
            prCurr_ex(hh2)=prCurr_ex(hh2)+(prCurr{ii}{t});
            end
            
        % posterior for the changed node k_LP 
        X_s{hh2}=mean_basisMulti(site_new_ex1{hh2},Zt_new_ex1{hh2},S); 
        prProp_ex(hh2)=0; 
        for t=1:S
           % prProp_ex11{t}=0;
        %[mu_Zt11, SIGMA11, beta_X11, prProp_ex11{t}]=Multi_mu_s(Zt_new_ex1{hh2}{t},AA_ex1{hh2}{t}, site_new_ex1{hh2}{t},X_s{hh2}{t});
        %prProp{ii}{t}=(prProp_ex11{t}); 
        %prProp_ex(hh2)=prProp_ex(hh2)+prProp_ex11{t};
        [mu_Zt11, SIGMA11, beta_X11, prProp{ii}{t}]=Multi_mu_s3(Zt_new_ex1{hh2}{t},AA_ex1{hh2}{t}, site_new_ex1{hh2}{t},X_s{hh2}{t},tau1,tau2);
        prProp_ex(hh2)=prProp_ex(hh2)+prProp{ii}{t};
        
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
%[MH_ratio, u];
if (u<MH_ratio)
% accept
    site=site_new1_1; Zt=Zt_new1_1;
    site_real=site_real_new; Zt_real=Zt_real_new; 
    nn=nn1; prCurr = prProp; lim_ALL = lim_ALL_new;
    NODES = NODES_new;
    
end
     end  
end
