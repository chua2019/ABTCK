function [Zt,site,prCurr,AA,a_new,lim_all, KrigS,Var_krigS,Krig_rho, Krig_W_etaMAT, Krig_site_eta]=...
    BTMulti(PR_CURR,SITE,ZT,n,AA_E,LIM_ALL,a_new,a_NODES,NODES,pprop, aalpha, bbeta,S,low_l,site_e,n_e,BB,tau1,tau2,Delta,funname)
KrigS.mean= zeros(BB,n_e); KrigS.Low= zeros(BB,n_e);
KrigS.Up= zeros(BB,n_e); Var_krigS= zeros(BB,n_e); 
Krig_site_eta=[]; Krig_W_etaMAT=zeros(BB,size(site_e,1)); Krig_rho{1}=zeros(BB,size(site_e,1));   Krig_rho{2}=zeros(BB,size(site_e,1)); 
Dim_x =size(SITE{1}{1},2); mm=sqrt(size(site_e,1))-1; 

for i=1:BB
    % ZT_NEW, SITE_NEW,  LIM_ALL, PR_CURR FOR ALL THE NODES. 
    % Zt_new AND site_new FOR EXTERNAL NODES.
    operation = randi(5);
    if (operation ==1)
    % Split 
    [SITE,ZT, n, PR_CURR, AA_E,LIM_ALL, a_new,NODES,a_NODES]= Split_Multi(PR_CURR,SITE,ZT,n,AA_E,S,low_l,LIM_ALL,a_new,a_NODES,NODES, pprop, aalpha, bbeta,tau1,tau2,funname); 
    elseif  (operation ==2)
            % Merge 
    if a_NODES>2 
        [SITE,ZT, n, PR_CURR, AA_E,a_new, LIM_ALL,NODES,a_NODES]=Merge_Multi(PR_CURR,SITE,ZT,n,AA_E,S,LIM_ALL,a_new,a_NODES,NODES, aalpha, bbeta,tau1,tau2,funname);   
    end
    elseif  (operation >2 && operation < 5)
    % Change 
    if a_NODES>2     
               
                                                          
        [SITE,ZT, n, LIM_ALL, PR_CURR,NODES]= Change_Multi(PR_CURR,SITE,ZT,n,AA_E,S,low_l,LIM_ALL, a_NODES,NODES,tau1,tau2,funname);

    end
    elseif  (operation ==5)
        %Swap-Rotate
    if a_NODES>2
        suboperation = randi(2);
        if (suboperation ==1)
            % Swap
            [SITE,ZT, n, LIM_ALL, PR_CURR,NODES]= Swap_Multi(PR_CURR,SITE,ZT,n,AA_E,S,low_l,LIM_ALL, a_NODES,NODES,tau1,tau2,funname);
        else
            %Rotate

        end
    end
    end 

   for inter=1:a_NODES

   for t=1:S
    [YY,II]=sort(SITE{inter}{t}(:,1));
    SITE{inter}{t}=SITE{inter}{t}(II,:); 
    
    ZT{inter}{t}=ZT{inter}{t}(II,:); 
   end 
   end
    
    
    % Find the external NODES so we can make the GP updates.  
   clear Zt; clear AA; clear beta_X; clear site; clear lim_all;
    hh2=0; 
    
    for ii=1:a_NODES
        if (NODES{ii}.inter==0)
            hh2=hh2+1; fin_nodes(hh2)=NODES{ii};
            depth_each(hh2)=floor(NODES{ii}.number/2);
            prCurr{hh2}= PR_CURR{ii}; Zt{hh2}= ZT{ii};           
            site{hh2}=SITE{ii}; 
            AA{hh2}= AA_E{ii}; 
            lim_all{hh2}=LIM_ALL{ii}; 
            %nn{hh2}=n{ii};
            % GP update  
           for jj=1:1
                % GP update 
                [AA{hh2}, sigma2{hh2}, prCurr{hh2}, mu{hh2}, beta_X{hh2}]= Multi_MH_s(Zt{hh2},site{hh2},AA{hh2},S,Delta,tau1,tau2);
           end
            PR_CURR{ii}=prCurr{hh2};  AA_E{ii}=AA{hh2};  LIM_ALL{ii}=lim_all{hh2}; 
            beta_X{hh2}{1}(2)=1;
            for t=1:S 
            beta{hh2}{t}=beta_X{hh2}{t}(1);
            rho{hh2}{t}=beta_X{hh2}{t}(2);
            %[hh2 beta_X{hh2}{t}(1) beta_X{hh2}{t}(2)]
            end
            
        end
    end     
       



M=1;

[site_e_t]= RE_site(site_e,lim_all,a_new);    
[Kriging_W, site_W, Var_krigW]=Krig_Multi_Tree(Zt,beta,rho,sigma2,AA,site,site_e_t,a_new,tau1,tau2,S);

XY_eta= [site_W, Kriging_W, Var_krigW]; 
XY_11_eta=sortrows(XY_eta,2);
XY_12_eta=sortrows(XY_11_eta,1);
Kriging_W_etaS=XY_12_eta(:,(Dim_x+1):(Dim_x+M));
Krig_W_etaMAT(i,:)=Kriging_W_etaS;%
Var_krigS(i,:)=XY_12_eta(:,(Dim_x+M+1));%

for kk=2:S
 [Rho_mat,site_W2]=Rho_Tree(rho,site_e_t,a_new,kk);
XY_rho= [site_W2, Rho_mat]; 
XY_11_rho=sortrows(XY_rho,2);
XY_12_rho=sortrows(XY_11_rho,1);
Krig_rho{kk}(i,:)=XY_12_rho(:,(Dim_x+1):(Dim_x+M));
end


Krig_site_eta=XY_12_eta(:,1:Dim_x);
    
%figure
if mod(i,100)<=0
    
    if (Dim_x==2)
XX1=[]; XX2=[]; KW_U=[];
for jj=1:(mm+1)
     KW_U=[KW_U Kriging_W_etaS((1+(jj-1)*(mm+1)):(jj*(mm+1)))];%KW_U=[KW_U Krig_rho((1+(jj-1)*(mm+1)):(jj*(mm+1)))];%
     XX1=[XX1 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),1)];
     XX2=[XX2 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),Dim_x)];    
end

imagesc(XX1(:),XX2(:),KW_U); colorbar; set(gca,'YDir','normal'); %caxis([-1.5 3]);     
 hold on;
for ii=1:a_new
    l1=lim_all{ii}(1,1); l2=lim_all{ii}(1,2);        
    %l3=lim_all{ii}((Dim_x-p),1); l4=lim_all{ii}((Dim_x-p),2);
    l3=lim_all{ii}((Dim_x),1); l4=lim_all{ii}((Dim_x),2);
%text(xcurr(1,i), xcurr(2,i),num2str(i), )
text((l1+l2)/2,(l3+l4)/2,num2str(ii), 'Color','r');
    plot([l1 l1],[l3 l4]);hold on; plot([l2 l2],[l3 l4]);    
    plot([l1 l2], [l3 l3]); plot([l1 l2],[l4 l4]);
end    
hold off; pause(.0001); % <- a time consuming OP
    end
fprintf(1,'iter %d\n',i);

end    
 





end

end
% hold off; 