function [Rho_mat,site_W]=Rho_Tree(rho,site_e,a_new,S)


for k=1:a_new
[Rho_mati{k}]=ones(size(site_e{k},1),1).*(rho{k}{S});

end

         [Rho_mat]=[Rho_mati{1}];
         [site_W]=[site_e{1}];
         for i=2:a_new
             [Rho_mat]=[Rho_mat; Rho_mati{i}];
             [site_W]=[site_W ; site_e{i}];
         end   

end