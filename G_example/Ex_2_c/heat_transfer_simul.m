


% The partial differential equation for transient conduction heat transfer
% is:

% $$
% \rho C_p \frac{\partial T}{\partial t} -\nabla \cdot (k \nabla T) = f
% $$

% where $T$ is the temperature, 
% $\rho$ is the material density, 
% $C_p$  is the specific heat, 
% $k$ is the thermal conductivity. 
% $f$ is the heat generated inside the body which is zero in this example.

% Corse model
%
% R = heat_transfer_simul(0.5, 0, 0)
% u_xy = interpolateSolution(heat_transfer_simul(0.5, 0, 0),x,y) ;

% Intermediate model
% 
% R = heat_transfer_simul(0.02, 1, 0)
% u_xy = interpolateSolution(heat_transfer_simul(0.02, 1, 0),x,y) ;


% Fine model
% 
% R = heat_transfer_simul(0.015, 2, 0)
% u_xy = interpolateSolution(heat_transfer_simul(0.015, 2, 0),x,y) ;

% in order to get the output at (x,y), you do this
% u_xy = interpolateSolution(R,x,y) ;


function R = heat_transfer_simul(Hsize, QCont, Qplot)

% Hsize From 0 to 0.5; The step of the grid
% QCont =0 constant contactivity; =1 variable contactivity



%% Create a PDE Model with a single dependent variable
numberOfPDE = 1;
pdem = createpde(numberOfPDE);

%% Geometry
r1 = [3 4 0 1 1 0 0 0 3 3];
r2 = [3 4 0.5 0.5001 0.5001 0.5  1 1 2.5 2.5];
gdm = [r1; r2]';
g = decsg(gdm,'R1-R2',['R1'; 'R2']');

%% Convert the decsg format into a geometry object 

geometryFromEdges(pdem,g);


%% PDE Coefficients and Boundary Conditions
% The rest are Neuman du/d(whatever...) = 0

uRight = applyBoundaryCondition(pdem,'neumann','Edge',1,'g',-20);
uLeft = applyBoundaryCondition(pdem,'dirichlet','Edge',6,'u',100);


%% Coefficient Definition
% Elliptic pde kind of thing ...
% Check the conditions, ktl ...
rho = 1; 
Cp = 1; 
 
if QCont == 0
    k = @(region,state)exp(1.5*sin(10*pi*(region.y)/3).*(region.y<=0)); 
elseif QCont == 1
    k = @(region,state)exp(1.5*sin(10*pi*(region.y)/3).*(region.y<1.8)); 
elseif QCont == 2
    k = @(region,state)exp(1.5*sin(10*pi*(region.y)/3).*(region.y<=3)); 
else
    k = 0 ;
end


c = k;
a = 0;
f = 1; % a few thin, previous heat ...
d = rho*Cp; 
specifyCoefficients(pdem,'m',0,'d',0,'c',c,'a',a,'f',f);

%% Mesh
% corse hmax = 0.3 ;
% fine hmax = 0.02 ;
hmax = Hsize; 
msh=generateMesh(pdem,'Hmax',hmax);
% figure 
% pdeplot(pdem); 
% axis equal
% title 'Block With Finite Element Mesh Displayed'

%% Steady State Solution

R = solvepde(pdem);
    
%% Plots

% Plot the Contour
if Qplot
    figure
    pdeplot(pdem,'XYData',R.NodalSolution,'Contour','on','ColorMap','jet'); 
    axis equal
    title 'Temperature, Steady State Solution'
    orient landscape
  caxis([53, 100])
      set(gcf, 'PaperPosition', [0 0 25 25]); 
    set(gcf, 'PaperSize', [25 25]); 
    set(gca,'FontSize',20)
    saveas(gcf,sprintf('plotcont_%d.pdf', QCont))
    saveas(gcf,sprintf('plotcont_%d.fig', QCont))
end

% Plot the 3d plot
if Qplot
    figure;
    NN = 100 ;
    x_1 = linspace(0,1,NN) ;
    x_2 = linspace(0,3,NN) ;
    u_x1y2 = nan(NN,NN) ;
    for i = 1:NN
        for j = 1:NN
            u_x1y2(i,j) = interpolateSolution(R,x_1(i),x_2(j)) ;
        end
    end
    meshc(x_2,x_1,u_x1y2)
    title 'Temperature, Steady State Solution'
    orient landscape
    xlabel('$x_2$','interpreter','latex')
    ylabel('$x_1$','interpreter','latex')
    zlabel('$u(x_1,x_2)$','interpreter','latex')
    
    set(gcf, 'PaperPosition', [0 0 25 25]); 
    set(gcf, 'PaperSize', [25 25]); 
    set(gca,'FontSize',20)
    
    zlim([53, 100])

    if 1
    view(130, 30);
    saveas(gcf,sprintf('plot3d_%d.pdf', QCont))
    saveas(gcf,sprintf('plot3d_%d.fig', QCont))
    end
end

end


