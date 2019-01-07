clear all
close all

% sgp12km ================================================================

% sgp12km_daily_rain_obs -------------------------------------------------

clear all ; close all ;

load ./timeaveraged_data/sgp12km_daily_rain_obs_mean.mat sgp12km_daily_rain_obs_mean J K
load ./timeaveraged_data/sgp12km_lat_lon2d.mat sgp12km_lat_lon2d J K

xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

y = sgp12km_daily_rain_obs_mean ;
x1 = squeeze(sgp12km_lat_lon2d(1,:,:)) ;
x2 = squeeze(sgp12km_lat_lon2d(2,:,:)) ;


[ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
    %
    colorbar('location','NorthOutside')
    %
    title('Percipitation (Experimental)')
    %
    set(gca_m,'FontSize',25);
    set(findall(gcf_m,'type','text'),'FontSize',25);
    %
    %orient landscape
    orient tall
    set(gca_m,'Position',[.05 .01 .9 .95])
    %
    thefile = 'sgp12km_daily_rain_obs_mean_plot' ;
    saveas(gcf_m,sprintf('./plots/%s.fig',thefile));
    saveas(gcf_m,sprintf('./plots/%s.pdf',thefile));

clear all ; 

% sgp12km_daily_rain_June -------------------------------------------------

clear all ; close all ;

load ./timeaveraged_data/sgp12km_daily_rain_June_mean.mat sgp12km_daily_rain_June_mean J K
load ./timeaveraged_data/sgp12km_lat_lon2d.mat sgp12km_lat_lon2d J K

xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

y = squeeze(mean(sgp12km_daily_rain_June_mean,1)) ;
x1 = squeeze(sgp12km_lat_lon2d(1,:,:)) ;
x2 = squeeze(sgp12km_lat_lon2d(2,:,:)) ;


[ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
    %
    colorbar('location','NorthOutside')
    %
    title('Percipitation (Simulations)')
    %
    set(gca_m,'FontSize',25);
    set(findall(gcf_m,'type','text'),'FontSize',25);
    %
    %orient landscape
    orient tall
    set(gca_m,'Position',[.05 .01 .9 .95])
    %
    thefile = 'sgp12km_daily_rain_June_mean_plot' ;
    saveas(gcf_m,sprintf('./plots/%s.fig',thefile));
    saveas(gcf_m,sprintf('./plots/%s.pdf',thefile));

clear all ; 

% sgp12km_daily_rain_June TRUNK ------------------------------------------

clear all ; close all ;

load ./timeaveraged_data/sgp12km_daily_rain_June_mean_trunc.mat sgp12km_daily_rain_June_mean J K
load ./timeaveraged_data/sgp12km_lat_lon2d.mat sgp12km_lat_lon2d J K

xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

y = squeeze(mean(sgp12km_daily_rain_June_mean,1)) ;
x1 = squeeze(sgp12km_lat_lon2d(1,:,:)) ;
x2 = squeeze(sgp12km_lat_lon2d(2,:,:)) ;


[ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
    %
    colorbar('location','NorthOutside')
    %
    title('Percipitation (Simulations)')
    %
    set(gca_m,'FontSize',25);
    set(findall(gcf_m,'type','text'),'FontSize',25);
    %
    %orient landscape
    orient tall
    set(gca_m,'Position',[.05 .01 .9 .95])
    %
    thefile = 'sgp12km_daily_rain_June_mean_trunc_plot' ;
    saveas(gcf_m,sprintf('./plots/%s.fig',thefile));
    saveas(gcf_m,sprintf('./plots/%s.pdf',thefile));

clear all ; 
    
% sgp25km ================================================================

% sgp25km_daily_rain_obs -------------------------------------------------

clear all ; close all ;

load ./timeaveraged_data/sgp25km_daily_rain_obs_mean.mat sgp25km_daily_rain_obs_mean J K
load ./timeaveraged_data/sgp25km_lat_lon2d.mat sgp25km_lat_lon2d J K

xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

y = sgp25km_daily_rain_obs_mean ;
x1 = squeeze(sgp25km_lat_lon2d(1,:,:)) ;
x2 = squeeze(sgp25km_lat_lon2d(2,:,:)) ;


[ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
    %
    colorbar('location','NorthOutside')
    %
    title('Percipitation (Experimental)')
    %
    set(gca_m,'FontSize',25);
    set(findall(gcf_m,'type','text'),'FontSize',25);
    %
    %orient landscape
    orient tall
    set(gca_m,'Position',[.05 .01 .9 .95])
    %
    thefile = 'sgp25km_daily_rain_obs_mean_plot' ;
    saveas(gcf_m,sprintf('./plots/%s.fig',thefile));
    saveas(gcf_m,sprintf('./plots/%s.pdf',thefile));
%
clear all ;

% sgp25km_daily_rain_June -------------------------------------------------

clear all ; close all ;

load ./timeaveraged_data/sgp25km_daily_rain_June_mean.mat sgp25km_daily_rain_June_mean J K
load ./timeaveraged_data/sgp25km_lat_lon2d.mat sgp25km_lat_lon2d J K

xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

y = squeeze(mean(sgp25km_daily_rain_June_mean,1)) ;
x1 = squeeze(sgp25km_lat_lon2d(1,:,:)) ;
x2 = squeeze(sgp25km_lat_lon2d(2,:,:)) ;


[ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
    %
    colorbar('location','NorthOutside')
    %
    title('Percipitation (Simulations)')
    %
    set(gca_m,'FontSize',25);
    set(findall(gcf_m,'type','text'),'FontSize',25);
    %
    %orient landscape
    orient tall
    set(gca_m,'Position',[.05 .01 .9 .95])
    %
    thefile = 'sgp25km_daily_rain_June_mean_plot' ;
    saveas(gcf_m,sprintf('./plots/%s.fig',thefile));
    saveas(gcf_m,sprintf('./plots/%s.pdf',thefile));
%
clear all

% sgp25km_daily_rain_June TRUNK ------------------------------------------

clear all ; close all ;

load ./timeaveraged_data/sgp25km_daily_rain_June_mean_trunc.mat sgp25km_daily_rain_June_mean J K
load ./timeaveraged_data/sgp25km_lat_lon2d.mat sgp25km_lat_lon2d J K

xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

y = squeeze(mean(sgp25km_daily_rain_June_mean,1)) ;
x1 = squeeze(sgp25km_lat_lon2d(1,:,:)) ;
x2 = squeeze(sgp25km_lat_lon2d(2,:,:)) ;


[ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
    %
    colorbar('location','NorthOutside')
    %
    title('Percipitation (Simulations)')
    %
    set(gca_m,'FontSize',25);
    set(findall(gcf_m,'type','text'),'FontSize',25);
    %
    %orient landscape
    orient tall
    set(gca_m,'Position',[.05 .01 .9 .95])
    %
    thefile = 'sgp25km_daily_rain_June_mean_trunc_plot' ;
    saveas(gcf_m,sprintf('./plots/%s.fig',thefile));
    saveas(gcf_m,sprintf('./plots/%s.pdf',thefile));

clear all ; 

    
% sgp50km ================================================================

% sgp50km_daily_rain_obs -------------------------------------------------

clear all ; close all ;

load ./timeaveraged_data/sgp50km_daily_rain_obs_mean.mat sgp50km_daily_rain_obs_mean J K
load ./timeaveraged_data/sgp50km_lat_lon2d.mat sgp50km_lat_lon2d J K

xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

y = sgp50km_daily_rain_obs_mean ;
x1 = squeeze(sgp50km_lat_lon2d(1,:,:)) ;
x2 = squeeze(sgp50km_lat_lon2d(2,:,:)) ;


[ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
    %
    colorbar('location','NorthOutside')
    %
    title('Percipitation (Experimental)')
    %
    set(gca_m,'FontSize',25);
    set(findall(gcf_m,'type','text'),'FontSize',25);
    %
    %orient landscape
    orient tall
    set(gca_m,'Position',[.05 .01 .9 .95])
    %
    thefile = 'sgp50km_daily_rain_obs_mean_plot' ;
    saveas(gcf_m,sprintf('./plots/%s.fig',thefile));
    saveas(gcf_m,sprintf('./plots/%s.pdf',thefile));
%
clear all

% sgp50km_daily_rain_June -------------------------------------------------

clear all ; close all ;

load ./timeaveraged_data/sgp50km_daily_rain_June_mean.mat sgp50km_daily_rain_June_mean J K
load ./timeaveraged_data/sgp50km_lat_lon2d.mat sgp50km_lat_lon2d J K

xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

y = squeeze(mean(sgp50km_daily_rain_June_mean,1)) ;
x1 = squeeze(sgp50km_lat_lon2d(1,:,:)) ;
x2 = squeeze(sgp50km_lat_lon2d(2,:,:)) ;


[ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
    %
    colorbar('location','NorthOutside')
    %
    title('Percipitation (Simulations)')
    %
    set(gca_m,'FontSize',25);
    set(findall(gcf_m,'type','text'),'FontSize',25);
    %
    %orient landscape
    orient tall
    set(gca_m,'Position',[.05 .01 .9 .95])
    %
    thefile = 'sgp50km_daily_rain_June_mean_plot' ;
    saveas(gcf_m,sprintf('./plots/%s.fig',thefile));
    saveas(gcf_m,sprintf('./plots/%s.pdf',thefile));
%
clear all

% sgp50km_daily_rain_June TRUNK ------------------------------------------

clear all ; close all ;

load ./timeaveraged_data/sgp50km_daily_rain_June_mean_trunc.mat sgp50km_daily_rain_June_mean J K
load ./timeaveraged_data/sgp50km_lat_lon2d.mat sgp50km_lat_lon2d J K

xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

y = squeeze(mean(sgp50km_daily_rain_June_mean,1)) ;
x1 = squeeze(sgp50km_lat_lon2d(1,:,:)) ;
x2 = squeeze(sgp50km_lat_lon2d(2,:,:)) ;


[ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
    %
    colorbar('location','NorthOutside')
    %
    title('Percipitation (Simulations)')
    %
    set(gca_m,'FontSize',25);
    set(findall(gcf_m,'type','text'),'FontSize',25);
    %
    %orient landscape
    orient tall
    set(gca_m,'Position',[.05 .01 .9 .95])
    %
    thefile = 'sgp50km_daily_rain_June_mean_trunc_plot' ;
    saveas(gcf_m,sprintf('./plots/%s.fig',thefile));
    saveas(gcf_m,sprintf('./plots/%s.pdf',thefile));

clear all ; 
    
% cam25km ================================================================

% cam25km_daily_rain_obs -------------------------------------------------

clear all

load ./timeaveraged_data/cam25km_daily_rain_obs_mean.mat cam25km_daily_rain_obs_mean J K
load ./timeaveraged_data/cam25km_lat_lon2d.mat cam25km_lat_lon2d J K

xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

y = cam25km_daily_rain_obs_mean ;
x1 = squeeze(cam25km_lat_lon2d(1,:,:)) ;
x2 = squeeze(cam25km_lat_lon2d(2,:,:)) ;


[ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
    %
    colorbar('location','NorthOutside')
    %
    title('Percipitation (Experimental)')
    %
    set(gca_m,'FontSize',25);
    set(findall(gcf_m,'type','text'),'FontSize',25);
    %
    %orient landscape
    orient tall
    set(gca_m,'Position',[.05 .01 .9 .95])
    %
    thefile = 'cam25km_daily_rain_obs_mean_plot' ;
    saveas(gcf_m,sprintf('./plots/%s.fig',thefile));
    saveas(gcf_m,sprintf('./plots/%s.pdf',thefile));
%
clear all

% cam25km_daily_rain_June -------------------------------------------------

clear all ; close all

load ./timeaveraged_data/cam25km_daily_rain_June_mean.mat cam25km_daily_rain_June_mean J K
load ./timeaveraged_data/cam25km_lat_lon2d.mat cam25km_lat_lon2d J K


xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

y = squeeze(mean(cam25km_daily_rain_June_mean(:,:,:),1)) ;
x1 = squeeze(cam25km_lat_lon2d(1,:,:)) ;
x2 = squeeze(cam25km_lat_lon2d(2,:,:)) ;


[ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
    %
    colorbar('location','NorthOutside')
    %
    title('Percipitation (Simulations)')
    %
    set(gca_m,'FontSize',25);
    set(findall(gcf_m,'type','text'),'FontSize',25);
    %
    %orient landscape
    orient tall
    set(gca_m,'Position',[.05 .01 .9 .95])
    %
    thefile = 'cam25km_daily_rain_June_mean_plot' ;
    saveas(gcf_m,sprintf('./plots/%s.fig',thefile));
    saveas(gcf_m,sprintf('./plots/%s.pdf',thefile));
%
clear all

% cam25km_daily_rain_June TRUNK ------------------------------------------

clear all ; close all ;

load ./timeaveraged_data/cam25km_daily_rain_June_mean_trunc.mat cam25km_daily_rain_June_mean J K
load ./timeaveraged_data/cam25km_lat_lon2d.mat cam25km_lat_lon2d J K

xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

y = squeeze(mean(cam25km_daily_rain_June_mean,1)) ;
x1 = squeeze(cam25km_lat_lon2d(1,:,:)) ;
x2 = squeeze(cam25km_lat_lon2d(2,:,:)) ;


[ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
    %
    colorbar('location','NorthOutside')
    %
    title('Percipitation (Simulations)')
    %
    set(gca_m,'FontSize',25);
    set(findall(gcf_m,'type','text'),'FontSize',25);
    %
    %orient landscape
    orient tall
    set(gca_m,'Position',[.05 .01 .9 .95])
    %
    thefile = 'cam25km_daily_rain_June_mean_trunc_plot' ;
    saveas(gcf_m,sprintf('./plots/%s.fig',thefile));
    saveas(gcf_m,sprintf('./plots/%s.pdf',thefile));

clear all ; 
    
    
    
    
    
    
