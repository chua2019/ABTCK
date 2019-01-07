clear all
close all

x1_min = 25.5 ;
x1_max = 41 ;
x2_min = -109.3 ;
x2_max = -93.31 ;
    
xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

% sgp12km ================================================================

% sgp12km_daily_rain -------------------------------------------------

R = 150;

load ./sgp12km_daily_rain_obs_mean.mat sgp12km_daily_rain_obs_mean J K
load ./sgp12km_daily_rain_June_mean.mat sgp12km_daily_rain_June_mean
load ./sgp12km_lat_lon2d.mat sgp12km_lat_lon2d
load ./sgp12km_calibr.mat sgp12km_calibr

y = sgp12km_daily_rain_obs_mean ;
eta = sgp12km_daily_rain_June_mean ;
x = sgp12km_lat_lon2d ;
x1 = squeeze(x(1,:,:)) ;
x2 = squeeze(x(2,:,:)) ;

ind = find( ...
        (isnan(y(:))==0) ...
        .*(x1(:) > x1_min) ...
        .*(x1(:) < x1_max) ...
        .*(x2(:) > x2_min) ...
        .*(x2(:) < x2_max) ...
        ) ;
    
for r = 1:R
    ee = squeeze(eta(r,:,:)) ;
    aa = abs(ee(ind)-y(ind))/mean(y(ind)) ;
    DH(r) =  mean(aa) ;
    DH_trim(r) =  trimmean(aa,5/100) ;
end

[DH_order,ind_order]=sort(DH) ;
[DH_order_trim,ind_order_trim]=sort(DH_trim) ;

figure(1)
    plot((1:R),DH_order_trim,(1:R),DH_order)
    legend 'trimmed' 'complete'
%    

% SHOW

for ir = 1:R
    r = ind_order_trim(ir) ; 
    fprintf('%3d %3d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n',r, ir, sgp12km_calibr(r,:)) ;
    y = squeeze((sgp12km_daily_rain_June_mean(r,:,:))) ;
    figure(100)
        subplot(121)
            plot((1:R),DH_order_trim,(1:R),DH_order)
            aa = sprintf('%d %g',ir, DH_order_trim(ir) ) ;
            text(ir,DH_order_trim(ir),['(' '$\leftarrow$  ' aa ')'] ,'FontSize',15,'interpreter','latex')
            legend 'trimmed' 'complete'
        subplot(122)
            [ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
            pause(0.5) ;
end


% sgp25km_daily_rain_obs -------------------------------------------------

clearvars -except x1_min x1_max x2_min x2_max xlims

R = 150;

load ./sgp25km_daily_rain_obs_mean.mat sgp25km_daily_rain_obs_mean J K
load ./sgp25km_daily_rain_June_mean.mat sgp25km_daily_rain_June_mean
load ./sgp25km_lat_lon2d.mat sgp25km_lat_lon2d
load ./sgp25km_calibr.mat sgp25km_calibr

y = sgp25km_daily_rain_obs_mean ;
eta = sgp25km_daily_rain_June_mean ;
x = sgp25km_lat_lon2d ;
x1 = squeeze(x(1,:,:)) ;
x2 = squeeze(x(2,:,:)) ;


ind = find( ...
        (isnan(y(:))==0) ...
        .*(x1(:) > x1_min) ...
        .*(x1(:) < x1_max) ...
        .*(x2(:) > x2_min) ...
        .*(x2(:) < x2_max) ...
        ) ;
    
for r = 1:R
    ee = squeeze(eta(r,:,:)) ;
    aa = abs(ee(ind)-y(ind))/mean(y(ind)) ;
    DH(r) =  mean(aa) ;
    DH_trim(r) =  trimmean(aa,5/100) ;
    
end

[DH_order,ind_order]=sort(DH) ;
[DH_order_trim,ind_order_trim]=sort(DH_trim) ;

figure(2)
    plot((1:R),DH_order_trim,(1:R),DH_order)
    legend 'trimmed' 'complete'
%    

% SHOW

for ir = 1:R
    r = ind_order_trim(ir) ; 
    fprintf('%3d %3d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n',r, ir, sgp25km_calibr(r,:)) ;
    y = squeeze((sgp25km_daily_rain_June_mean(r,:,:))) ;
    figure(100)
        subplot(121)
            plot((1:R),DH_order_trim,(1:R),DH_order)
            aa = sprintf('%d %g',ir, DH_order_trim(ir) ) ;
            text(ir,DH_order_trim(ir),['(' '$\leftarrow$  ' aa ')'] ,'FontSize',15,'interpreter','latex')
            legend 'trimmed' 'complete'
        subplot(122)
            [ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
            pause(0.5) ;
end

% sgp50km_daily_rain_obs -------------------------------------------------

clearvars -except x1_min x1_max x2_min x2_max xlims

R = 140;

load ./sgp50km_daily_rain_obs_mean.mat sgp50km_daily_rain_obs_mean J K
load ./sgp50km_daily_rain_June_mean.mat sgp50km_daily_rain_June_mean
load ./sgp50km_lat_lon2d.mat sgp50km_lat_lon2d
load ./sgp50km_calibr.mat sgp50km_calibr

y = sgp50km_daily_rain_obs_mean ;
eta = sgp50km_daily_rain_June_mean ;
x = sgp50km_lat_lon2d ;
x1 = squeeze(x(1,:,:)) ;
x2 = squeeze(x(2,:,:)) ;


ind = find( ...
        (isnan(y(:))==0) ...
        .*(x1(:) > x1_min) ...
        .*(x1(:) < x1_max) ...
        .*(x2(:) > x2_min) ...
        .*(x2(:) < x2_max) ...
        ) ;
    
for r = 1:R
    ee = squeeze(eta(r,:,:)) ;
    aa = abs(ee(ind)-y(ind))/mean(y(ind)) ;
    DH(r) =  mean(aa) ;
    DH_trim(r) =  trimmean(aa,5/100) ;
    
end

[DH_order,ind_order]=sort(DH) ;
[DH_order_trim,ind_order_trim]=sort(DH_trim) ;

figure(3)
    plot((1:R),DH_order_trim,(1:R),DH_order)
    legend 'trimmed' 'complete'
%  

for ir = 1:R
    r = ind_order_trim(ir) ; 
    fprintf('%3d %3d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n',r, ir, sgp50km_calibr(r,:)) ;
    y = squeeze((sgp50km_daily_rain_June_mean(r,:,:))) ;
    figure(100)
        subplot(121)
            plot((1:R),DH_order_trim,(1:R),DH_order)
            aa = sprintf('%d %g',ir, DH_order_trim(ir) ) ;
            text(ir,DH_order_trim(ir),['(' '$\leftarrow$  ' aa ')'] ,'FontSize',15,'interpreter','latex')
            legend 'trimmed' 'complete'
        subplot(122)
            [ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
            pause(0.5) ;
end

% cam25km_daily_rain_obs -------------------------------------------------

clearvars -except x1_min x1_max x2_min x2_max xlims

R = 150;

load ./cam25km_daily_rain_obs_mean.mat cam25km_daily_rain_obs_mean J K
load ./cam25km_daily_rain_June_mean.mat cam25km_daily_rain_June_mean
load ./cam25km_lat_lon2d.mat cam25km_lat_lon2d
load ./sgp12km_calibr.mat sgp25km_calibr

y = cam25km_daily_rain_obs_mean ;
eta = cam25km_daily_rain_June_mean ;
x = cam25km_lat_lon2d ;
x1 = squeeze(x(1,:,:)) ;
x2 = squeeze(x(2,:,:)) ;


ind = find( ...
        (isnan(y(:))==0) ...
        .*(x1(:) > x1_min) ...
        .*(x1(:) < x1_max) ...
        .*(x2(:) > x2_min) ...
        .*(x2(:) < x2_max) ...
        ) ;
    
for r = 1:R
    ee = squeeze(eta(r,:,:)) ;
    aa = abs(ee(ind)-y(ind))/mean(y(ind)) ;
    DH(r) =  mean(aa) ;
    DH_trim(r) =  trimmean(aa,5/100) ;
    
end

[DH_order,ind_order]=sort(DH) ;
[DH_order_trim,ind_order_trim]=sort(DH_trim) ;

figure(4)
    plot((1:R),DH_order_trim,(1:R),DH_order)
    legend 'trimmed' 'complete'
%    
    
for ir = 1:R
    r = ind_order_trim(ir) ; 
    fprintf('%3d %3d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n',r, ir, cam25km_calibr(r,:)) ;
    y = squeeze((cam25km_daily_rain_June_mean(r,:,:))) ;
    figure(100)
        subplot(121)
            plot((1:R),DH_order_trim,(1:R),DH_order)
            aa = sprintf('%d %g',ir, DH_order_trim(ir) ) ;
            text(ir,DH_order_trim(ir),['(' '$\leftarrow$  ' aa ')'] ,'FontSize',15,'interpreter','latex')
            legend 'trimmed' 'complete'
        subplot(122)
            [ax_m,gca_m,gcf_m]=PrintMap(x1,x2,y,xlims) ;
            pause(0.5) ;
end


    