
% Generate training data set

% ========================================================================

clear all
close all

x1_min = 25.5 ;
x1_max = 41 ;
x2_min = -109.3 ;
x2_max = -93.31 ;
    
xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;
data_out = './data_WRF' ;

% experimental ===========================================================

clearvars -except x1_min x1_max x2_min x2_max data_out xlims

n_sys = 100 ;

% y_exp

load ./sgp12km_daily_rain_obs_mean.mat sgp12km_daily_rain_obs_mean J K
load ./sgp12km_lat_lon2d.mat sgp12km_lat_lon2d

y = sgp12km_daily_rain_obs_mean ;
x1 = squeeze(sgp12km_lat_lon2d(1,:,:)) ;
x2 = squeeze(sgp12km_lat_lon2d(2,:,:)) ;

% find 
ind = find( ... 
        (isnan(y(:))==0) ...
        .*(x1(:) > x1_min) ...
        .*(x1(:) < x1_max) ...
        .*(x2(:) > x2_min) ...
        .*(x2(:) < x2_max) ...
        ) ;
    
n_sys = min(n_sys,length(ind)) ;
ind = ind(randperm(length(ind),n_sys)) ;
    
y = y(ind) ;
x1 = x1(ind) ;
x2 = x2(ind) ;

% save

fileID = fopen(sprintf('%s/y_exp.dat',data_out),'w') ;
for i=1:n_sys
    fprintf(fileID,'%f \n', y(i) );
end
fprintf(fileID,'\n');
fclose(fileID);

fileID = fopen(sprintf('%s/x_exp.dat',data_out),'w') ;
for i=1:n_sys
    fprintf(fileID,'%f %f \n', x1(i), x2(i) );
end
fprintf(fileID,'\n');
fclose(fileID);

% sim 1 (param_sgp12km.dat) ==============================================

clearvars -except x1_min x1_max x2_min x2_max data_out xlims

load ./sgp12km_daily_rain_June_mean_trunc.mat sgp12km_daily_rain_June_mean R J K
load ./sgp12km_lat_lon2d.mat sgp12km_lat_lon2d
load ./sgp12km_calibr_trunc.mat sgp12km_calibr R 

n_sim = 300 ;
nx_sim = ceil(max(1,n_sim/R)) ;

eta_raw = sgp12km_daily_rain_June_mean ;
x1_raw = squeeze(sgp12km_lat_lon2d(1,:,:)) ;
x2_raw = squeeze(sgp12km_lat_lon2d(2,:,:)) ;
t_raw = sgp12km_calibr ;

eta = [] ;
x = [] ;
t = [] ;

for r = 1:R
    
   eta_r =  squeeze(eta_raw(r,:,:)) ;
   
   ind = find( ... 
        (isnan(eta_r(:))==0) ...
        .*(x1_raw(:) > x1_min) ...
        .*(x1_raw(:) < x1_max) ...
        .*(x2_raw(:) > x2_min) ...
        .*(x2_raw(:) < x2_max) ...
        ) ;
    n_sim = min(n_sim,length(ind)) ;
    ind = ind(randperm(length(ind),nx_sim)) ;
   
    eta = [eta ; eta_r(ind)] ;
    x = [x ; [x1_raw(ind) x2_raw(ind)]] ; 
    t = [t ; repmat(t_raw(r,:),nx_sim,1)] ;
    
end

n_sim = length(eta) ;

% save

fileID = fopen(sprintf('%s/eta_sim1.dat',data_out),'w') ;
for i=1:n_sim
    fprintf(fileID,'%f \n', eta(i) );
end
fprintf(fileID,'\n');
fclose(fileID);

fileID = fopen(sprintf('%s/x_sim1.dat',data_out),'w') ;
for i=1:n_sim
    fprintf(fileID,'%f %f \n',x(i,1), x(i,2) );
end
fprintf(fileID,'\n');
fclose(fileID);

fileID = fopen(sprintf('%s/t_sim1.dat',data_out),'w') ;
for i=1:n_sim
    fprintf(fileID,'%f %f %f %f %f %f\n', ...
        t(i,1), t(i,2), t(i,3), t(i,4), t(i,5), t(i,6) );
end
fprintf(fileID,'\n');
fclose(fileID);

% sim 2 (param_sgp25km.dat) ==============================================

clearvars -except x1_min x1_max x2_min x2_max data_out xlims

load ./sgp25km_daily_rain_June_mean_trunc.mat sgp25km_daily_rain_June_mean R J K
load ./sgp25km_lat_lon2d.mat sgp25km_lat_lon2d
load ./sgp25km_calibr_trunc.mat sgp25km_calibr R 

n_sim = 300 ;
nx_sim = ceil(max(1,n_sim/R)) ;

eta_raw = sgp25km_daily_rain_June_mean ;
x1_raw = squeeze(sgp25km_lat_lon2d(1,:,:)) ;
x2_raw = squeeze(sgp25km_lat_lon2d(2,:,:)) ;
t_raw = sgp25km_calibr ;

eta = [] ;
x = [] ;
t = [] ;

for r = 1:R
    
   eta_r =  squeeze(eta_raw(r,:,:)) ;
   
   ind = find( ... 
        (isnan(eta_r(:))==0) ...
        .*(x1_raw(:) > x1_min) ...
        .*(x1_raw(:) < x1_max) ...
        .*(x2_raw(:) > x2_min) ...
        .*(x2_raw(:) < x2_max) ...
        ) ;
    n_sim = min(n_sim,length(ind)) ;
    ind = ind(randperm(length(ind),nx_sim)) ;
   
    eta = [eta ; eta_r(ind)] ;
    x = [x ; [x1_raw(ind) x2_raw(ind)]] ; 
    t = [t ; repmat(t_raw(r,:),nx_sim,1)] ;
    
end

n_sim = length(eta) ;

% save

fileID = fopen(sprintf('%s/eta_sim2.dat',data_out),'w') ;
for i=1:n_sim
    fprintf(fileID,'%f \n', eta(i) );
end
fprintf(fileID,'\n');
fclose(fileID);

fileID = fopen(sprintf('%s/x_sim2.dat',data_out),'w') ;
for i=1:n_sim
    fprintf(fileID,'%f %f \n',x(i,1), x(i,2) );
end
fprintf(fileID,'\n');
fclose(fileID);

fileID = fopen(sprintf('%s/t_sim2.dat',data_out),'w') ;
for i=1:n_sim
    fprintf(fileID,'%f %f %f %f %f %f\n', ...
        t(i,1), t(i,2), t(i,3), t(i,4), t(i,5), t(i,6) );
end
fprintf(fileID,'\n');
fclose(fileID);

% sim 3 (param_sgp50km.dat) ==============================================

clearvars -except x1_min x1_max x2_min x2_max data_out xlims

load ./sgp50km_daily_rain_June_mean_trunc.mat sgp50km_daily_rain_June_mean R J K
load ./sgp50km_lat_lon2d.mat sgp50km_lat_lon2d
load ./sgp50km_calibr_trunc.mat sgp50km_calibr R 

n_sim = 300 ;
nx_sim = ceil(max(1,n_sim/R)) ;

eta_raw = sgp50km_daily_rain_June_mean ;
x1_raw = squeeze(sgp50km_lat_lon2d(1,:,:)) ;
x2_raw = squeeze(sgp50km_lat_lon2d(2,:,:)) ;
t_raw = sgp50km_calibr ;

eta = [] ;
x = [] ;
t = [] ;

for r = 1:R
    
   eta_r =  squeeze(eta_raw(r,:,:)) ;
   
   ind = find( ... 
        (isnan(eta_r(:))==0) ...
        .*(x1_raw(:) > x1_min) ...
        .*(x1_raw(:) < x1_max) ...
        .*(x2_raw(:) > x2_min) ...
        .*(x2_raw(:) < x2_max) ...
        ) ;
    n_sim = min(n_sim,length(ind)) ;
    ind = ind(randperm(length(ind),nx_sim)) ;
   
    eta = [eta ; eta_r(ind)] ;
    x = [x ; [x1_raw(ind) x2_raw(ind)]] ; 
    t = [t ; repmat(t_raw(r,:),nx_sim,1)] ;
    
end

n_sim = length(eta) ;

% save

fileID = fopen(sprintf('%s/eta_sim3.dat',data_out),'w') ;
for i=1:n_sim
    fprintf(fileID,'%f \n', eta(i) );
end
fprintf(fileID,'\n');
fclose(fileID);

fileID = fopen(sprintf('%s/x_sim3.dat',data_out),'w') ;
for i=1:n_sim
    fprintf(fileID,'%f %f \n',x(i,1), x(i,2) );
end
fprintf(fileID,'\n');
fclose(fileID);

fileID = fopen(sprintf('%s/t_sim3.dat',data_out),'w') ;
for i=1:n_sim
    fprintf(fileID,'%f %f %f %f %f %f\n', ...
        t(i,1), t(i,2), t(i,3), t(i,4), t(i,5), t(i,6) );
end
fprintf(fileID,'\n');
fclose(fileID);

% sim 4 (param_cam25km.dat) ==============================================

clearvars -except x1_min x1_max x2_min x2_max data_out xlims

load ./cam25km_daily_rain_June_mean_trunc.mat cam25km_daily_rain_June_mean R J K
load ./cam25km_lat_lon2d.mat cam25km_lat_lon2d
load ./cam25km_calibr_trunc.mat cam25km_calibr R 

n_sim = 300 ;
nx_sim = ceil(max(1,n_sim/R)) ;

eta_raw = cam25km_daily_rain_June_mean ;
x1_raw = squeeze(cam25km_lat_lon2d(1,:,:)) ;
x2_raw = squeeze(cam25km_lat_lon2d(2,:,:)) ;
t_raw = cam25km_calibr ;

eta = [] ;
x = [] ;
t = [] ;

for r = 1:R
    
   eta_r =  squeeze(eta_raw(r,:,:)) ;
   
   ind = find( ... 
        (isnan(eta_r(:))==0) ...
        .*(x1_raw(:) > x1_min) ...
        .*(x1_raw(:) < x1_max) ...
        .*(x2_raw(:) > x2_min) ...
        .*(x2_raw(:) < x2_max) ...
        ) ;
    n_sim = min(n_sim,length(ind)) ;
    ind = ind(randperm(length(ind),nx_sim)) ;
   
    eta = [eta ; eta_r(ind)] ;
    x = [x ; [x1_raw(ind) x2_raw(ind)]] ; 
    t = [t ; repmat(t_raw(r,:),nx_sim,1)] ;
    
end

n_sim = length(eta) ;

% save

fileID = fopen(sprintf('%s/eta_sim4.dat',data_out),'w') ;
for i=1:n_sim
    fprintf(fileID,'%f \n', eta(i) );
end
fprintf(fileID,'\n');
fclose(fileID);

fileID = fopen(sprintf('%s/x_sim4.dat',data_out),'w') ;
for i=1:n_sim
    fprintf(fileID,'%f %f \n',x(i,1), x(i,2) );
end
fprintf(fileID,'\n');
fclose(fileID);

fileID = fopen(sprintf('%s/t_sim4.dat',data_out),'w') ;
for i=1:n_sim
    fprintf(fileID,'%f %f %f %f %f %f\n', ...
        t(i,1), t(i,2), t(i,3), t(i,4), t(i,5), t(i,6) );
end
fprintf(fileID,'\n');
fclose(fileID);









