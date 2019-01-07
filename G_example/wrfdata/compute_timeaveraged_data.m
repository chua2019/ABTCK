
% COmpute time averaged data set

% ========================================================================

clear all
close all

% sgp12km ================================================================

% sgp12km_daily_rain_obs -------------------------------------------------

clear all ; close all ;

% y_exp

I = 30;
J = 182 ;
K = 187 ;

data = load('../ansi/sgp12km_daily_rain_obs.dat') ;
% y_exp = zeros(I,J,K) ;
% ll = 0 ;
% for i = 1:I
%     for j = 1:J
%         for k = 1:K
%            ll = ll + 1 ;
%            y_exp(i,j,k) = data(ll) ; 
%         end
%     end
% end
y_exp = permute(reshape(data,K,J,I),[3 2 1]) ;
y_exp = squeeze(mean(y_exp,1)) ;
y_exp(find(y_exp<0)) = nan ;

sgp12km_daily_rain_obs_mean = y_exp ;

save sgp12km_daily_rain_obs_mean.mat sgp12km_daily_rain_obs_mean J K

clear all

% sgp12km_lat_lon2d ------------------------------------------------------

clear all ; close all ;

J = 182 ;
K = 187 ;
data = load('../ansi/sgp12km_lat_lon2d.dat') ;
% x_exp = zeros(2,J,K) ;
% ll = 0 ;
% for ii = 1:2
%     for j = 1:J
%         for k = 1:K
%            ll = ll + 1 ;
%            x_exp(ii,j,k) = data(ll) ; 
%         end
%     end
% end
x_exp = permute(reshape(data,K,J,2),[3 2 1]) ;
sgp12km_lat_lon2d = x_exp ;

save sgp12km_lat_lon2d.mat sgp12km_lat_lon2d J K

clear all

% sgp12km_daily_rain_June ------------------------------------------------

clear all ; close all ;

load sgp12km_daily_rain_obs_mean.mat sgp12km_daily_rain_obs_mean
MYIND = find(isnan(sgp12km_daily_rain_obs_mean)) ;
clearvars -except MYIND ;

R = 150;
I = 30;
J = 182 ;
K = 187 ;

data = load('../ansi/sgp12km_daily_rain_June.dat') ;
% data_mat = zeros(R,I,J,K) ;
% ll = 0 ;
% for r = 1:R
%     for i = 1:I
%         for j = 1:J
%             for k = 1:K
%                ll = ll + 1 ;
%                data_mat(r,i,j,k) = data(ll) ; 
%             end
%         end
%     end
% end
data_mat = permute(reshape(data,K,J,I,R),[4 3 2 1]) ;
data_mat_mean = squeeze(mean(data_mat,2)) ;
for r = 1:R
    a = squeeze(data_mat_mean(r,:,:)) ;
    a(MYIND) = nan ;
    data_mat_mean(r,:,:) = a ;
end
sgp12km_daily_rain_June_mean = data_mat_mean ;

save sgp12km_daily_rain_June_mean.mat sgp12km_daily_rain_June_mean R J K

clear all

% sgp12km_calibration ------------------------------------------------

clear all ; close all ;

sgp12km_calibr = get_calibpar('sgp12km') ;

R = size(sgp12km_calibr, 1) ;

save sgp12km_calibr.mat sgp12km_calibr R

clear all

% sgp25km ================================================================

% sgp25km_daily_rain_obs -------------------------------------------------

clear all ; close all ;

% y_exp

I = 30;
J = 87 ;
K = 89 ;

data = load('../ansi/sgp25km_daily_rain_obs.dat') ;
% y_exp = zeros(I,J,K) ;
% ll = 0 ;
% for i = 1:I
%     for j = 1:J
%         for k = 1:K
%            ll = ll + 1 ;
%            y_exp(i,j,k) = data(ll) ; 
%         end
%     end
% end
y_exp = permute(reshape(data,K,J,I),[3 2 1]) ;
y_exp = squeeze(mean(y_exp,1)) ;
y_exp(find(y_exp<0)) = nan ;

sgp25km_daily_rain_obs_mean = y_exp ;

save sgp25km_daily_rain_obs_mean.mat sgp25km_daily_rain_obs_mean J K

clear all

% sgp25km_lat_lon2d ------------------------------------------------------

clear all ; close all ;

J = 87 ;
K = 89 ;

data = load('../ansi/sgp25km_lat_lon2d.dat') ;
% x_exp = zeros(2,J,K) ;
% ll = 0 ;
% for ii = 1:2
%     for j = 1:J
%         for k = 1:K
%            ll = ll + 1 ;
%            x_exp(ii,j,k) = data(ll) ; 
%         end
%     end
% end
x_exp = permute(reshape(data,K,J,2),[3 2 1]) ;
sgp25km_lat_lon2d = x_exp ;

save sgp25km_lat_lon2d.mat sgp25km_lat_lon2d J K

clear all

% sgp25km_daily_rain_June ------------------------------------------------

clear all ; close all ;

load sgp25km_daily_rain_obs_mean.mat sgp25km_daily_rain_obs_mean
MYIND = find(isnan(sgp25km_daily_rain_obs_mean)) ;
clearvars -except MYIND ;

R = 150;
I = 30;
J = 87 ;
K = 89 ;

data = load('../ansi/sgp25km_daily_rain_June.dat') ;
% data_mat = zeros(R,I,J,K) ;
% ll = 0 ;
% for r = 1:R
%     for i = 1:I
%         for j = 1:J
%             for k = 1:K
%                ll = ll + 1 ;
%                data_mat(r,i,j,k) = data(ll) ; 
%             end
%         end
%     end
% end
data_mat = permute(reshape(data,K,J,I,R),[4 3 2 1]) ;
data_mat_mean = squeeze(mean(data_mat,2)) ;
for r = 1:R
    a = squeeze(data_mat_mean(r,:,:)) ;
    a(MYIND) = nan ;
    data_mat_mean(r,:,:) = a ;
end
sgp25km_daily_rain_June_mean = data_mat_mean ;

save sgp25km_daily_rain_June_mean.mat sgp25km_daily_rain_June_mean R J K

clear all

% sgp25km_calibration ------------------------------------------------

clear all ; close all ;

sgp25km_calibr = get_calibpar('sgp25km') ;

R = size(sgp25km_calibr, 1) ;

save sgp25km_calibr.mat sgp25km_calibr R

clear all

% sgp50km ================================================================

% sgp50km_daily_rain_obs -------------------------------------------------

clear all ; close all ;

% y_exp

I = 30;
J = 43 ;
K = 44 ;

data = load('../ansi/sgp50km_daily_rain_obs.dat') ;
% y_exp = zeros(I,J,K) ;
% ll = 0 ;
% for i = 1:I
%     for j = 1:J
%         for k = 1:K
%            ll = ll + 1 ;
%            y_exp(i,j,k) = data(ll) ; 
%         end
%     end
% end
y_exp = permute(reshape(data,K,J,I),[3 2 1]) ;
y_exp = squeeze(mean(y_exp,1)) ;
y_exp(find(y_exp<0)) = nan ;

sgp50km_daily_rain_obs_mean = y_exp ;

save sgp50km_daily_rain_obs_mean.mat sgp50km_daily_rain_obs_mean J K

clear all

% sgp50km_lat_lon2d ------------------------------------------------------

clear all ; close all ;

J = 43 ;
K = 44 ;
data = load('../ansi/sgp50km_lat_lon2d.dat') ;
% x_exp = zeros(2,J,K) ;
% ll = 0 ;
% for ii = 1:2
%     for j = 1:J
%         for k = 1:K
%            ll = ll + 1 ;
%            x_exp(ii,j,k) = data(ll) ; 
%         end
%     end
% end
x_exp = permute(reshape(data,K,J,2),[3 2 1]) ;
sgp50km_lat_lon2d = x_exp ;

save sgp50km_lat_lon2d.mat sgp50km_lat_lon2d J K

clear all

% sgp50km_daily_rain_June ------------------------------------------------

clear all ; close all ;

load sgp50km_daily_rain_obs_mean.mat sgp50km_daily_rain_obs_mean
MYIND = find(isnan(sgp50km_daily_rain_obs_mean)) ;
clearvars -except MYIND ;

R = 141; % ACTUALLY 140 = 141-1
I = 30;
J = 43 ;
K = 44 ;

data = load('../ansi/sgp50km_daily_rain_June.dat') ;
% data_mat = zeros(R,I,J,K) ;
% ll = 0 ;
% for r = 1:R
%     for i = 1:I
%         for j = 1:J
%             for k = 1:K
%                ll = ll + 1 ;
%                data_mat(r,i,j,k) = data(ll) ; 
%             end
%         end
%     end
% end
data_mat = permute(reshape(data,K,J,I,R),[4 3 2 1]) ;
data_mat_mean = squeeze(mean(data_mat,2)) ;
for r = 1:R
    a = squeeze(data_mat_mean(r,:,:)) ;
    a(MYIND) = nan ;
    data_mat_mean(r,:,:) = a ;
end
sgp50km_daily_rain_June_mean = data_mat_mean(1:end-1,:,:) ;
R = size(sgp50km_daily_rain_June_mean,1) ;

save sgp50km_daily_rain_June_mean.mat sgp50km_daily_rain_June_mean R J K

clear all

% sgp50km_calibration ------------------------------------------------

clear all ; close all ;

sgp50km_calibr = get_calibpar('sgp50km') ;

sgp50km_calibr = sgp50km_calibr(11:end-1,:) ;

R = size(sgp50km_calibr, 1) ;

save sgp50km_calibr.mat sgp50km_calibr R

clear all

% cam25km ================================================================

% cam25km_daily_rain_obs -------------------------------------------------

clear all ; close all ;

% y_exp

I = 30;
J = 87 ;
K = 89 ;

data = load('../ansi/cam25km_daily_rain_obs.dat') ;
% y_exp = zeros(I,J,K) ;
% ll = 0 ;
% for i = 1:I
%     for j = 1:J
%         for k = 1:K
%            ll = ll + 1 ;
%            y_exp(i,j,k) = data(ll) ; 
%         end
%     end
% end
y_exp = permute(reshape(data,K,J,I),[3 2 1]) ;
y_exp = squeeze(mean(y_exp,1)) ;
y_exp(find(y_exp<0)) = nan ;

cam25km_daily_rain_obs_mean = y_exp ;

save cam25km_daily_rain_obs_mean.mat cam25km_daily_rain_obs_mean J K

clear all

% cam25km_lat_lon2d ------------------------------------------------------

clear all ; close all ;

J = 87 ;
K = 89 ;

data = load('../ansi/cam25km_lat_lon2d.dat') ;
% x_exp = zeros(2,J,K) ;
% ll = 0 ;
% for ii = 1:2
%     for j = 1:J
%         for k = 1:K
%            ll = ll + 1 ;
%            x_exp(ii,j,k) = data(ll) ; 
%         end
%     end
% end
x_exp = permute(reshape(data,K,J,2),[3 2 1]) ;
cam25km_lat_lon2d = x_exp ;

save cam25km_lat_lon2d.mat cam25km_lat_lon2d J K

clear all

% cam25km_daily_rain_June ------------------------------------------------

clear all ; close all ;

load cam25km_daily_rain_obs_mean.mat cam25km_daily_rain_obs_mean
MYIND = find(isnan(cam25km_daily_rain_obs_mean)) ;
clearvars -except MYIND ;

R = 150;
I = 30;
J = 87 ;
K = 89 ;

data = load('../ansi/cam25km_daily_rain_June.dat') ;
% data_mat = zeros(R,I,J,K) ;
% ll = 0 ;
% for r = 1:R
%     for i = 1:I
%         for j = 1:J
%             for k = 1:K
%                ll = ll + 1 ;
%                data_mat(r,i,j,k) = data(ll) ; 
%             end
%         end
%     end
% end
data_mat = permute(reshape(data,K,J,I,R),[4 3 2 1]) ;
data_mat_mean = squeeze(mean(data_mat,2)) ;
for r = 1:R
    a = squeeze(data_mat_mean(r,:,:)) ;
    a(MYIND) = nan ;
    data_mat_mean(r,:,:) = a ;
end
cam25km_daily_rain_June_mean = data_mat_mean ;

save cam25km_daily_rain_June_mean.mat cam25km_daily_rain_June_mean R J K

clear all

% cam25km_calibration ------------------------------------------------

clear all ; close all ;

cam25km_calibr = get_calibpar('cam25km') ;

R = size(cam25km_calibr, 1) ;

save cam25km_calibr.mat cam25km_calibr R

clear all



% ========================================================================
% ========================================================================

% truncate the samples ...

% ========================================================================
% ========================================================================

clear all
close all

x1_min = 25.5 ;
x1_max = 41 ;
x2_min = -109.3 ;
x2_max = -93.31 ;
    
xlims = [ [25.50 41] ; [-109.3 -93.31] ] ;

% sgp12km_daily_rain -------------------------------------------------

clearvars -except x1_min x1_max x2_min x2_max xlims

load ./sgp12km_daily_rain_obs_mean.mat sgp12km_daily_rain_obs_mean J K 
load ./sgp12km_daily_rain_June_mean.mat sgp12km_daily_rain_June_mean
load ./sgp12km_lat_lon2d.mat sgp12km_lat_lon2d
load ./sgp12km_calibr.mat sgp12km_calibr R

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

R_trunc = 77 ;
R = R_trunc ;

sgp12km_daily_rain_June_mean = ...
    sgp12km_daily_rain_June_mean(ind_order_trim(1:R),:,:) ;

save ./sgp12km_daily_rain_June_mean_trunc.mat ...
            sgp12km_daily_rain_June_mean J K R ...
            ind_order

% sgp25km_daily_rain -------------------------------------------------

clearvars -except x1_min x1_max x2_min x2_max xlims

load ./sgp25km_daily_rain_obs_mean.mat sgp25km_daily_rain_obs_mean J K
load ./sgp25km_daily_rain_June_mean.mat sgp25km_daily_rain_June_mean
load ./sgp25km_lat_lon2d.mat sgp25km_lat_lon2d
load ./sgp25km_calibr.mat sgp25km_calibr R

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

R_trunc = 85 ;
R = R_trunc ;

sgp25km_daily_rain_June_mean = ...
    sgp25km_daily_rain_June_mean(ind_order_trim(1:R),:,:) ;

save ./sgp25km_daily_rain_June_mean_trunc.mat ...
            sgp25km_daily_rain_June_mean J K R ...
            ind_order

% sgp50km_daily_rain -------------------------------------------------

clearvars -except x1_min x1_max x2_min x2_max xlims

load ./sgp50km_daily_rain_obs_mean.mat sgp50km_daily_rain_obs_mean J K
load ./sgp50km_daily_rain_June_mean.mat sgp50km_daily_rain_June_mean
load ./sgp50km_lat_lon2d.mat sgp50km_lat_lon2d
load ./sgp50km_calibr.mat sgp50km_calibr R

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

R_trunc = 85 ;
R = R_trunc ;

sgp50km_daily_rain_June_mean = ...
    sgp50km_daily_rain_June_mean(ind_order_trim(1:R),:,:) ;

save ./sgp50km_daily_rain_June_mean_trunc.mat ...
            sgp50km_daily_rain_June_mean J K R ...
            ind_order

% cam25km_daily_rain -------------------------------------------------

clearvars -except x1_min x1_max x2_min x2_max xlims

load ./cam25km_daily_rain_obs_mean.mat cam25km_daily_rain_obs_mean J K
load ./cam25km_daily_rain_June_mean.mat cam25km_daily_rain_June_mean
load ./cam25km_lat_lon2d.mat cam25km_lat_lon2d
load ./cam25km_calibr.mat cam25km_calibr R

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

R_trunc = 73 ;
R = R_trunc ;

cam25km_daily_rain_June_mean = ...
    cam25km_daily_rain_June_mean(ind_order_trim(1:R),:,:) ;

save ./cam25km_daily_rain_June_mean_trunc.mat ...
            cam25km_daily_rain_June_mean J K R ...
            ind_order


% sgp12km_calibration trunc ----------------------------------------------

clear all ; close all ;

load ./sgp12km_daily_rain_June_mean_trunc.mat R ind_order

sgp12km_calibr = get_calibpar('sgp12km') ;

sgp12km_calibr = sgp12km_calibr(ind_order(1:R),:) ;

save sgp12km_calibr_trunc.mat sgp12km_calibr R 

clear all

% sgp25km_calibration trunc ----------------------------------------------

clear all ; close all ;

load ./sgp25km_daily_rain_June_mean_trunc.mat R ind_order

sgp25km_calibr = get_calibpar('sgp25km') ;
sgp25km_calibr = sgp25km_calibr(ind_order(1:R),:) ;

save sgp25km_calibr_trunc.mat sgp25km_calibr R 

clear all

% sgp50km_calibration trunc ----------------------------------------------

clear all ; close all ;

load ./sgp50km_daily_rain_June_mean_trunc.mat R ind_order

sgp50km_calibr = get_calibpar('sgp50km') ;
sgp50km_calibr = sgp50km_calibr(ind_order(1:R),:) ;

save sgp50km_calibr_trunc.mat sgp50km_calibr R 

clear all

% cam25km_calibration trunc ----------------------------------------------

clear all ; close all ;

load ./cam25km_daily_rain_June_mean_trunc.mat R ind_order

cam25km_calibr = get_calibpar('cam25km') ;
cam25km_calibr = cam25km_calibr(ind_order(1:R),:) ;

save cam25km_calibr_trunc.mat cam25km_calibr R 

clear all



