%% Estimate RMS tidal displacement
%  Pang & Nikurashin (2026)
%
%  Purpose:
%  Estimate the RMS vertical displacement associated with selected
%  barotropic tidal constituents interacting with bathymetric slopes.
%
%  Main output:
%  displacement_rms: RMS vertical displacement [m]
%
%  Important:
%  Check all input-grid dimension conventions before running.

clear
addpath('~/bin');

%% 1. Define the calculation domain and physical constants

% Calculation domain: Indonesian Seas
lon_min_deg = 110;
lon_max_deg = 140;
lat_min_deg = -15;
lat_max_deg = 9;

% Earth radius
earth_radius_m = 6.371e6;                 % [m]

% Minimum water depth used for tidal velocities
minimum_depth_m = -50;                    % [m]

% Number of tidal constituents
n_constituents = 8;

% Tidal constituent periods [s]
% Order: M2, S2, N2, K2, K1, O1, P1, Q1
tidal_periods = [44714.16, 43200.00, 45569.88, 43081.92,...
    86164.20, 92949.48, 86637.24, 96726.24];

constit_names = {'M2','S2','N2','K2','K1','O1','P1','Q1'};

%% 2. Read bathymetry

topography_file = './STRM30PLUS.nc';

lon_bathy_deg = ncread(topography_file, 'x');
lat_bathy_deg = ncread(topography_file, 'y');
topography_m  = ncread(topography_file, 'z');

% Ocean depth: negative downward; land set to zero.
% Positive topographic elevations represent land.
depth_m = topography_m;
depth_m(depth_m > 0) = 0;

ii = find(lon_bathy_deg > lon_min_deg & ...
          lon_bathy_deg <= lon_max_deg);

jj = find(lat_bathy_deg > lat_min_deg & ...
          lat_bathy_deg <= lat_max_deg);

lon_bathy_deg = lon_bathy_deg(ii);
lat_bathy_deg = lat_bathy_deg(jj);

depth_m = depth_m(ii,jj);

%% 3. Read TPXO tidal velocities

tpxo_file = './tpxo9.v4.nc';

lon_u_deg = ncread(tpxo_file, 'lon_u');
lon_v_deg = ncread(tpxo_file, 'lon_v');
lat_u_deg = ncread(tpxo_file, 'lat_u');
lat_v_deg = ncread(tpxo_file, 'lat_v');

% Tidal velocity [m s^-1]
u_tide_raw = ncread(tpxo_file, 'ua');
v_tide_raw = ncread(tpxo_file, 'va');

% Select the calculation domain
% Raw data dimensions: [lat_deg, lon_deg, n_constituents] 
ii = find(lon_u_deg(1,:) > lon_min_deg & ...
          lon_u_deg(1,:) <= lon_max_deg);

jj = find(lat_u_deg(:,1) > lat_min_deg & ...
          lat_u_deg(:,1) <= lat_max_deg);

lon_u_deg = lon_u_deg(jj,ii);
lat_u_deg = lat_u_deg(jj,ii);
lon_v_deg = lon_v_deg(jj,ii);
lat_v_deg = lat_v_deg(jj,ii);

u_tide_raw = u_tide_raw(jj,ii,:);
v_tide_raw = v_tide_raw(jj,ii,:);

u_tide_raw(u_tide_raw == 0) = NaN;
v_tide_raw(v_tide_raw == 0) = NaN;

%% 4. Interpolate tidal velocities onto the bathymetry grid

% Interpolated velocity arrays:
% Dimensions: [n_lon_bathy, n_lat_bathy, n_constituents]

[nx,ny] = size(depth_m);
ua = zeros(nx,ny,n_constituents);
va = zeros(nx,ny,n_constituents);

for n = 1:8
ua(:,:,n) = interp2(lat_u_deg(:,1),lon_u_deg(1,:),...
    u_tide_raw(:,:,n)',lat_bathy_deg,lon_bathy_deg');
va(:,:,n) = interp2(lat_v_deg(:,1),lon_v_deg(1,:),...
    v_tide_raw(:,:,n)',lat_bathy_deg,lon_bathy_deg');
end

[yc,xc] = meshgrid(lon_bathy_deg',lat_bathy_deg);

ws = [0.1 0.1 10];
hc = depth_m; hc(hc>0)=0; 

for n = 1:8
disp(n)
ua(:,:,n) = my_interpolation(lon_u_deg',lat_u_deg',0*lon_u_deg',...
    u_tide_raw(:,:,n)',xc,yc,0*xc,ua(:,:,n),hc,ws);
va(:,:,n) = my_interpolation(lon_v_deg',lat_v_deg',0*lon_v_deg',...
    v_tide_raw(:,:,n)',xc,yc,0*xc,va(:,:,n),hc,ws);
end

ua(depth_m>minimum_depth_m) = 0;
va(depth_m>minimum_depth_m) = 0;
% Remove the shallow regions

%% 5. Calculate bathymetric gradients

% Geographic coordinates [rad]
lon_grid_rad = deg2rad(xc);
lat_grid_rad = deg2rad(yc);

% Grid-cell distances [m]
% These expressions follow the original coordinate orientation.
delta_lon_rad = diff(lon_grid_rad, 1, 1);
delta_lat_rad = diff(lat_grid_rad, 1, 2);

dx_m = earth_radius_m .* cos(lat_grid_rad(1:end-1,:)) ...
       .* delta_lon_rad;

dy_m = earth_radius_m .* delta_lat_rad;

% Bathymetric gradients [m/m], dimensionless.
dz_dx = diff(depth_m, 1, 1) ./ dx_m;
dz_dy = diff(depth_m, 1, 2) ./ dy_m;

% Restore the original array dimensions using edge replication.
dz_dx = cat(1, dz_dx, dz_dx(end,:));
dz_dy = cat(2, dz_dy, dz_dy(:,end));

%% 6. Define the vertical grid and weighting

% layer_depth here be defined by the model vertical grid [m].
z_layer_m = layer_depth;

if any(z_layer_m(:) > 0)
    error('layer_depth must use negative depths below the surface.');
end

nz = numel(z_layer_m);
% nz is 100 for our regional model

% Replicate bottom depth and vertical coordinates.
depth_m = repmat(depth_m, 1, 1, nz);
z_layer_m  = repmat(z_layer_m, nx, ny, 1);

% Vertical weighting factor from the original formulation.
% This factor is zero below the bottom and decreases with depth.
dummy = abs(z_layer_m) ./ depth_m;
dummy(dummy > 1) = 0;

%% 7. Calculate the constituent displacement contributions

% Accumulated variance [m^2]
displace_var = zeros(nx, ny, nz);

for k = 1:n_constituents

    % Squared bottom-velocity: [m^2 s^-2]
    wbot2 = ua(:,:,k).^2 .* dz_dx.^2 + ...
            va(:,:,k).^2 .* dz_dy.^2;

    % Replicate over the vertical dimension.
    wbot2 = repmat(wbot2, 1, 1, nz);

    % Contribution from this constituent [m^2].
    constit_var = 0.5 .* wbot2 .* dummy.^2 ...
        .*(tidal_periods(k) ./ (2*pi)).^2;

    displace_var = displace_var + constit_var;
end

%% 8. Calculate RMS displacement

displacement_rms = sqrt(displace_var);
