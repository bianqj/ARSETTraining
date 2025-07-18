function tilemap3(varargin)
% TILEMAP3 - Maps between geographic coordinates and MODIS grids
% converted from https://landweb.modaps.eosdis.nasa.gov/tilecalc
% Usage: tilemap3(projection, pixel_size, direction, type, point_data)
%
% Parameters:
%   projection: 'is_k', 'is_h', 'is_q', 'gh', 'np', 'sp', 'sn', 'ha'
%   pixel_size: 'k', 'h', 'q', 's', 'm', 'l', 'x'
%   direction: 'fwd' (forward) or 'inv' (inverse)
%   type: 'tp' (tile/pixel), 'gm' (global map), 'gp' (global pixel)
%   point_data: varies based on direction and type
%
% Examples:
%   tilemap3('is_q', 'k', 'fwd', 'tp', 38.53, -77.0)
%   tilemap3('is_q', 'k', 'inv', 'tp', 5, 11, 175.90, 1171.23)

if nargin < 5
    print_usage();
    return;
end

% Parse input arguments
projection = varargin{1};
pixel_size = varargin{2};
direction = varargin{3};
type = varargin{4};

% Initialize projection parameters
proj_info = get_projection_info(projection);
pixel_ratio = get_pixel_ratio(pixel_size);

% Initialize tile structure
tile = tile_init(proj_info, pixel_ratio);

% Process based on direction
if strcmp(direction, 'fwd')
    if nargin ~= 6
        error('Forward mapping requires latitude and longitude');
    end
    lat = varargin{5};
    lon = varargin{6};
    
    % Convert to radians
    lat_rad = lat * pi / 180;
    lon_rad = lon * pi / 180;
    
    % Forward mapping
    [itile_line, itile_samp, line, samp] = tile_fwd(tile, lon_rad, lat_rad);
    
    fprintf('lat %.6f  long %.6f  =>', lat, lon);
    
    switch type
        case 'tp'
            fprintf('  vert tile %d  horiz tile %d  line %.2f  samp %.2f\n', ...
                    itile_line, itile_samp, line, samp);
        case 'gm'
            [x, y] = tile_fwd_map(tile, itile_line, itile_samp, line, samp);
            fprintf('  x %.2f  y %.2f\n', x, y);
        case 'gp'
            [line_global, samp_global] = tile_fwd_pix(tile, itile_line, itile_samp, line, samp);
            fprintf('  line %.2f  samp %.2f\n', line_global, samp_global);
    end
    
elseif strcmp(direction, 'inv')
    switch type
        case 'tp'
            if nargin ~= 8
                error('Inverse tile mapping requires 4 parameters');
            end
            itile_line = varargin{5};
            itile_samp = varargin{6};
            line = varargin{7};
            samp = varargin{8};
            fprintf('vert tile %d  horiz tile %d  line %.2f  samp %.2f  =>', ...
                    itile_line, itile_samp, line, samp);
            
        case 'gm'
            if nargin ~= 6
                error('Inverse global map mapping requires x and y coordinates');
            end
            x = varargin{5};
            y = varargin{6};
            fprintf('x %.2f  y %.2f  =>', x, y);
            [itile_line, itile_samp, line, samp] = tile_inv_map(tile, x, y);
            
        case 'gp'
            if nargin ~= 6
                error('Inverse global pixel mapping requires line and sample');
            end
            line_global = varargin{5};
            samp_global = varargin{6};
            fprintf('line %.2f  samp %.2f  =>', line_global, samp_global);
            [itile_line, itile_samp, line, samp] = tile_inv_pix(tile, line_global, samp_global);
    end
    
    % Inverse mapping to lat/lon
    [lon_rad, lat_rad] = tile_inv(tile, itile_line, itile_samp, line, samp);
    lat = lat_rad * 180 / pi;
    lon = lon_rad * 180 / pi;
    fprintf('  lat %.6f  long %.6f\n', lat, lon);
else
    error('Invalid direction. Use "fwd" or "inv"');
end

end

function proj_info = get_projection_info(projection)
% Get projection information structure

proj_names = {'is_k', 'is_h', 'is_q', 'gh', 'np', 'sp', 'sn', 'ha'};
proj_idx = find(strcmp(projection, proj_names));

if isempty(proj_idx)
    error('Unknown projection: %s', projection);
end

% Projection constants (simplified structure)
proj_info.name = projection;
proj_info.idx = proj_idx - 1; % Convert to 0-based indexing
proj_info.sphere = 6371007.181; % Earth radius in meters

% Define projection-specific parameters
switch projection
    case {'is_k', 'is_h', 'is_q', 'sn'}
        proj_info.ul_x = -20015109.354;
        proj_info.ul_y = 10007554.677;
        proj_info.pixel_size = 926.62543305;
        proj_info.nl_tile = 1200;
        proj_info.ns_tile = 1200;
        proj_info.ntile_line = 18;
        proj_info.ntile_samp = 36;
        
    case 'gh'
        proj_info.ul_x = -20015500.0;
        proj_info.ul_y = 8673500.0;
        proj_info.pixel_size = 1000.0;
        proj_info.nl_tile = 964;
        proj_info.ns_tile = 1112;
        proj_info.ntile_line = 18;
        proj_info.ntile_samp = 36;
        
    case {'np', 'sp'}
        proj_info.ul_x = -9058902.1845;
        proj_info.ul_y = 9058902.1845;
        proj_info.pixel_size = 1002.701;
        proj_info.nl_tile = 951;
        proj_info.ns_tile = 951;
        proj_info.ntile_line = 19;
        proj_info.ntile_samp = 19;
        
    case 'ha'
        proj_info.ul_x = -18020554.088;
        proj_info.ul_y = 9010277.044;
        proj_info.pixel_size = 997.15328068;
        proj_info.nl_tile = 1004;
        proj_info.ns_tile = 1004;
        proj_info.ntile_line = 18;
        proj_info.ntile_samp = 36;
end

end

function pixel_ratio = get_pixel_ratio(pixel_size)
% Convert pixel size character to ratio

switch pixel_size
    case 'k'
        pixel_ratio = 1;
    case 'h'
        pixel_ratio = 2;
    case 'q'
        pixel_ratio = 4;
    case 's'
        pixel_ratio = -1;
    case 'm'
        pixel_ratio = -2;
    case 'l'
        pixel_ratio = -3;
    case 'x'
        pixel_ratio = -4;
    otherwise
        error('Unknown pixel size: %s', pixel_size);
end

end

function tile = tile_init(proj_info, pixel_ratio)
% Initialize tile structure

tile.proj_info = proj_info;
tile.pixel_ratio = pixel_ratio;

if pixel_ratio > 0
    tile.nl = proj_info.nl_tile * proj_info.ntile_line * pixel_ratio;
    tile.ns = proj_info.ns_tile * proj_info.ntile_samp * pixel_ratio;
    tile.nl_tile = proj_info.nl_tile * pixel_ratio;
    tile.ns_tile = proj_info.ns_tile * pixel_ratio;
    tile.siz_x = proj_info.pixel_size / pixel_ratio;
    tile.siz_y = proj_info.pixel_size / pixel_ratio;
    tile.nl_offset = 0;
    tile.ns_offset = 0;
else
    % Global grid ratios (simplified)
    global_ratios = [240, 40, 20, 5];
    ratio_idx = abs(pixel_ratio);
    factor = global_ratios(ratio_idx);
    
    tile.nl = proj_info.nl_tile * proj_info.ntile_line / factor;
    tile.ns = proj_info.ns_tile * proj_info.ntile_samp / factor;
    tile.nl_tile = tile.nl;
    tile.ns_tile = tile.ns;
    tile.siz_x = (2.0 * abs(proj_info.ul_x)) / tile.ns;
    tile.siz_y = (2.0 * proj_info.ul_y) / tile.nl;
    tile.nl_offset = 0;
    tile.ns_offset = 0;
end

tile.ul_x = proj_info.ul_x + (0.5 * tile.siz_x);
tile.ul_y = proj_info.ul_y - (0.5 * tile.siz_y);
tile.nl_p = tile.nl - (2 * tile.nl_offset);
tile.ns_p = tile.ns - (2 * tile.ns_offset);

% Initialize integerized sinusoidal parameters if needed
if contains(proj_info.name, 'is_')
    switch proj_info.name
        case 'is_k'
            nrow_half = 180 * 60;
        case 'is_h'
            nrow_half = 180 * 60 * 2;
        case 'is_q'
            nrow_half = 180 * 60 * 4;
    end
    tile.isinu = isinu_init(nrow_half, proj_info.sphere);
end

end

function isinu = isinu_init(nrow_half, sphere)
% Initialize integerized sinusoidal structure

isinu.nrow_half = nrow_half;
isinu.nrow = nrow_half * 2;
isinu.sphere_inv = 1.0 / sphere;
isinu.ang_size_inv = double(isinu.nrow) / pi;

% Pre-compute column information for each row
isinu.icol_cen = zeros(nrow_half, 1);
isinu.ncol_inv = zeros(isinu.nrow, 1);

for irow = 1:nrow_half
    % Calculate latitude at center of row
    clat = (pi/2) * (1.0 - (irow - 1 + 0.5) / nrow_half);
    
    % Calculate number of columns per row
    ncol = round(2.0 * cos(clat) * isinu.nrow);
    if ncol < 1
        ncol = 1;
    end
    
    % Save center column and inverse of number of columns
    isinu.icol_cen(irow) = floor((ncol + 1) / 2);
    isinu.ncol_inv(irow) = 1.0 / double(ncol);
end

% Set inverse column distance at equator
ncol_equator = round(2.0 * isinu.nrow);
isinu.col_dist_inv = ncol_equator / (2 * pi * sphere);

end

function [itile_line, itile_samp, line, samp] = tile_fwd(tile, lon, lat)
% Forward mapping from lat/lon to tile coordinates

% Map projection (simplified - actual implementation would use specific projections)
[x, y] = forward_project(tile.proj_info, lon, lat);

% Convert to tile coordinates
[itile_line, itile_samp, line, samp] = tile_inv_map(tile, x, y);

end

function [x, y] = forward_project(proj_info, lon, lat)
% Forward map projection (simplified implementation)

switch proj_info.name
    case {'is_k', 'is_h', 'is_q'}
        % Integerized Sinusoidal projection
        x = proj_info.sphere * lon * cos(lat);
        y = proj_info.sphere * lat;
        
    case 'sn'
        % Sinusoidal projection
        x = proj_info.sphere * lon * cos(lat);
        y = proj_info.sphere * lat;
        
    case 'gh'
        % Goode's Homolosine (simplified)
        if abs(lat) > 40.44 * pi/180
            % Use Mollweide for high latitudes
            x = proj_info.sphere * lon * cos(lat);
            y = proj_info.sphere * lat;
        else
            % Use Sinusoidal for low latitudes
            x = proj_info.sphere * lon * cos(lat);
            y = proj_info.sphere * lat;
        end
        
    case {'np', 'sp'}
        % Lambert Azimuthal Equal Area
        if strcmp(proj_info.name, 'np')
            % North polar
            lat0 = pi/2;
        else
            % South polar
            lat0 = -pi/2;
        end
        
        k = sqrt(2 / (1 + sin(lat0) * sin(lat) + cos(lat0) * cos(lat) * cos(lon)));
        x = proj_info.sphere * k * cos(lat) * sin(lon);
        y = proj_info.sphere * k * (cos(lat0) * sin(lat) - sin(lat0) * cos(lat) * cos(lon));
        
    case 'ha'
        % Hammer-Aitoff projection
        z = sqrt(1 + cos(lat) * cos(lon/2));
        x = proj_info.sphere * 2 * sqrt(2) * cos(lat) * sin(lon/2) / z;
        y = proj_info.sphere * sqrt(2) * sin(lat) / z;
        
    otherwise
        error('Unsupported projection: %s', proj_info.name);
end

end

function [lon, lat] = inverse_project(proj_info, x, y)
% Inverse map projection (simplified implementation)

switch proj_info.name
    case {'is_k', 'is_h', 'is_q', 'sn'}
        % Sinusoidal projection
        lat = y / proj_info.sphere;
        if abs(cos(lat)) < 1e-10
            lon = 0;
        else
            lon = x / (proj_info.sphere * cos(lat));
        end
        
    case 'gh'
        % Goode's Homolosine (simplified)
        lat = y / proj_info.sphere;
        if abs(cos(lat)) < 1e-10
            lon = 0;
        else
            lon = x / (proj_info.sphere * cos(lat));
        end
        
    case {'np', 'sp'}
        % Lambert Azimuthal Equal Area
        rho = sqrt(x^2 + y^2);
        if rho < 1e-10
            if strcmp(proj_info.name, 'np')
                lat = pi/2;
            else
                lat = -pi/2;
            end
            lon = 0;
        else
            c = 2 * asin(rho / (2 * proj_info.sphere));
            if strcmp(proj_info.name, 'np')
                lat = asin(cos(c) * sin(pi/2) + y * sin(c) * cos(pi/2) / rho);
            else
                lat = asin(cos(c) * sin(-pi/2) + y * sin(c) * cos(-pi/2) / rho);
            end
            lon = atan2(x * sin(c), rho * cos(c));
        end
        
    case 'ha'
        % Hammer-Aitoff projection
        z = sqrt(1 - (x/(2*proj_info.sphere))^2 - (y/(2*proj_info.sphere))^2);
        lat = asin(z * y / proj_info.sphere);
        lon = 2 * atan2(z * x, 2 * proj_info.sphere^2 - x^2 - y^2);
        
    otherwise
        error('Unsupported projection: %s', proj_info.name);
end

end

function [itile_line, itile_samp, line, samp] = tile_inv_map(tile, x, y)
% Convert map coordinates to tile coordinates

line_global = ((tile.ul_y - y) / tile.siz_y) - tile.nl_offset;
samp_global = ((x - tile.ul_x) / tile.siz_x) - tile.ns_offset;

[itile_line, itile_samp, line, samp] = tile_inv_pix(tile, line_global, samp_global);

end

function [itile_line, itile_samp, line, samp] = tile_inv_pix(tile, line_global, samp_global)
% Convert global pixel coordinates to tile coordinates

% Check bounds
if line_global < -0.5 || line_global > (tile.nl - 0.5) || ...
   samp_global < -0.5 || samp_global > (tile.ns - 0.5)
    error('Coordinates out of bounds');
end

% Convert to integer coordinates
iline_global = round(line_global + 0.5);
if iline_global >= tile.nl_p
    iline_global = tile.nl - 1;
end
if iline_global < 0
    iline_global = 0;
end

isamp_global = round(samp_global + 0.5);
if isamp_global >= tile.ns_p
    isamp_global = tile.ns - 1;
end
if isamp_global < 0
    isamp_global = 0;
end

% Calculate tile indices
itile_line = floor(iline_global / tile.nl_tile);
itile_samp = floor(isamp_global / tile.ns_tile);

% Calculate within-tile coordinates
line = line_global - (itile_line * tile.nl_tile);
samp = samp_global - (itile_samp * tile.ns_tile);

end

function [lon, lat] = tile_inv(tile, itile_line, itile_samp, line, samp)
% Inverse mapping from tile coordinates to lat/lon

% Convert to map coordinates
[x, y] = tile_fwd_map(tile, itile_line, itile_samp, line, samp);

fprintf('  x %.6f  y %.6f\n', x, y);

% Inverse project to lat/lon
[lon, lat] = inverse_project(tile.proj_info, x, y);

end

function [x, y] = tile_fwd_map(tile, itile_line, itile_samp, line, samp)
% Convert tile coordinates to map coordinates

[line_global, samp_global] = tile_fwd_pix(tile, itile_line, itile_samp, line, samp);

y = tile.ul_y - ((line_global + tile.nl_offset) * tile.siz_y);
x = tile.ul_x + ((samp_global + tile.ns_offset) * tile.siz_x);

end

function [line_global, samp_global] = tile_fwd_pix(tile, itile_line, itile_samp, line, samp)
% Convert tile coordinates to global pixel coordinates

% Check bounds
if line < -0.5 || line > (tile.nl_tile - 0.5) || ...
   samp < -0.5 || samp > (tile.ns_tile - 0.5)
    error('Tile coordinates out of bounds');
end

line_global = line + (tile.nl_tile * itile_line);
samp_global = samp + (tile.ns_tile * itile_samp);

if line_global < -0.5 || line_global > (tile.nl_p - 0.5) || ...
   samp_global < -0.5 || samp_global > (tile.ns_p - 0.5)
    error('Global coordinates out of bounds');
end

end

function print_usage()
% Print usage information

fprintf('Usage:\n');
fprintf('  tilemap3(projection, pixel_size, direction, type, point_data)\n');
fprintf('     where projection: is_k, is_h, is_q, gh, np, sp, sn, ha\n');
fprintf('           pixel_size: k, h, q, s, m, l, x\n');
fprintf('           direction: fwd, inv\n');
fprintf('           type: tp, gm, gp\n');
fprintf('           point_data: lat, lon (forward mapping)\n');
fprintf('                      vert_tile, horiz_tile, line, samp (inverse tile pixel)\n');
fprintf('                      x, y (inverse global map)\n');
fprintf('                      line, samp (inverse global pixel)\n');
fprintf('\n');
fprintf('Examples:\n');
fprintf('  tilemap3(''is_q'', ''k'', ''fwd'', ''tp'', 38.53, -77.0)\n');
fprintf('  tilemap3(''is_q'', ''k'', ''inv'', ''tp'', 5, 11, 175.90, 1171.23)\n');

end
