function [phi, theta] = calculate_wasp_pattern(S, M, rpa_range)
% Inputs
%   S            number of half-radial projections per "full" interleaf ("full" = Up+down Archimedean spiral)
%   M            number of randomly tilted interleaves
%   rpa_range    range of random polar angles in degree
% Outputs
%   phi          list of azimuthal angles [rad]
%   theta        list of polar angles [rad]
% This is the correct version of WASP trajectory!


%% Calculate the number of half-radial projections per "half" interleaf ("full" = Up+down Archimedean spiral)
N = S / 2;

%% Calculate sample points on a unit sphere (up Archimedean spiral)
n = (1:N).';
z_up = (2 * n - N - 1) / N; % = cos(theta)
phi_up = sqrt(N * pi) * asin(z_up);
x_up = cos(phi_up) .* sqrt(1 - z_up.^2);
y_up = sin(phi_up) .* sqrt(1 - z_up.^2);

%% Calculate sample points on a unit sphere (down Archimedean spiral)
z_down = (N + 1 - 2 * n) / N; % = cos(theta)
phi_down = sqrt(N * pi) * asin(z_down) + pi;
x_down = cos(phi_down) .* sqrt(1 - z_down.^2);
y_down = sin(phi_down) .* sqrt(1 - z_down.^2);

%% Calculate a full interleaf
x_full_interleaf = cat(1, x_up, x_down);
y_full_interleaf = cat(1, y_up, y_down);
z_full_interleaf = cat(1, z_up, z_down);

%% Calculate a small random polar angle [degree]
rng('default');
polar_angles = rpa_range(1) + (rpa_range(2) - rpa_range(1)) * rand(M,1); % [degree]

azimuthal_range = [0 360]; % range of random azimuthal angles in degree 
azimuthal_angles = azimuthal_range(1) + (azimuthal_range(2) - azimuthal_range(1)) * rand(M,1); % [degree]

x_polar_axis = cos(azimuthal_angles * pi / 180) .* sin(polar_angles * pi / 180);
y_polar_axis = sin(azimuthal_angles * pi / 180) .* sin(polar_angles * pi / 180);
z_polar_axis = cos(polar_angles * pi / 180);

%% Calculate WASP trajectories (wobbling Archemedean spiral pole)
phi_gold = pi * (3 + sqrt(5)); % golden angle [rad]

x = zeros(S, M, 'double');
y = zeros(S, M, 'double');
z = zeros(S, M, 'double');

xyz_full_interleaf = cat(1, x_full_interleaf.', y_full_interleaf.', z_full_interleaf.'); % 3 x N

for m = 1:M
    %----------------------------------------------------------------------
    % Calculate a rotation matrix by phi about the axis u = (ux,uy,uz)
    %----------------------------------------------------------------------
    ux = x_polar_axis(m);
    uy = y_polar_axis(m);
    uz = z_polar_axis(m);

    phi = (m - 1) * phi_gold;

    u_cross = [ 0   -uz   uy ;
                uz   0   -ux ;
               -uy   ux   0 ];

    uuT = [  ux^2    ux * uy   ux * uz ;
           ux * uy     uy^2    uy * uz ;
           ux * uz   uy * uz     uz^2 ];

    R = cos(phi) * eye(3,3) + sin(phi) * u_cross + (1 - cos(phi)) * uuT;

    %----------------------------------------------------------------------
    % Rotate the whole system by a golden angle about the the tilted polar axis
    %----------------------------------------------------------------------
    xyz_rotated = R * xyz_full_interleaf; % 3 x N

    %----------------------------------------------------------------------
    % Update the spatial coordinates
    %----------------------------------------------------------------------
    x(:,m) = xyz_rotated(1,:);
    y(:,m) = xyz_rotated(2,:);
    z(:,m) = xyz_rotated(3,:);
end

%% Convert the Cartesian coordinates (x,y,z) to the spherical coordinates (r, theta, phi) = (radius, polar, azimuth)
%--------------------------------------------------------------------------
% phi  : azimuthal angle
% theta: polar angle
%
%              ^ z
%              |
%              |\
%              | \
%              |  \
%              |   \  spherical coordinates (r,theta,phi)        
%              |    + Cartesian coordinates (x,y,z)
%        theta |\  /|
%              | \/ |
%              | /  |
%              |/   |
%              +----|-------------------> y
%             / \   |   /
%            /   \  |  /
%           /---->\ | /              r = sqrt(x^2 + y^2 + z^2)
%          /  phi  \|/               phi = atan(y / x)
%         /---------+                theta = atan(sqrt(x^2 + y^2) / z)
%        /
%       v x
%        
%--------------------------------------------------------------------------
phi   = atan2(y(:), x(:));                    % azimuth angle [rad]
theta = atan2(sqrt(x(:).^2 + y(:).^2), z(:)); % polar angle [rad]

end