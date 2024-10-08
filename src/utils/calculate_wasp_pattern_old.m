function [phi, theta] = calculate_wasp_pattern(N, M, rpa_range)
% Inputs
%   N            number of half-radial projections per interleaf
%   M            number of randomly tilted interleaves
%   rpa_range    range of random polar angles in degree
% Outputs
%   phi          list of azimuthal angles [rad]
%   theta        list of polar angles [rad]


%% Calculate WASP trajectories (wobbling Archemedean spiral pole)
phi_gold = pi * (3 + sqrt(5)); % golden angle [rad]

x = zeros(N, M, 'double');
y = zeros(N, M, 'double');
z = zeros(N, M, 'double');
n = (1:N).';

for m = 1:M
    %----------------------------------------------------------------------
    % Calculate the z coordinate of a readout
    %----------------------------------------------------------------------
    z(:,m) = (2 * n - N - 1) / N; % = cos(theta)

    %----------------------------------------------------------------------
    % Calculate the azimuthal angle of a readout [rad]
    %----------------------------------------------------------------------
    phi = (sqrt(N * pi / M) * asin(z(:,m)) + 2 * m * pi / M); % [rad]

    %----------------------------------------------------------------------
    % Rotate a spiral iterleaf about the polar axis by the golden angle
    %----------------------------------------------------------------------
    phi = phi + (m - 1) * phi_gold;

    %----------------------------------------------------------------------
    % Calculate the magnitude of a sample in the azimuthal plane
    %----------------------------------------------------------------------
    sin_theta = sqrt(1 - z(:,m).^2);

    %----------------------------------------------------------------------
    % Calculate the x coordinate of a readout
    %----------------------------------------------------------------------
    x(:,m) = cos(phi) .* sin_theta;

    %----------------------------------------------------------------------
    % Calculate the y coordinate of a readout
    %----------------------------------------------------------------------
    y(:,m) = sin(phi) .* sin_theta;
end

%% Calculate a small random polar angle [degree]
rng('default');
random_angles = rpa_range(1) + (rpa_range(2) - rpa_range(1)) * rand(M,1); % [degree]

%% Tilt a spiral interleaf by a small random polar angle
%--------------------------------------------------------------------------
%         ^ z(3)                       ^ z(3) (theta on the yz-plane)
%         |                     theta  |
%         |                          \-|
%         |                           \|
%         +---------> y(2) =>          +---------> y(2)
%        /                            /
%       /                            /
%      v x(1)                       v x(1)
%--------------------------------------------------------------------------
for m = 1:M
    %----------------------------------------------------------------------
    % Calculate a rotation matrix by theta about the x-axis
    %----------------------------------------------------------------------
    theta = pi / 180 * random_angles(m);
    Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];

    %----------------------------------------------------------------------
    % Rotate a spiral interleaf about the x-axis
    %----------------------------------------------------------------------
    xyz = cat(2, x(:,m), y(:,m), z(:,m)); % N x 3
    xyz_rotated = (Rx * xyz.').'; % N x 3

    %----------------------------------------------------------------------
    % Update the spatial coordinates
    %----------------------------------------------------------------------
    x(:,m) = xyz_rotated(:,1);
    y(:,m) = xyz_rotated(:,2);
    z(:,m) = xyz_rotated(:,3);
end

%% Reverse the traversing direction of even-numbered interleaves
for m = 1:M
    if mod(m,2) == 0 % even
        x(:,m) = x(end:-1:1,m);
        y(:,m) = y(end:-1:1,m);
        z(:,m) = z(end:-1:1,m);
    end
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