function [phi,theta] = calculate_spiral_phyllotaxis_pattern(nr_readouts, nr_interleaves, flag_self)
% calculate_spiral_phyllotaxis_pattern -- Calculate an adapted spiral phyllotaxis pattern for center-out readouts
%
%  Usage
%    [phi,theta] = calculate_spiral_phyllotaxis_pattern(nr_readouts, nr_interleaves, flag_self)
%  Inputs
%   nr_readouts        number of readouts (spokes) per interleave
%   nr_interleaves     number of interleaves
%   flag_self          flag for SI projections
%  Outputs
%   phi                a list of azimuthal angles [rad]
%   theta              a list of polar angles [rad]
%
%  Description
%    Spiral phyllotaxis pattern described by Delacoste et al., 2018 MRM
%    "A Double Echo Ultra Short Echo Time (UTE) Acquisition for Respiratory
%    Motion-Suppressed High Resolution Imaging of the Lung"
%    https://onlinelibrary.wiley.com/doi/10.1002/mrm.26891
%
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/08/2022, Last modified: 08/08/2022

%% Calculate the number of projections (spokes + SI projection)
if flag_self
    nr_projections_per_frame = nr_readouts + 1;
else
    nr_projections_per_frame = nr_readouts;
end
N = nr_readouts * nr_interleaves; % total number of readouts

%% Calculate the adapted spiral phyllotaxis trajectory for center-out readouts
phi_gold = 137.51; % golden angle increment [degree]

phi   = zeros(nr_projections_per_frame * nr_interleaves, 1, 'double'); % azimuthal angle [rad]
theta = zeros(nr_projections_per_frame * nr_interleaves, 1, 'double'); % polar angle [rad]

n = 1;
for p = 1:nr_projections_per_frame % number of spokes + one SI projection
    for m = 1:nr_interleaves % M interleaves
        %------------------------------------------------------------------
        % Calculate the linear index
        %------------------------------------------------------------------
        index = p + (m - 1) * nr_projections_per_frame;

        if (flag_self && (p == 1)) % the first projection of every interleaf
            phi(index) = 0;
            theta(index) = 0;
        else
            %--------------------------------------------------------------
            % Calculate the azimuthal angle (phi) of a readout
            %--------------------------------------------------------------
            phi(index) = (2 * pi) / 360 * n * phi_gold;

            %--------------------------------------------------------------
            % Calculate the polar angle (theta) of a readout
            %--------------------------------------------------------------
            if n < N / 2
                theta(index) = pi / sqrt(2) * sqrt(n / N);
            else
                theta(index) = pi - pi / sqrt(2) * sqrt((N - n) / N);
            end
            n = n + 1;
        end
    end
end

end