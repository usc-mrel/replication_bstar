function resp = estimate_resp(dc, tr, n, fl, fh, fw, flag_self)
% Original source code from
% Dr.Frank Ong's extreme_mri Python package: https://github.com/mikgroup/extreme_mri
% Dr.Li Feng's XD-GRASP MATLAB package: https://cai2r.net/resources/xd-grasp-matlab-code
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 05/27/2023, Last modified: 05/27/2023

%% Perform trajectory-dependent frequency component filtering
%--------------------------------------------------------------------------
% Dr.Frank Ong's extreme_mri Python package: https://github.com/mikgroup/extreme_mri
%
% Estimate respiratory signal from DC.
%
% The function performs:
% 1) Filter DC with a band-pass filter with symmetric extension.
% 2) Normalize each channel by a robust estimation of mean and variance.
% 3) Return the channel with the maximum variance.
%
% Args:
%     dc (array): multi-channel DC array of shape [num_coils, num_tr].
%     tr (float): TR in seconds.
%     n (int): length of band-pass filter
%     fl (float): lower cut-off of band-pass filter in Hz.
%     fh (float): higher cut-off of band-pass filter in Hz.
%     fw (float): transition width of band-pass filter in Hz.
%
% Returns:
%     array: respiratory signal of length num_tr
%--------------------------------------------------------------------------
% For digital filters, the cutoff frequencies must lie between 0 and 1, 
% where 1 corresponds to the Nyquist rate—half the sample rate or pi rad/sample.
if flag_self
    [NiNs,Nc] = size(dc);
    dc_filtered = dc;
else
    [NiNs,Nc] = size(dc);
    dc = abs(dc);
    fs = 1 / tr;
    bands = [0, fl - fw, fl, fh, fh + fw, fs / 2] / (fs / 2);
    desired = [0, 0, 1, 1, 0, 0];

    filt = firls(n-1, bands, desired);
    dc_filtered = zeros(NiNs, Nc, 'single');
    for c = 1:Nc
        dc_pad = cat(1, flip(dc(1+(1:floor(n/2)),c),1), dc(:,c), dc(end-(1:floor(n/2)),c));
        resp_c = conv(dc_pad, filt, 'valid');
        sigma_c = 1.4826 * median(abs(resp_c - median(resp_c)));
        dc_filtered(:,c) = (resp_c - median(resp_c)) / sigma_c;
    end    
end

%% Respiratory motion detection
%--------------------------------------------------------------------------
% Source code obtained from Dr.Li Feng's XD-GRASP MATLAB package
% https://cai2r.net/resources/xd-grasp-matlab-code/
%--------------------------------------------------------------------------
ZIP = reshape(dc_filtered, [1 NiNs Nc]); % # of readout points x # of TRs x # of channels

%--------------------------------------------------------------------------
% Perform PCA on each coil element
%--------------------------------------------------------------------------
kk = 1;
for c = 1:Nc
    tmp = permute(ZIP(:,:,c), [1 3 2]); % # of readout points x 1 x # of TRs
    tmp = reshape(tmp, [1 NiNs]).'; % # of TRs x # of readout points

    covariance = cov(tmp);
    [tmp2, V] = eig(covariance);
    V = diag(V);
    [junk, rindices] = sort(-1 * V);
    V = V(rindices);
    tmp2 = tmp2(:,rindices);
    PC = (tmp2.' * tmp.').';

    % Take the first two principal components from each coil element
    for jj = 1:1
        tmp3 = smooth(PC(:,jj), 6, 'lowess'); % do some moving average smoothing

        % Normalize the signal for display
        tmp3 = tmp3 - min(tmp3(:));
        tmp3 = tmp3 ./ max(tmp3(:));
        PCs(:,kk) = tmp3;
        kk = kk + 1;
    end
end

%--------------------------------------------------------------------------
% Do coil clustering to find the respiratory motion signal
% Function obtained from Tao Zhang (http://web.stanford.edu/~tzhang08/software.html)
%--------------------------------------------------------------------------
thresh = 0.95;
[resp, cluster] = CoilClustering(PCs, thresh);

%--------------------------------------------------------------------------
% Change the sign of a signal
%--------------------------------------------------------------------------
if resp(2) < resp(1) % descending
    resp = -resp;
end

%--------------------------------------------------------------------------
% Normalize the signal for display
%--------------------------------------------------------------------------
resp = resp - min(resp(:));
resp = resp ./ max(resp(:));

end