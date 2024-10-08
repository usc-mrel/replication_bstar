% demo_remove_bulk_motion.m

%% Clean slate
close all; clear all; clc;

%% Read a .cfl file
cfl_file = 'D:\replication_bstar_v2\data\vol0977_20240919\meas_MID02046_FID23506_pulseq_bstar_fb_@_iso\k0';
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
k0 = real(readcfl(cfl_file));
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));



%% Perform data binning for respiratory motion
margin = 10;
bins = prctile(resp, linspace(0 + margin, 100 - margin, B + 1).');
bidx = cell(B,1);
for b = 1:B
    idx = find((resp >= bins(b)) & (resp < bins(b + 1)));
    bidx{b} = idx;
end

%%
nr_radial_views = size(resp);
TR = 2.3e-3;
t = (0:nr_radial_views-1).' * TR;

diff_resp = diff(cat(1, resp, resp(end)));
positive_slope = (diff_resp > 0);
nagative_slope = (diff_resp <= 0);

figure;
hold on;
plot(t(positive_slope), resp(positive_slope), '.');
plot(t(nagative_slope), resp(nagative_slope), '.');