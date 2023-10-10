% display_trj_file.m

figure('Color', 'w', 'Position', [-2 246 1604 552]);
subplot(2,6,1); hold on;
plot(trj(:,1), 'LineWidth', 1); grid on;
set(gca, 'Box', 'On');
xlabel('Sample index');
ylabel('Amplitude [cycle/m]');
title('+X, echo 1');
legend('USC', 'Location', 'northwest');

subplot(2,6,2); hold on;
plot(trj(:,2), 'LineWidth', 1); grid on;
set(gca, 'Box', 'On');
xlabel('Sample index');
ylabel('Amplitude [cycle/m]');
title('-X, echo 1');
legend('USC', 'Location', 'northwest');

subplot(2,6,3); hold on;
plot(trj(:,3), 'LineWidth', 1); grid on;
set(gca, 'Box', 'On');
xlabel('Sample index');
ylabel('Amplitude [cycle/m]');
title('+Y, echo 1');
legend('USC', 'Location', 'northwest');

subplot(2,6,4); hold on;
plot(trj(:,4), 'LineWidth', 1); grid on;
set(gca, 'Box', 'On');
xlabel('Sample index');
ylabel('Amplitude [cycle/m]');
title('-Y, echo 1');
legend('USC', 'Location', 'northwest');

subplot(2,6,5); hold on;
plot(trj(:,5), 'LineWidth', 1); grid on;
set(gca, 'Box', 'On');
xlabel('Sample index');
ylabel('Amplitude [cycle/m]');
title('+Z, echo 1');
legend('USC', 'Location', 'northwest');

subplot(2,6,6); hold on;
plot(trj(:,6), 'LineWidth', 1); grid on;
set(gca, 'Box', 'On');
xlabel('Sample index');
ylabel('Amplitude [cycle/m]');
title('-Z, echo 1');
legend('USC', 'Location', 'northwest');

subplot(2,6,7); hold on;
plot(trj(:,7), 'LineWidth', 1); grid on;
set(gca, 'Box', 'On');
xlabel('Sample index');
ylabel('Amplitude [cycle/m]');
title('+X, echo 2');
legend('USC', 'Location', 'northeast');

subplot(2,6,8); hold on;
plot(trj(:,8), 'LineWidth', 1); grid on;
set(gca, 'Box', 'On');
xlabel('Sample index');
ylabel('Amplitude [cycle/m]');
title('-X, echo 2');
legend('USC', 'Location', 'northeast');

subplot(2,6,9); hold on;
plot(trj(:,9), 'LineWidth', 1); grid on;
set(gca, 'Box', 'On');
xlabel('Sample index');
ylabel('Amplitude [cycle/m]');
title('+Y, echo 2');
legend('USC', 'Location', 'northeast');

subplot(2,6,10); hold on;
plot(trj(:,10), 'LineWidth', 1); grid on;
set(gca, 'Box', 'On');
xlabel('Sample index');
ylabel('Amplitude [cycle/m]');
title('-Y, echo 2');
legend('USC', 'Location', 'northeast');

subplot(2,6,11); hold on;
plot(trj(:,11), 'LineWidth', 1); grid on;
set(gca, 'Box', 'On');
xlabel('Sample index');
ylabel('Amplitude [cycle/m]');
title('+Z, echo 2');
legend('USC', 'Location', 'northeast');

subplot(2,6,12); hold on;
plot(trj(:,12), 'LineWidth', 1); grid on;
set(gca, 'Box', 'On');
xlabel('Sample index');
ylabel('Amplitude [cycle/m]');
title('-Z, echo 2');
legend('USC', 'Location', 'northeast');

% sgtitle({sprintf('Measured k-space trajectories, bandwith = %4.0f [Hz/Px], samples = %d', bandwidth, nr_samples), ...
%          sprintf('RUTime/TotalTime/ADC[0]/Resolution = %3.0f/%3.0f/%6.1f/%4.2f', ramp_time * 1e6, total_time * 1e6, adc_duration_bstar * 1e6, resolution * 1e3), ...
%          sprintf('scope/bstar amplitude = %3.1f/%3.1f [mT/m], scale factor = %6.4f', amplitude_scope, amplitude_bstar, scale_factor)});

output_fullpath = fullfile(output_directory, sprintf('traj_b%4.0f_n%d', bandwidth, nr_samples));
export_fig(output_fullpath, '-r300', '-tif');

close(gcf);
