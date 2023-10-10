% Demostration of coil clustering
% (c) Tao Zhang, 2015, Stanford University

load ButterflyNavigator.mat;

thresh = 0.95;
[navCC, cluster] = CoilClustering(navSI, thresh);

navAve = mean(navSI,2);
navManual = navSI(:,31);

figure,plot(t, navAve, 'r', t, navCC, 'b', t, navManual, 'k', 'LineWidth',2),axis([30 60 -2 4]),xlabel('time (s)'),ylabel('Motion (pixels)'),legend('Averaging','Coil Clustering','Manual Selection')

