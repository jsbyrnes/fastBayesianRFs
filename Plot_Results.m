%%%%%
%plot receiver functions
clear
close all
clc

files = dir('./3J_localevents100Hz/*.rjMCMC.mat');

scale = 50;

figure(1)
hold on
set(gca, 'YDir', 'reverse')
ylim([ -0.01 0.15 ])

for k = 1:length(files)

    load([ './3J_localevents100Hz/' files(k).name ], 'rf_rj', 'rf_std', 'allWfs')

    if max(rf_std)/max(abs(rf_rj)) > 8

        continue

    end

    t = (1:length(rf_rj))/300 - 1;

    %rf = rf_td(:, k);
    %rf = rf_mt(:, k);
    %rf = rf_rj(k, :);
    %rf = 3*(1/length(allWfs))*scale*rf;

    amp = max(abs(rf_rj));

    plot(scale*rf_rj/amp + allWfs.delta*111.12*1e3 - 291694, t, 'k')
    plot(scale*(rf_rj + rf_std/2)/amp + allWfs.delta*111.12*1e3 - 291694, t, 'k--')
    plot(scale*(rf_rj - rf_std/2)/amp + allWfs.delta*111.12*1e3 - 291694, t, 'k--')

end
