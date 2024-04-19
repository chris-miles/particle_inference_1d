clear all
close all

load('postareas_src.mat')

bins = linspace(30,110,35);
subplot(4,1,1);
histogram(areas0,bins,'normalization','pdf','EdgeAlpha',0);
xline(mean(areas0))
xline(mean(areas0)-std(areas0))
xline(mean(areas0)+std(areas0))
ylim([0, 0.1])

box off;
xlim([30,110]);


subplot(4,1,2);
histogram(areas,bins,'normalization','pdf','EdgeAlpha',0);
xline(mean(areas))
xline(mean(areas)-std(areas0))
xline(mean(areas)+std(areas))

box off;
xlim([30,110]);
ylim([0, 0.1])


load('postareas_boundary.mat')

subplot(4,1,3);
histogram(areas,bins,'normalization','pdf','EdgeAlpha',0);
xline(mean(areas))
xline(mean(areas)-std(areas0))
xline(mean(areas)+std(areas))
box off;
xlim([30,110]);
ylim([0, 0.1])



load('postareas_both.mat')

subplot(4,1,4);
histogram(areas,bins,'normalization','pdf','EdgeAlpha',0);
xline(mean(areas))
xline(mean(areas)-std(areas0))
xline(mean(areas)+std(areas))
box off;
xlim([30,110]);
ylim([0, 0.1])

