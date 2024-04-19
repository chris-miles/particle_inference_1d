clear all
close all

addpath(genpath('cbrewer'));

sigmaplot = 7;
sigmavals = linspace(0.1, 0.9,sigmaplot);

zvals = linspace(0,1,301);
Lvals = linspace(0,3,301);

cc= cbrewer('seq','Greys',9);

figure;
subplot(1,2,1)


for ss=1:sigmaplot
    
  betavals = betapdf(zvals,1/sigmavals(ss),1/sigmavals(ss));
  plot(zvals, betavals,'LineWidth',3,'Color',cc(ss+2,:));
  hold on;
end 
box off
axis square;
axis tight; xlim([0 1]);

subplot(1,2,2)

for ss=1:sigmaplot
    
    gammavals = gampdf(Lvals,1/(sigmavals(ss)^2), (sigmavals(ss)^2)/1);
  plot(Lvals, gammavals,'LineWidth',3,'Color',cc(ss+2,:));
  hold on;
end 

box off
axis square;
axis tight; xlim([0 3]);