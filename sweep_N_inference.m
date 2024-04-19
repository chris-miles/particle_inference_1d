clear all
close all

rng(1)


    D_true =1;%
    kbirth_true = 500;
    kdeath_true = 10;% decay rate

    Tfinal = 100; % experiment length
    Xmax=1; % size of domain [0,R];


      Nexperiments = 100; % do 1000 experiments


    Xsource_mean = 0.5; % position of source, 1d position
    Xs_noise = 0;%0.1;

    Xmax_noise=0;%; 

    Nparticles = zeros(1,Nexperiments);
    data = cell(1, Nexperiments);
    zvals= zeros(1,Nexperiments);
    Lvals = zeros(1,Nexperiments);

    tic
    parfor n = 1:Nexperiments

        Xmax_n = Xmax+Xmax_noise*randn;
        if Xmax_n<0 
            Xmax_n=-Xmax_n;
        end 
        
        Xs_n = Xsource_mean + Xs_noise*randn;

        if Xs_n<0
            Xs_n = -Xs_n;
        elseif Xs_n>Xmax_n
            Xs_n = Xs_n-Xmax_n;
        end

        Lvals(n) = Xmax_n;
        zvals(n) = Xs_n;
    
        data{n} = montecarlosim_1d(Tfinal, Xmax_n, Xs_n, kbirth_true, kdeath_true, D_true);
    
        positions = data{n};
        Nparticles(n) = length(positions);
    end
toc

data10 = data(1:10);
Lvals10 = Lvals(1:10);
zvals10 = zvals(1:10);


data25 = data(1:25);
Lvals25 = Lvals(1:25);
zvals25 = zvals(1:25);


data50 = data(1:50);
Lvals50 = Lvals(1:50);
zvals50 = zvals(1:50);


data100 = data(1:100);
Lvals100 = Lvals(1:100);
zvals100 = zvals(1:100);


data250 = data(1:250);
Lvals250 = Lvals(1:250);
zvals250 = zvals(1:250);


data500 = data(1:500);
Lvals500 = Lvals(1:500);
zvals500 = zvals(1:500);


data1000 = data(1:1000);
Lvals1000 = Lvals(1:1000);
zvals1000 = zvals(1:1000);


sigmaInit = 1;
meanLambda = 250;
varLambda = meanLambda^2;

alphaLambda = meanLambda^2/varLambda;
betaLambda  =  meanLambda/varLambda;
muMax = 3*(kdeath_true)/D_true;

acceptThresh1=0.3;
acceptThresh2=0.5;

sWarmup=250;
steps=5e4;



[outsamps10,acceptProb10,postEvals10] = mhGibbsWithWarmup(data10,zvals10,...,
    Lvals10,alphaLambda,betaLambda, muMax, sigmaInit, steps,sWarmup,...
acceptThresh1,acceptThresh2);

[outsamps10_2,acceptProb10_2] = mhGibbsWithWarmup(data10,zvals10,...,
    Lvals10,alphaLambda,betaLambda, muMax, sigmaInit, steps,sWarmup,...
acceptThresh1,acceptThresh2);



[outsamps25,acceptProb25,postEvals25] = mhGibbsWithWarmup(data25,zvals25,...,
    Lvals25,alphaLambda,betaLambda, muMax, sigmaInit, steps,sWarmup,...
acceptThresh1,acceptThresh2);


[outsamps50,acceptProb50,postEvals50] = mhGibbsWithWarmup(data50,zvals50,...,
    Lvals50,alphaLambda,betaLambda, muMax, sigmaInit, steps,sWarmup,...
acceptThresh1,acceptThresh2);


[outsamps100,acceptProb100,postEvals100] = mhGibbsWithWarmup(data100,zvals100,...,
    Lvals100,alphaLambda,betaLambda, muMax, sigmaInit, steps,sWarmup,...
acceptThresh1,acceptThresh2);



[outsamps250,acceptProb250,postEvals250] = mhGibbsWithWarmup(data250,zvals250,...,
    Lvals250,alphaLambda,betaLambda, muMax, sigmaInit, steps,sWarmup,...
acceptThresh1,acceptThresh2);



[outsamps500,acceptProb500,postEvals500] = mhGibbsWithWarmup(data500,zvals500,...,
    Lvals500,alphaLambda,betaLambda, muMax, sigmaInit, steps,sWarmup,...
acceptThresh1,acceptThresh2);


[outsamps1000,acceptProb1000,postEvals1000] = mhGibbsWithWarmup(data1000,zvals1000,...,
    Lvals1000,alphaLambda,betaLambda, muMax, sigmaInit, steps,sWarmup,...
acceptThresh1,acceptThresh2);



addpath(genpath('cbrewer'));


alpha =.11 ;



q10 = quantile(postEvals10,[alpha, 1]);
hi10=outsamps10(postEvals10>q10(1),:);
[~,area10]=boundary(hi10);
bdry10 = hi10(boundary(hi10),:);

q25 = quantile(postEvals25,[alpha, 1]);
hi25=outsamps25(postEvals25>q25(1),:);
[~,area25]=boundary(hi25);

bdry25 = hi25(boundary(hi25),:);

q50 = quantile(postEvals50,[alpha, 1]);
hi50=outsamps50(postEvals50>q50(1),:);
[~,area50]=boundary(hi50);

bdry50 = hi50(boundary(hi50),:);

q100 = quantile(postEvals100,[alpha, 1]);
hi100=outsamps100(postEvals100>q100(1),:);
[~,area100]=boundary(hi100);
bdry100 = hi100(boundary(hi100),:);

q250 = quantile(postEvals250,[alpha, 1]);
hi250=outsamps250(postEvals250>q250(1),:);
[~,area250]=boundary(hi250);
bdry250 = hi250(boundary(hi250),:);

q500 = quantile(postEvals500,[alpha, 1]);
hi500=outsamps500(postEvals500>q500(1),:);
[~,area500]=boundary(hi500);
bdry500 = hi500(boundary(hi500),:);

q1000 = quantile(postEvals1000,[alpha, 1]);
hi1000=outsamps1000(postEvals1000>q1000(1),:);
[~,area1000]=boundary(hi1000);
bdry1000 = hi1000(boundary(hi1000),:);

Nvals =[10, 25, 50, 100, 250, 500, 1000];
areas = [area10, area25,area50,area100,area250,area500,area1000];
%%%%



addpath('cbrewer')
cc = flipud(cbrewer('seq','YlGnBu',9));


figure('position',[641 326 250 250]);


plot(Nvals,20000*Nvals.^(-1))
hold on;
for nnn=1:7
scatter(Nvals(nnn),areas(nnn),100,'filled','MarkerEdgeAlpha',0,'MarkerFaceColor',cc(nnn,:),'MarkerFaceAlpha',1)
hold on;
end 

set(gca,'YScale','log')
set(gca,'XScale','log')
hold on
plot(Nvals,20000*Nvals.^(-1))
axis square;
set(gca,'LineWidth',1.25)




figure;

subplot(1,3,1);




plot(bdry10(:,1),bdry10(:,2),'LineWidth',2.5,'Color',cc(1,:));
hold on;
%plot(bdry25(:,1),bdry25(:,2),'LineWidth',2.5,'Color',cc(2,:));

plot(bdry50(:,1),bdry50(:,2),'LineWidth',2.5,'Color',cc(3,:));
%plot(bdry100(:,1),bdry100(:,2),'LineWidth',2.5,'Color',cc(3,:));
%plot(bdry250(:,1),bdry250(:,2),'LineWidth',2.5,'Color',cc(4,:));
plot(bdry500(:,1),bdry500(:,2),'LineWidth',2.5,'Color',cc(6,:));
%plot(bdry1000(:,1),bdry1000(:,2),'LineWidth',2.5,'Color',cc(6,:));
scatter(kbirth_true/D_true,(kdeath_true/D_true),150,'+','LineWidth',2,'MarkerEdgeColor',[161,217,155]/256)

axis square;
ylim([0 17])
%xlim([0.5*kbirth_true/D_true 2*kbirth_true/D_true])
xlim([200 700])
%ylim([0.5*kdeath_true/D_true 2*kdeath_true/D_true])
xlabel('\lambda')
ylabel('\mu')
set(gca,'LineWidth',1.25)
%legend({'N=10','N=25','N=50','N=100','N=250','N=500','N=1000'},'position','bestoutside')
set(gca,'LineWidth',1.25)

subplot(1,3,2);

xline(kbirth_true,'LineWidth',2,'color',[161,217,155]/256)
hold on;

xvals = 0:.1:2*kbirth_true/D_true;
plot(xvals, gampdf(xvals,alphaLambda,betaLambda^(-1)),...,
    '--','LineWidth',2,'color',[.6,.6,.6])


hold on;
lambdabins = linspace(0,2*kbirth_true/D_true,23);


%histogram(outsamps25(:,1),lambdabins,'normalization','pdf','EdgeAlpha',0.0,...
%    'FaceColor',cc(2,:),'DisplayStyle','bar','FaceAlpha',0.3);
%histogram(outsamps50(:,1),lambdabins,'normalization','pdf','EdgeAlpha',0.0,...
%    'FaceColor',cc(3,:),'DisplayStyle','bar','FaceAlpha',0.3);

%histogram(outsamps250(:,1),lambdabins,'normalization','pdf','EdgeAlpha',0.0,...
%    'FaceColor',cc(5,:),'DisplayStyle','bar','FaceAlpha',0.3);
%histogram(outsamps500(:,1),lambdabins,'normalization','pdf','EdgeAlpha',0.0,...
%    'FaceColor',cc(6,:),'DisplayStyle','bar','FaceAlpha',0.3);
histogram(outsamps500(:,1),lambdabins,'normalization','pdf','EdgeAlpha',0.0,...
    'FaceColor',cc(6,:),'DisplayStyle','bar','FaceAlpha',.7);
histogram(outsamps50(:,1),lambdabins,'normalization','pdf','EdgeAlpha',0.0,...
    'FaceColor',cc(3,:),'DisplayStyle','bar','FaceAlpha',.7);
histogram(outsamps10(:,1),lambdabins,'normalization','pdf','EdgeAlpha',0.0,...
    'FaceColor',cc(1,:),'DisplayStyle','bar','FaceAlpha',.7);hold on;
%set(gca,'YScale','log')
axis tight;
axis square;
xlabel('\lambda')
ylabel('pdf')
legend({'true','prior','N=10','N=25','N=50','N=100','N=250','N=500','N=1000'},...
    'location','eastoutside')
legend box off;
set(gca,'LineWidth',1.25)


subplot(1,3,3);


xvals = .01:.1:2*(kdeath_true/D_true);
xline((kdeath_true/D_true),'LineWidth',2,'color',[161,217,155]/256)
hold on
plot(xvals, ones(size(xvals))/muMax,'--','LineWidth',2,'color',[.6,.6,.6])

hold on;
mubins = linspace(0,2*kdeath_true/D_true,23);

histogram(outsamps500(:,2),mubins,'normalization','pdf','EdgeAlpha',0.0,...
    'FaceColor',cc(6,:),'DisplayStyle','bar','FaceAlpha',0.7);

histogram(outsamps50(:,2),mubins,'Normalization','pdf','EdgeAlpha',0.0,...
    'FaceColor',cc(3,:),'DisplayStyle','bar','FaceAlpha',0.7);
histogram(outsamps10(:,2),mubins,'normalization','pdf','EdgeAlpha',0.0,...
    'FaceColor',cc(1,:),'DisplayStyle','bar','FaceAlpha',0.7);
xlim([0 2*kdeath_true/D_true])
%set(gca,'YScale','log')
axis tight;
axis square;

xlabel('\mu')
ylabel('pdf')
set(gca,'LineWidth',1.25)



%%%%%%%%%%%%%%%%%%%%
% analytical figures start here!!

muvals = linspace(0,2*kdeath_true/D_true,101);
lambdavals = linspace(0,2*kbirth_true/D_true,101);
for m = 1:length(muvals)

logL10mu(m) = logL(data10,zvals10,Lvals10,kbirth_true,muvals(m));
logL25mu(m) = logL(data25,zvals25,Lvals25,kbirth_true,muvals(m));
logL50mu(m) = logL(data50,zvals50,Lvals50,kbirth_true,muvals(m));
logL100mu(m) = logL(data100,zvals100,Lvals100,kbirth_true,muvals(m));
logL250mu(m) = logL(data250,zvals250,Lvals250,kbirth_true,muvals(m));
logL500mu(m) = logL(data500,zvals500,Lvals500,kbirth_true,muvals(m));
logL1000mu(m) = logL(data1000,zvals1000,Lvals1000,kbirth_true,muvals(m));

logL10lamb(m) = logL(data10,zvals10,Lvals10,lambdavals(m),(kdeath_true/D_true));
logL25lamb(m) = logL(data25,zvals25,Lvals25,lambdavals(m),(kdeath_true/D_true));
logL50lamb(m) = logL(data50,zvals50,Lvals50,lambdavals(m),(kdeath_true/D_true));
logL100lamb(m) = logL(data100,zvals100,Lvals100,lambdavals(m),(kdeath_true/D_true));
logL250lamb(m) = logL(data250,zvals250,Lvals250,lambdavals(m),(kdeath_true/D_true));
logL500lamb(m) = logL(data500,zvals500,Lvals500,lambdavals(m),(kdeath_true/D_true));
logL1000lamb(m) = logL(data1000,zvals1000,Lvals1000,lambdavals(m),(kdeath_true/D_true));

end 

figure('position',[411 549 525 225]);
subplot(1,2,1);
xline(kbirth_true,'LineWidth',2,'color',[161,217,155]/256)
hold on;
plot(lambdavals, logL10lamb-max(logL10lamb),'LineWidth',2,'color',cc(1,:))
hold on;
plot(lambdavals, logL25lamb-max(logL25lamb),'LineWidth',2,'color',cc(2,:))
plot(lambdavals, logL50lamb-max(logL50lamb),'LineWidth',2,'color',cc(3,:))
plot(lambdavals, logL100lamb-max(logL100lamb),'LineWidth',2,'color',cc(4,:))
plot(lambdavals, logL250lamb-max(logL250lamb),'LineWidth',2,'color',cc(5,:))
plot(lambdavals, logL500lamb-max(logL500lamb),'LineWidth',2,'color',cc(6,:))
plot(lambdavals, logL1000lamb-max(logL1000lamb),'LineWidth',2,'color',cc(7,:))
ylim([-100,0])
axis square;
box off
ylabel('$\log L(\lambda,\mu_\mathrm{true})$','interpreter','latex')
xlabel('\lambda')
set(gca,'LineWidth',1.25)

subplot(1,2,2);
xline((kdeath_true/D_true),'LineWidth',2,'color',[161,217,155]/256)
hold on;

plot(muvals, logL10mu-max(logL10mu),'LineWidth',2,'color',cc(1,:))
plot(muvals, logL25mu-max(logL25mu),'LineWidth',2,'color',cc(2,:))
plot(muvals, logL50mu-max(logL50mu),'LineWidth',2,'color',cc(3,:))
plot(muvals, logL100mu-max(logL100mu),'LineWidth',2,'color',cc(4,:))
plot(muvals, logL250mu-max(logL250mu),'LineWidth',2,'color',cc(5,:))
plot(muvals, logL500mu-max(logL500mu),'LineWidth',2,'color',cc(6,:))
plot(muvals, logL1000mu-max(logL1000mu),'LineWidth',2,'color',cc(7,:))
ylim([-100,0])
axis square;
box off
ylabel('$\log L(\lambda_\mathrm{true},\mu)$','interpreter','latex')
xlabel('\mu')
legend({'true','N=10','N=25','N=50','N=100','N=250','N=500','N=1000'},...
    'location','eastoutside')
legend box off;
set(gca,'LineWidth',1.25)





muvals = linspace(.001,2*kdeath_true/D_true,101);
lambdavals = linspace(0,2*kbirth_true/D_true,101);
for m = 1:length(muvals)
for mm =1:length(lambdavals)
    Utheory(m,mm)= Uint(mean(zvals),mean(Lvals),lambdavals(m),muvals(mm));
    Ljoint10(m,mm) = logL(data10,zvals10,Lvals10,lambdavals(m),muvals(mm));
    Ljoint25(m,mm) = logL(data25,zvals25,Lvals25,lambdavals(m),muvals(mm));
    Ljoint50(m,mm) = logL(data50,zvals50,Lvals50,lambdavals(m),muvals(mm));
    Ljoint100(m,mm) = logL(data100,zvals100,Lvals100,lambdavals(m),muvals(mm));
    Ljoint250(m,mm) = logL(data250,zvals250,Lvals250,lambdavals(m),muvals(mm));
    Ljoint500(m,mm) = logL(data500,zvals500,Lvals500,lambdavals(m),muvals(mm));
    Ljoint1000(m,mm) = logL(data1000,zvals1000,Lvals1000,lambdavals(m),muvals(mm));
end 
end 


Uintval = Uint(0.5,1,kbirth_true/D_true,kdeath_true/D_true);
muvalsline = linspace(0,20,100);
z=0.5;
L=1;
lambdavalsline = Uintval.*muvalsline.*(1+(-1).*cosh(z.*muvalsline.^(1/2))+sinh(z.*muvalsline.^(1/2)).*tanh((1/2) ...
.*L.*muvalsline.^(1/2))).^(-1);

figure('position',[60 897 890 111]);
subplot(1,7,1);
imagesc(lambdavals,muvals, Ljoint10'-max(Ljoint10,[],'all'))
hold on;
plot(lambdavalsline,muvalsline);
scatter(kbirth_true/D_true,kdeath_true/D_true)
set(gca,'YDir','normal')
axis square;
clim([-100, 0])
title('N=10','fontweight','normal')
set(gca,'LineWidth',1.25)


subplot(1,7,2);
imagesc(lambdavals,muvals, Ljoint25'-max(Ljoint25,[],'all'))
hold on;
plot(lambdavalsline,muvalsline);
scatter(kbirth_true/D_true,kdeath_true/D_true)
set(gca,'YDir','normal')
axis square;
clim([-100, 0])
title('N=25','fontweight','normal')
set(gca,'LineWidth',1.25)


subplot(1,7,3);
imagesc(lambdavals,muvals, Ljoint50'-max(Ljoint50,[],'all'))
hold on;
plot(lambdavalsline,muvalsline);
scatter(kbirth_true/D_true,kdeath_true/D_true)
set(gca,'YDir','normal')
axis square;
clim([-100, 0])
set(gca,'LineWidth',1.25)

title('N=50','fontweight','normal')

subplot(1,7,4);
imagesc(lambdavals,muvals, Ljoint100'-max(Ljoint100,[],'all'))
hold on;
plot(lambdavalsline,muvalsline);
scatter(kbirth_true/D_true,kdeath_true/D_true)
set(gca,'YDir','normal')
axis square;
clim([-100, 0])
title('N=100','fontweight','normal')
set(gca,'LineWidth',1.25)


subplot(1,7,5);
imagesc(lambdavals,muvals, Ljoint250'-max(Ljoint250,[],'all'))
hold on;
plot(lambdavalsline,muvalsline);
scatter(kbirth_true/D_true,kdeath_true/D_true)
set(gca,'YDir','normal')
axis square;
clim([-100, 0])
title('N=250','fontweight','normal')

set(gca,'LineWidth',1.25)

subplot(1,7,6);
imagesc(lambdavals,muvals, Ljoint500'-max(Ljoint500,[],'all'))
hold on;
plot(lambdavalsline,muvalsline);
scatter(kbirth_true/D_true,kdeath_true/D_true)
set(gca,'YDir','normal')
axis square;
clim([-100, 0])
title('N=500','fontweight','normal')
set(gca,'LineWidth',1.25)


subplot(1,7,7);
imagesc(lambdavals,muvals, Ljoint1000'-max(Ljoint1000,[],'all'))
hold on;
plot(lambdavalsline,muvalsline);
scatter(kbirth_true/D_true,kdeath_true/D_true)
set(gca,'YDir','normal')
axis square;
clim([-100, 0])
title('N=1000','fontweight','normal')

colormap  plasma;
set(gca,'LineWidth',1.25)

%%
% test rubin gelman Rhat

figure;

 addpath('mcmcdiag');

 chains(:,:,1)=outsamps10;
 chains(:,:,2)=outsamps10_2;
  rhat = cpsrf(chains,100);
 plot(abs(rhat(:,1)-1),'LineWidth',2,'color',[27,158,119]/256);

 set(gca,'YScale','log')
 set(gca,'XScale','log')
 axis square;
 box off;
 xlabel('MCMC step')
ylabel('Rubin-Gelman $|\hat{R}-1|$','interpreter','latex')
pbaspect([4 2 1])
axis tight
