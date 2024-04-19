clear all
close all


lambda=500;
Lfix=1;

addpath('cbrewer')
n_mus_to_plot = 20;


muvals = logspace(-1,2,51); 

mus_to_plot = round(linspace(1,length(muvals),n_mus_to_plot));



zvals=linspace(.0,1,301);
Lvals = linspace(0.01,10,length(zvals));


IoutsZ = zeros(length(muvals),length(zvals));
IoutsL = zeros(length(muvals),length(zvals));


sigmavals = linspace(.01,.99,31);
IoutavgZ = zeros(length(muvals),length(sigmavals));
IoutavgL = zeros(length(muvals),length(sigmavals));




for mm = 1:length(muvals)
     mu=muvals(mm);

    parfor zz = 1:length(zvals)
         mu=muvals(mm);

    z = zvals(zz);
    IoutsZ(mm,zz)=I(z,Lfix,lambda,mu);
    end


    parfor ll = 1:length(Lvals)
             mu=muvals(mm);

    LL = Lvals(ll);
    IoutsL(mm,ll)=I(0.5*LL,LL,lambda,mu);
    end


    %Ioutmid(mm)=I(0.5*Lfix,Lfix,lambda,mu);

    
    parfor ss=1:length(sigmavals)
    sigma_source = sigmavals(ss);
    betavals = betapdf(zvals,1/sigma_source,1/sigma_source);
    gammavals = gampdf(Lvals,1/(sigma_source^2), (sigma_source^2)/1);

    intZ=IoutsZ(mm,:).*betavals;
    intZ(isnan(intZ))=0;

    intL=  IoutsL(mm,:).*gammavals;
    intL(isnan(intL))=0;
    
    IoutavgZ(mm,ss) = trapz(zvals, intZ);
    IoutavgL(mm,ss) = trapz(Lvals, intL);
    end 
end 


%%%%%%%%%%%%%
figure

addpath('cbrewer')
cc=flipud(cbrewer('seq','YlGn',31,'linear'));
%cc=flipud(cbrewer('seq','YlGn',31,'linear'));
cc = cmasher('bubblegum',n_mus_to_plot)
%%
subplot(1,2,1);
s=surf(zvals, muvals,IoutsZ);
s.EdgeAlpha=0;
view(2);

colormap viridis;
axis tight;
%set(gca,'ColorScale','log')
set(gca,'YScale','log')
xlabel('source location, z')
ylabel('degradation rate \mu')
axis square;
cbar=colorbar;
set(gca,'ColorScale','log')
clim([1e-8,1e-5])
cbar.Label.String = 'information, I';



subplot(1,2,2)
for m=1:n_mus_to_plot
plot(zvals,IoutsZ(mus_to_plot(m),:),'color',cc(m,:),'LineWidth',2)
hold on;
end 
set(gca,'YScale','log')
axis square;
xlabel('source loc, z')
ylabel('information, I')

cbar=colorbar;
colormap(viridis)
clim([muvals(1),muvals(end)])
cbar.Label.String='degradation rate, \mu'

box off;


figure;
for m=1:length(muvals)
[~,maxind]= max(IoutsZ(m,:));
maxZ(m)= zvals(maxind);
end 
maxZfix = abs(maxZ - 0.5)+0.5;
maxZfix_neg = abs(abs( 0.5-maxZ)-0.5);
plot(muvals,maxZfix)
hold on
plot(muvals,maxZfix_neg)

%set(gca,'YScale','log')
axis square;
xlabel('degradation rate, z')
ylabel('information, arg max I')
ylim([0 1])
box off;



figure;


for m=1:n_mus_to_plot
plot(sigmavals,IoutavgZ(mus_to_plot(m),:),'color',cc(m,:),'LineWidth',2)

hold on;
end
set(gca,'YScale','log')
axis square;
xlabel('source loc noise, \sigma')
ylabel('information, I')
box off;

cbar=colorbar;
colormap(cmasher('bubblegum'))
clim([muvals(1),muvals(end)])
cbar.Label.String='degradation rate, \mu'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure

%%
subplot(1,2,1);
s=surf(Lvals, muvals,IoutsL);
s.EdgeAlpha=0;
view(2);
colormap viridis;
axis tight;
%set(gca,'ColorScale','log')
set(gca,'YScale','log')
xlabel('bdry size, L')
ylabel('degradation rate \mu')
axis square;
cbar=colorbar;
set(gca,'ColorScale','log')
clim([1e-10,1])
cbar.Label.String = 'information, I';



subplot(1,2,2)
for m=1:n_mus_to_plot
plot(Lvals,IoutsL(mus_to_plot(m),:),'color',cc(m,:),'LineWidth',2)
hold on;
end 
set(gca,'YScale','log')
axis square;
xlabel('bdry size, L')
ylabel('information, I')

cbar=colorbar;
colormap(viridis)
clim([muvals(1),muvals(end)])
cbar.Label.String='degradation rate, \mu'

box off;


figure;
for m=1:n_mus_to_plot
plot(sigmavals,IoutavgL(mus_to_plot(m),:),'color',cc(m,:),'LineWidth',2)
hold on;
end
set(gca,'YScale','log')
axis square;
xlabel('bdry size noise, \sigma')
ylabel('information, I')
box off;

cbar=colorbar;
colormap(viridis)
clim([muvals(1),muvals(end)])
cbar.Label.String='degradation rate, \mu'




for m=1:length(muvals)
[~,maxind]= max(IoutavgL(m,:));
maxsigmaL(m)= sigmavals(maxind);
[~,maxind]= max(IoutavgZ(m,:));
maxsigmaZ(m)= sigmavals(maxind);
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Iout = I(z,L,lambda,mu)
xvals = linspace(0,L,5e4);
umuvals = umu(xvals,z,L,lambda,mu);
%ulambdavals = ulambda(xvals,z,L,lambda,mu);
uvals = u(xvals,z,L,lambda,mu);


%integrand1 = (ulambdavals.*ulambdavals)./(uvals);
%integrand2 = (ulambdavals.*umuvals)./(uvals);
integrand3 = (umuvals.*umuvals)./(uvals);


%integrand1(uvals==0)=0;
%integrand2(uvals==0)=0;
integrand3(uvals==0)=0;

%i11 = trapz(xvals,integrand1);
%i12 = trapz(xvals,integrand2);
i22 = trapz(xvals, integrand3);

%i11_mathematica_old = lambda.^(-1).*mu.^(-2).*(1+(-1).*cosh((L+(-1).*z).*mu)+sinh((L+( ...
%-1).*z).*mu).*tanh((1/2).*L.*mu));

%i12_mathematica_old = (1/2).*mu.^(-3).*sech((1/2).*L.*mu).^2.*((-2)+(-2).*cosh(L.*mu)+ ...
%2.*cosh((L+(-1).*z).*mu)+2.*cosh(z.*mu)+z.*mu.*sinh((L+(-1).*z).* ...
%mu)+(L+(-1).*z).*mu.*sinh(z.*mu));

i11_mathematica=lambda.^(-1).*mu.^(-1).*(1+(-1).*cosh((1/2).*(L+(-2).*z).*mu.^( ...
1/2)).*sech((1/2).*L.*mu.^(1/2)));

i12_mathematica = (1/2).*(exp(1).^(L.*mu.^(1/2))+exp(1).^(2.*L.*mu.^(1/2))).^(-2).*( ...
(-2).*exp(1).^(2.*L.*mu.^(1/2))+(-4).*exp(1).^(3.*L.*mu.^(1/2))+( ...
-2).*exp(1).^(4.*L.*mu.^(1/2))+exp(1).^((3.*L+z).*mu.^(1/2)).*(2+( ...
L+(-1).*z).*mu.^(1/2))+exp(1).^((2.*L+z).*mu.^(1/2)).*(2+(-1).*z.* ...
mu.^(1/2))+exp(1).^((4.*L+(-1).*z).*mu.^(1/2)).*(2+z.*mu.^(1/2))+ ...
exp(1).^((3.*L+(-1).*z).*mu.^(1/2)).*(2+((-1).*L+z).*mu.^(1/2))).* ...
mu.^(-2);

%disp([i11,i11_mathematica,i12,i12_mathematica])
Iout = det([i11_mathematica, i12_mathematica;i12_mathematica, i22]);
end 


