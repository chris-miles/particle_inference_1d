clear all
close all

rng(4)


D_true =1;%
kbirth_true = 500;
kdeath_true = 10;% decay rate

Tfinal = 100; % experiment length
Xmax=1; % size of domain [0,R];


Ncells = 250; % do 1000

Nexperiments=2500;


areas = zeros(1,Nexperiments);
areas0 = zeros(1,Nexperiments);


estimates = zeros(2,Nexperiments);
estimates0 = zeros(2,Nexperiments);

    Xsource_mean = 0.5; % position of source, 1d position
    Xs_noise = 0.75;

    Xmax_noise=.25;%;


    alphaL = Xmax^2/Xmax_noise^2;
    betaL  =  Xmax/Xmax_noise^2;
    
parfor sss=1:Nexperiments



    Nparticles = zeros(1,Ncells);
    data = cell(1, Ncells);
    zvals= zeros(1,Ncells);
    Lvals = zeros(1,Ncells);


    Nparticles0 = zeros(1,Ncells);
    data0 = cell(1, Ncells);
    zvals0= zeros(1,Ncells);
    Lvals0 = zeros(1,Ncells);

    %xagg =[];
    %  xagg0=[];
    % tic

    for n = 1:Ncells

        Xmax_n = gamrnd(alphaL,1/betaL,1,1);%Xmax+Xmax_noise*randn;

       % if Xmax_n<0
         %   Xmax_n=-Xmax_n;
        %end

        Xs_n = betarnd(1/Xs_noise,1/Xs_noise)*Xmax_n;%Xsource_mean + Xs_noise*randn;

       % if Xs_n<0
      %      Xs_n = -Xs_n;
   % elseif Xs_n>Xmax_n
    %        Xs_n = 2*Xmax_n-Xs_n;
     %   end

        Lvals0(n) = Xmax;
        zvals0(n) = Xsource_mean;

        Lvals(n) = Xmax_n;
        zvals(n) = Xs_n;

        data0{n} = montecarlosim_1d(Tfinal, Xmax, Xsource_mean, kbirth_true, kdeath_true, D_true);
        data{n} = montecarlosim_1d(Tfinal, Xmax_n, Xs_n, kbirth_true, kdeath_true, D_true);


        positions0 = data0{n};

        positions = data{n};
        % xagg = [xagg;positions];
        % xagg0=[xagg0;positions0];

        Nparticles0(n) = length(positions0);
        Nparticles(n) = length(positions);

    end


    %%%

    sigmaInit = 1;
    meanLambda = 250;%0.5*kbirth_true/D_true;
    varLambda = meanLambda^2;

    alphaLambda = meanLambda^2/varLambda;
    betaLambda  =  meanLambda/varLambda;
    muMax = 5*(kdeath_true/D_true);

    acceptThresh1=0.3;
    acceptThresh2=0.5;

    sWarmup=250;
    steps=2e4;


    [outsamps,acceptProb,postEvals] = mhGibbsWithWarmup(data,zvals,...,
        Lvals,alphaLambda,betaLambda, muMax, sigmaInit, steps,sWarmup,...
        acceptThresh1,acceptThresh2);

    [outsamps0,acceptProb0,postEvals0] = mhGibbsWithWarmup(data0,zvals0,...,
        Lvals0,alphaLambda,betaLambda, muMax, sigmaInit, steps,sWarmup,...
        acceptThresh1,acceptThresh2);


    alpha = .11;

    q = quantile(postEvals,[alpha, 1]);
    hi=outsamps(postEvals>q(1),:);
    [bdry,area]  = boundary(hi);

    q0 = quantile(postEvals0,[alpha, 1]);
    hi0=outsamps0(postEvals0>q0(1),:);
    [bdry0,area0]  = boundary(hi0);

 
        estimates(:,sss) = mean(outsamps(end/2:end,:));
    estimates0(:,sss) = mean(outsamps0(end/2:end,:));


    areas(sss)=area;
    areas0(sss)=area0;
    disp(sss);

end

save('postareas_both')
%figure;
%histogram(areas,'normalization','pdf')
%hold on;
%histogram(areas0,'normalization','pdf')

