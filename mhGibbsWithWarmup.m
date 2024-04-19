function [outsamps,acceptProb,postEvals] = mhGibbsWithWarmup(data,zvals,...,
    Lvals,alphaLambda,betaLambda, muMax, sigmaInit, steps,sWarmup,...
acceptThresh1,acceptThresh2)


verbose=0;

sigma = sigmaInit;


currentAccept = 0;


warmUps = 1;
currentInit = [];



while currentAccept<acceptThresh1||currentAccept>acceptThresh2

[outWarm,acceptProbWarm] = mhGibbs(data,zvals,Lvals,alphaLambda,betaLambda, muMax, sigma, sWarmup,currentInit,verbose);
currentAccept = acceptProbWarm;

if acceptProbWarm<acceptThresh1
    sigma= sigma/2;
else
    sigma=sigma*2;
end 

currentInit = outWarm(end,:);
warmUps=warmUps+1;

end 

%disp('num warmups')
%disp(warmUps)


%steps=5e4;

tic
[outsamps,acceptProb] = mhGibbs(data,zvals,Lvals,alphaLambda,betaLambda, muMax, sigma, steps,currentInit,verbose);
toc
%disp(acceptProb)
postEvals = zeros(1,steps);

lgampdf = @(x,alpha,beta) alpha.*log(beta) + (alpha-1).*log(x)-beta.*x-log(gamma(alpha));

for s=1:steps
postEvals(s) = logL(data,zvals,Lvals,outsamps(s,1),outsamps(s,2))...
    +lgampdf(outsamps(s,1),alphaLambda,betaLambda);
end 
