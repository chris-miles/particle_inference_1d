function [outsamps,acceptProb] = mhGibbs(data,zvals,Lvals,alphaLambda,betaLambda, muMax, sigma, steps,init, verbose)


lambdaVals = zeros(1,steps);
muVals = zeros(1,steps);

if isempty(init)
    lambdaVals(1) = gamrnd(alphaLambda,betaLambda^(-1),1);
    muVals(1) = rand*muMax;
else
    lambdaVals(1) = init(1);
    muVals(1) = init(2);
end

Vint = @(z,L,mu) Uint(z,L,1,mu);

M = length(data);
nVals = zeros(1,M);
for m = 1:M
    nVals(m) = length(data{m});
end
nSum = sum(nVals);

acceptCount = 0;
for s = 1:steps-1

    % load previous vals
    lambda_s = lambdaVals(s);
    mu_s = muVals(s);


    %%%%%%%%%%%%%%%%%%%
    % MH step for mu

    muProp = mu_s + randn*sigma;
    piProp = logL(data,zvals,Lvals,lambda_s,muProp);%-log(muProp);
    piOld = logL(data,zvals,Lvals,lambda_s,mu_s);%-log(mu_s);

    if muProp <0 || muProp>muMax
        p = 0;
    else
        p = min(exp(piProp-piOld),1);
    end

    accept = rand<p;

    if accept
        acceptCount = acceptCount+1;
        muNew = muProp;
    else
        muNew = mu_s;
    end


    %%%%%%%%%%%%%%%%%%%
    % Gibbs step for lambda
    alphaNew  = alphaLambda + nSum;


    VintTot = sum(Vint(zvals,Lvals,muNew));

    betaNew = betaLambda+ VintTot;
    lambdaNew = gamrnd(alphaNew,betaNew^(-1),1);

    

    lambdaVals(s+1) = lambdaNew;
    muVals(s+1) = muNew;

    if verbose
        disp([s,acceptCount/s])
    end

end

outsamps = [lambdaVals; muVals]';
acceptProb = acceptCount/steps;
end 