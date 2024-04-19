function out =  logL(data,zvals,Lvals,lambda,mu)

M = length(data);
outVals = zeros(1,M);


for m = 1:M
 

   
    z_m =zvals(m);
    L_m = Lvals(m);

    out_m = -Uint(z_m,L_m,lambda,mu);

    x_m = data{m};

    %N_m = length(x_m);

    out_m = out_m+ sum(log(u(x_m,z_m,L_m,lambda,mu)));

  %  for n=1:N_m
       %logpx = log(u(x_m(n),z_m,L_m,lambda,mu));
  %     out_m = out_m+log(u(x_m(n),z_m,L_m,lambda,mu));
  %  end

outVals(m) = out_m;

end

out = sum(outVals);