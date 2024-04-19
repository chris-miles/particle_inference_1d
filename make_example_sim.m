clear all
close all

rng(33333) % 

dt=1e-2;

kon = 5;
koff = 0.75;
Tfinal= 10;
R=3;
D=.5;
z=R/3;

nT = round(Tfinal/dt);

nparticles_pre = poissrnd(Tfinal*kon); % total number of particles produced in [0,Tfinal]
birth_times = rand(1,nparticles_pre)*Tfinal; % generate random birth times

survived = ones(1,nparticles_pre); % vector to store whether they will survive at the end
life_times = exprnd(1/koff,1,nparticles_pre); % generate lifetimes

death_times = birth_times + life_times; % death is birth + lifetime

death_times = min(death_times,Tfinal);


positions = zeros(nparticles_pre,nT);


for n=1:nparticles_pre
    moves_n = [0;sqrt(2*D*dt)*randn(nT-1,1)];
    position_n = cumsum(moves_n,1)+z;
    pos_norm = min([position_n, R-position_n]);
    sb=sqrt(2*D*dt);
    % escape_probs = exp(-2*(R-pos_norm(1:end-1)).*(R-pos_norm(2:end))/(sb^2));
    %andrewsbray_escape_fix = find(randvals<escape_probs,1);
    %randvals = rand(size(escape_probs));
    exits = find(position_n>R|position_n<0,1);
    if ~isempty(exits)
        kill_time = exits(1);
        death_times(n) = min(birth_times(n)+dt*kill_time,death_times(n));
    end
    positions(n,:) = position_n;
end


figure;
hold on;

for n=1:nparticles_pre
    start_time = birth_times(n);
    end_time = death_times(n);
    if end_time == Tfinal
        survived = 1;
    else
        survived=0;
    end

    times = start_time:dt:end_time;
    if survived==1
        plot(D*(times-Tfinal)/R^2, positions(n,1:length(times))/R,'color',[1,0,0]);
        scatter(D*(start_time-Tfinal)/R^2, z/R,'filled','MarkerFaceColor','blue');
        scatter(D*(end_time-Tfinal)/R^2, positions(n,length(times))/R,'filled','MarkerFaceColor','red');
    else
        if min(abs(positions(n,length(times))/R), abs(positions(n,length(times))/R - 1))<.1
            bdry_exit = 1;
        else
            bdry_exit=0;
        end
        if bdry_exit==1

            plot(D*(times-Tfinal)/R^2, positions(n,1:length(times))/R,'color',[0.9,0.9,0.9]);
            scatter(D*(start_time-Tfinal)/R^2, z/R,'filled','MarkerFaceColor','magenta');
            scatter(D*(end_time-Tfinal)/R^2, positions(n,length(times))/R,'filled','cyan');
        else
            plot(D*(times-Tfinal)/R^2, positions(n,1:length(times))/R,'color',[0.7,0.7,0.7]);
            scatter(D*(start_time-Tfinal)/R^2, z/R,'filled','MarkerFaceColor','yellow');
            scatter(D*(end_time-Tfinal)/R^2, positions(n,length(times))/R,'filled','black');

        end

    end

end
yline(0);
yline(R/R)
yline(z/R);

xlim([-D*Tfinal/(2*R^2),0])
ylim([-0.1*R/R,1.1*R/R])
pbaspect([1 .5 1])
box on;
grid on;