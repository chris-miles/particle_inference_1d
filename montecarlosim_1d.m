function particles_at_end = montecarlosim_1d(Tfinal, R, Xsource, kon, koff, D)



dt=5e-6;

nparticles_pre = poissrnd(Tfinal*kon); % total number of particles produced in [0,Tfinal]

birth_times = rand(1,nparticles_pre)*Tfinal; % generate random birth times

survived = ones(1,nparticles_pre); % vector to store whether they will survive at the end 

life_times = exprnd(1/koff,1,nparticles_pre); % generate lifetimes 
death_times = birth_times + life_times; % death is birth + lifetime

survived(death_times<Tfinal) = 0; % if the death time is before T, they die

final_positions = zeros(nparticles_pre,1); % now we must move the remaining ones to see if they exit nucleus


for n=1:nparticles_pre
    if survived(n) % if they survived 
        steps_to_move = ceil((Tfinal-birth_times(n))/dt); % number of steps to generate
 
        moves_n = sqrt(2*D*dt)*randn(steps_to_move,1);
        position_n = cumsum(moves_n,1)+Xsource; % add up total moves to the end 


        pos_norm = min([position_n, R-position_n]);
        
        
        sb=sqrt(2*D*dt);
        escape_probs = exp(-2*(R-pos_norm(1:end-1)).*(R-pos_norm(2:end))/(sb^2));
        randvals = rand(size(escape_probs));
        andrewsbray_escape_fix = find(randvals<escape_probs,1);

       
        exits = find(position_n>R|position_n<0,1);  % check if the particle movement is ever NOT in the mask 
        
       % if isempty(exits)&&~isempty(andrewsbray_escape_fix)
       %     disp('fix helped');
      %  end 

      if ~isempty(andrewsbray_escape_fix)&&isempty(exits)
          disp('andrewsbrayfix mattered')
      end
        
        if ~isempty(exits)||~isempty(andrewsbray_escape_fix) % if it leaves, kill it 
            survived(n)=0;
        else % otherwise store its final position
            final_positions(n,:) = position_n(end);
        end
    end

end

% output remaining particles

particles_at_end = final_positions(find(survived));

