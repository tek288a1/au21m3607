% script m-file: bdayMult.m
%
% Calculates the probability of having a birthday match in a randomly
% chosen group of n people by running multiple simulations.
% 
% inputs: n = number of people
%         nr_sims = number of simulations
%
% From pp.453-4 of Learning MATLAB by Overman

% n = input(' n = ');
% nr_sims = input(' number of simulations = ');
Nrs = randi(365, n, nr_sims);
Nrs = sort(Nrs);
D_Nrs = diff(Nrs);
Match = zeros(size(D_Nrs));    % preallocate: Match = 0 for no match
iD_Nrs = find(D_Nrs == 0 );    %              Match = 1 for a match
Match(iD_Nrs) = 1;
matches = sum(Match);
i_matches = find(matches > 0);
ratio = length(i_matches)/nr_sims;
ratio_ext = 1 - prod((365 - (0:n-1))./365);
disp(' ')
disp([' n = ', num2str(n), ', nr_sim = ', num2str(nr_sims)])
disp('------------------------------------------')
disp([' fraction of matches (experiment): ', num2str(ratio, 4)])
disp([' fraction of matches (analytical): ', num2str(ratio_ext, 4)])

