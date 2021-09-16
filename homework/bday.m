% script m-file: bday.m
%
% Calculates the number of birthday match in a single
% randomly chosen group of n people.

% n = input(' n = ');
nrs =  randi(365, 1, n);
nrs = sort(nrs);
d_nrs = diff(nrs);
id_nrs = find(d_nrs == 0);
nr_matches = length(id_nrs)
