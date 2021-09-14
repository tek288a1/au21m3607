% script m-file: oblateSpheroid.m
% This script calculates the surface area of an oblate spheroid. It also
% provides an approximate surface area.

% r1 = input('Enter equatorial radius (r1): ');
% r2 = input('Enter polar radius (r2 < r1): ');
gamma = acos(r2/r1);
A_exact = 2*pi*( r1^2 + ...
    r2^2/sin(gamma) * log( cos(gamma)/(1-sin(gamma)) ));
A_approx = 4*pi*( (r1+r2)/2 )^2;
disp('The surface area of the given spheroid:')
disp(['Exact:  ', num2str(A_exact)])
disp(['Approx: ', num2str(A_approx)])