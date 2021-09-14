% script m-file: coordConvert.m
% converts cartesian coordinates to spherical coordinates.

% x = input('Enter x: ');
% y = input('Enter y: ');
% z = input('Enter z: ');
fprintf('Given input: (x,y,z) = (%g, %g, %g)\n', x, y, z)
rho = sqrt(x^2 + y^2 + z^2);
phi = acos(z/rho);
if x == 0 && y == 0     % (x,y) = (0,0)
    theta = 0;
elseif x < 0            % 2nd and 3rd quad.
    theta = atan(y/x) + pi;
elseif y < 0            % 4th quad.
    theta = atan(y/x) + 2*pi;
else                    % 1st quad.
    theta = atan(y/x);
end
disp(' ')
fprintf('  rho:   %10.4f\n', rho)
fprintf('  phi:   %10.4f\n', phi)
fprintf('  theta: %10.4f\n', theta)
