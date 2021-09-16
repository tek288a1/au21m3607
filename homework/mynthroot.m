% script m-file: mynthroot.m
% Calculates all n-th roots of negative unity given n.

% n = input('n = ');
x = exp(pi*1i*(1:2:2*n-1)'/n);       % all roots computed in one statement
disp(' ')
fprintf('** n = %d\n', n)
for j = 1:n
    fprintf('%8.4f + %8.4fi\n', real(x(j)), imag(x(j)))
end
