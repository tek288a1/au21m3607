% script m-file: isLeapYear.m
% 
% input: year
% output: l_LY (logical variable)
%          = true  (logical 1) if leap year
%          = false (logical 0) if not.

% year = input('Enter year: ');
if mod(year,4)~=0 || ( mod(year,100)==0 && mod(year,400)~=0 )
    l_LY = false;
else
    l_LY = true;
end

