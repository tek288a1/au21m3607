% script m-file: leapYear.m
% Determines whether a given year is a leap year.

% year = input('Enter year: ');
if mod(year,4)~=0 || ( mod(year,100)==0 && mod(year,400)~=0 )
    fprintf('Year %4d is a common year.\n', year)
else
    fprintf('Year %4d is a leap year.\n', year)
end

