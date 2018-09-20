%% Linear programming assignment
% Course: SC42055 Optimization in Systems and Control
% Jacob Lont, 4424409 and Casper van Engelenburg, 4237080

%% Always first pull from Github before making any changes

% E1, E2 and E3 are parameters changing from 0 to 18 for each group according to the sum of
% the last three numbers of the student IDs:
% E1 = Da1 + Db1, E2 = Da2 + Db2, E3 = Da3 + Db3,
% where Da;3 is the right-most digit of one student and Db;3 is the right-most digit of the other
% student.

E1 = 4+0; E2 = 0+8; E3 = 9+0;

num_employees = 100 + E2;
hours_per_month = 160;
month_salary = 3000 + 50*E3;
salary_costs = num_employees*month_salary; % Company costs on employees

storage_space = (15+E3)*10^3; % m^2 storage stock space available
space_R = 10; % m^2 taken by a model R car
space_W = 12; % m^2 taken by a model W car

% The model R and modelWare respectively being sold for 55000e and 75000e.
price_R = 55*10^3; % Euro
price_W = 75*10^3; % Euro

% the cost of manufacturing the cars without considering the worker salaries is 30000e for the
% model R and 45000e for the model W.
costs_R = 30*10^3; % Euro without worker salaries
costs_W = 45*10^3; % Euro without worker salaries