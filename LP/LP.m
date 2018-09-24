%% Linear programming assignment
% Course: SC42055 Optimization in Systems and Control
% Jacob Lont, 4424409 and Casper van Engelenburg, 4237080

%% Always first pull from Github before making any changes

% E1, E2 and E3 are parameters changing from 0 to 18 for each group according to the sum of
% the last three numbers of the student IDs:
% E1 = Da1 + Db1, E2 = Da2 + Db2, E3 = Da3 + Db3,
% where Da;3 is the right-most digit of one student and Db;3 is the right-most digit of the other
% student.

EE1 = 4+0; EE2 = 0+8; EE3 = 9+0;

num_employees = 100 + EE2;
hours_per_month = 160;
month_salary = 3000 + 50*EE3;
salary_costs = num_employees*month_salary; % Company costs on employees

storage_space = (15+EE3)*10^3; % m^2 storage stock space available
space_R = 10; % m^2 taken by a model R car
space_W = 12; % m^2 taken by a model W car

% The model R and model W are respectively being sold for 55000e and 75000e.
price_R = 55*10^3; % Euro
price_W = 75*10^3; % Euro

% the cost of manufacturing the cars without considering the worker salaries is 30000e for the
% model R and 45000e for the model W.
costs_R = 30*10^3; % Euro without worker salaries
costs_W = 45*10^3; % Euro without worker salaries

%% Formulating LP problem using linprog.m

% CONCISE FORMULATION:
% -----------------------------------------
%  min  -fx(x1,x2)
% x1,x2
%
% st    x1,x2 integer (including 0) (or: x1,x2 in R+)
%       Ax <= b : inequality constraint
% -----------------------------------------
%
% EQUATIONS EXPLICITLY: 
% fx = P1*x1 + P1*x2 - Cm : objective value
% A = [c1 c2; H1 H2; A1 A2]
% x = [x1 x2]'
% b = [ctot Htot Atot]'
%
% PARAMETERS
% x1(x2): amount of model R(W) sold per month
% P1(P2): price of model R(W) (excluding manufacturing costs)
% Cm: manufacturing costs (solid number in this question)
% c1(c2): battery cells needed to manufacture 1 model R (W)
% H1(H2): hours needed to manufacture 1 model R (W)
% A1(A2): storage needed to store 1 model R (W)
% ctot: maximum amount of battery cell production p/month
% Htot: maximum amount of available hours p/month
% Atot: maximum amount of available storage space

%parameter set

%price
P1 = price_R - costs_R;
P2 = price_W - costs_W;
P = -[P1 P2]';
Cm = salary_costs;

%battery cells
c1 = 4E3;
c2 = 6E3;
ctot = (5+EE1)*1E6;

%manufacturing hours
H1 = 10;
H2 = 15;
Htot = hours_per_month*num_employees;

%storage
A1 = space_R;
A2 = space_W;
Atot = storage_space;

%objective, inequality constraints into concise matrix notation
A = [c1 c2; H1 H2; A1 A2];
b = [ctot Htot Atot]';

%bounds on x and initial condition and equality constraints
lb = [0 0]';
ub = [ ]';
x0 = [0 0]';
Ae = [];
be = [];

%linprog
options = optimoptions('linprog','Algorithm','dual-simplex');
[x fval flag] = linprog(P, A, b, Ae, be, lb, ub, x0, options);

profit1 = fval - Cm; 











