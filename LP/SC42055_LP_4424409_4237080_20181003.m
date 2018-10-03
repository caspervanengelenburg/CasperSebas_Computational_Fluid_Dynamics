%% Linear programming assignment
% Course: SC42055 Optimization in Systems and Control
% Jacob Lont, 4424409 and Casper van Engelenburg, 4237080
% Deadline: 03-10-2018

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

%% 2 Formulating LP problem using linprog.m

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
ub = []';
x0 = [0 0]';
Ae = [];
be = [];

%linprog
options = optimoptions('linprog','Algorithm','dual-simplex');
[x, fval, flag] = linprog(P, A, b, Ae, be, lb, ub, x0, options);

profit2 = -fval - Cm  % Minus because it was a maximization problem initially


% RESULTS:

%  Optimal benefit:
% profit2 =
%    4.3573e+07
%  Number of cars:
%  x =
%         1728
%            0
% Limiting constraint: Htot = 17280, which is equal to 1728*10 hours



%% 3 A change in the market
ub = [1000 inf]';
[x, fval, flag] = linprog(P, A, b, Ae, be, lb, ub, x0, options);
profit3 = -fval - Cm  % Minus because it was a maximization problem initially

% RESULTS:
%  Optimal benefit:
% profit3 =
%     39187400
%  Number of cars:
%  x =
%     1000
%     485.3333 = 485



%% 5 further increase the profits of the company.
% for each additional worker that the company hires, the time required to 
% manufacture each car is reduced by 5 minutes. 
% this time reduction is limited to 6 hours, i.e. the maximum number of
% extra workers is limited to 72
ctot = (8+EE1)*1E6;
Atot = (22+EE3)*10^3;
ub = [1000 inf]';

for extra_employees = 0:72
    
    num_employees_total = num_employees + extra_employees;
    salary_costs = num_employees_total*month_salary; % Company costs on employees
    Cm = salary_costs;
    
    %manufacturing hours
    H1 = 10-(1/12)*extra_employees; 
    H2 = 15-(1/12)*extra_employees;
    Htot = hours_per_month*num_employees_total;

    %objective, inequality constraints into concise matrix notation
    A = [c1 c2; H1 H2; A1 A2];
    b = [ctot Htot Atot]';

    [x, fval, flag] = linprog(P, A, b, Ae, be, lb, ub, x0, options);
    
    % Store the data for analysis
    results(1,extra_employees+1) = num_employees_total;
    results(2,extra_employees+1) = -fval - Cm;
    results(3,extra_employees+1) = x(1);
    results(4,extra_employees+1) = x(2);
    
end

% Find the minimum of all these cases
[profit5, max_index] = max(results(2,:))  % Minus because it was a maximization problem initially
opt_num_workers = results(1,max_index);   % Optimal number of workers
x_max = [results(3,max_index),results(4,max_index)]; % Number of cars at the optimum

% RESULTS:
%  Optimal number of workers:
%       opt_num_workers = 144
%   Optimal benefit:
%       profit5 = 64503200
%  Number of cars:
%       x =
%          1000 model R
%          1333 model W



%% 6 Succes and a new model car!
% Because of the success, the 1000 limit on the model R is no longer in place 
% and they have signed contracts to produce at least 1250 model R and 
% 1000 model W per month.

c3 = 2E3;   % Battery cells for model V
A3 = 8;     % m^2 space taken by model V

extra_employees = opt_num_workers-num_employees;
H1 = 10-(1/12)*extra_employees; 
H2 = 15-(1/12)*extra_employees;
H3 = 8;     % 8 hours independently of the number of workers

price_V = 45E3; % Euro
costs_V = 45E3; % Euro

P1 = price_R - costs_R;
P2 = price_W - costs_W;
P3 = price_V - costs_V;
P = -[P1 P2 P3]';

Htot = hours_per_month*opt_num_workers;
ctot = (8+EE1)*1E6;
Atot = (22+EE3)*10^3;

%objective, inequality constraints into concise matrix notation
A = [c1 c2 c3; H1 H2 H3; A1 A2 A3];
b = [ctot Htot Atot]';

%bounds on x and initial condition and equality constraints
lb = [1250 1000 0]';    % implementing the contract constraints
ub = []';
x0 = [0 0]';
Ae = [];
be = [];

%linprog
options = optimoptions('linprog','Algorithm','dual-simplex');
[x, fval, flag] = linprog(P, A, b, Ae, be, lb, ub, x0, options);

profit6 = -fval - Cm  % Minus because it was a maximization problem initially

% RESULTS:
%   Optimal benefit:
%       profit5 = 66879000
%  Number of cars:
%       x =
%          1500 model R
%          1000 model W
%          0    model V


% is it economically beneficial to build the new model V?
%  No, it is not. It is better to stick to the model R and W

% If not, can Edison Automotive at least satisfy the new two contracts on 
% the model W and model R?
% Yes, Edison Automotive can satisfy the new two contracs. See the results.


