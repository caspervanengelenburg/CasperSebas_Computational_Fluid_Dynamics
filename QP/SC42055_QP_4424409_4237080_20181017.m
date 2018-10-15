%% Quadratic programming assignment
% Course: SC42055 Optimization in Systems and Control
% Jacob Lont, 4424409 and Casper van Engelenburg, 4237080
% Deadline: 17-10-2018


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ALWAYS PULL FIRST FROM GITHUB BEFORE MAKING ANY CHANGES %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Determine where this m-file's folder is and add that folder plus all subfolders to the path.
% Mainly important for importing the datafiles.
addpath(genpath(fileparts(which(mfilename))));

%% Initial parameter set

% E1, E2 and E3 are parameters changing from 0 to 18 for each group according to the sum of
% the last three numbers of the student IDs:
% E1 = Da1 + Db1, E2 = Da2 + Db2, E3 = Da3 + Db3,
% where Da;3 is the right-most digit of one student and Db;3 is the right-most digit of the other
% student.

EE1 = 4+0; EE2 = 0+8; EE3 = 9+0;

num_meas = 101+EE1; %total number of measurements 
dt = 3600; %time step


%% 1.



%% 2.

%loading data
% read_meas = readtable('measurements.csv');
% read_demand = readtable('demand.csv');
% read_prices = readtable('inputprices.csv');
read_meas = readtable(strcat(pwd, '\Data\measurements.csv'));
read_demand = readtable(strcat(pwd, '\Data\heatDemand.csv'));
read_prices = readtable(strcat(pwd, '\Data\inputPrices.csv'));



% optimoptions(’quadprog’,’Algorithm’,’active-set’)
% optimoptions(’quadprog’,’Algorithm’,’interior-point-convex’)
% optimoptions(’quadprog’,’MaxIter’,100)

% x = quadprog(H,c,A,b,Aeq,Beq,lb,ub,x0,options)
 
%heat transfer data
Qin = table2array(read_meas(:,2)); %heat transfer tank in
Qout = table2array(read_meas(:,3)); %heat transfer tank out
Q = [Qin Qout]'; %array of both heat transfer 
B = [1 -1]; %B matrix

%temperature data
T = table2array(read_meas(:,4)); %tank temperature T_k
Tamb = table2array(read_meas(:,5)); %ambient temperature T_amb
Tk1 = table2array(read_meas(2:end,4)); %tank temperature time step ahead T_{k+1}

%temperature differences
DTamb = T - Tamb; %difference between tank temperature and ambient temperature T_k - T_amb
DTk1 = Tk1 - T(1:end-1); %difference between tank temperatures per time step T_{k+1) - T_k

%determine H and c matrices 
H = zeros(2,2);
Hd = zeros(2,2);
c = zeros(2,1);
cd = zeros(2,1);


for idx = 1: num_meas
    
    Hd = dt^2* [ (DTamb(idx))^2            B*Q(:,idx)*DTamb(idx);
                 B*Q(:,idx)*DTamb(idx)          (B*Q(:,idx))^2    ];
    H = H + Hd;
    
    cd = dt*[2*DTk1(idx)*DTamb(idx);
             -2*DTk1(idx)*B*Q(:,idx)];
    c = c + cd;
end

H = H*2;

%optimize the variables
A = [];
b = [];
Aeq = [];
Beq = [];
lb = [0 0];
ub = [inf inf];
a0 = [];
o = optimoptions('quadprog', 'Algorithm', 'interior-point-convex');

a = quadprog(H,c,A,b,Aeq,Beq,lb,ub,a0,o);

% a =
% 
%    1.0e-06 *
% 
%     0.1349
%     0.0037


%% 3
N = 360; % 15 days horizon expressed in hours
Qin_max = (100 + EE2)*1E3;    % W (converted from kW using E3)
T1 = 330 + EE3;         % K(elvin)
Tamb = 275 + EE1;       % K(elvin)
Tmin = 315;             % K(elvin)
a = [1.96E-7; 3.8E-9];  % Given parameters


% (a) Is the problem quadratic? Justify.
% (b) Transform the prices to the corrects units so that they can be used in the optimi-
% zation problem. Which factor did you used to multiply them by?

% where lambda-k-in is the price of buying one unit of input heat at time step
% lambda_k_in is imported as [euro/MWh]
% Qout is imported as [W]

price = read_prices.Price * 1E-6; % Euro/Wh


% (c) What is the optimal cost of buying the input energy?
A_term = 1-a(1)*dt; %The term describing A in equation (3) of the QP assignment
T1_term = [A_term*T1 zeros(1,N-1)]';
c_k = a(1)*dt*Tamb*ones(N,1);

P = [dt*price(1:N)' zeros(N,1)']; % N input prices, N times 0
Ae = [a(2)*dt*eye(N) full(gallery('tridiag',N,A_term,-1,0))]; %part1: 1-diagonal, part2: -1 on diagonal, -A_term underthe diagonal
be = a(2)*dt*read_demand.Heat_demand(1:N) - T1_term - c_k;



% Define min and max values
lb = [zeros(1,N) Tmin*ones(1,N)]';  
ub = [Qin_max*ones(1,N) Inf*ones(1,N)]';



options = optimoptions('linprog','Algorithm','dual-simplex');
% options = optimoptions('linprog','Algorithm','interior-point');
[x, fval, flag] = linprog(P, [], [], Ae, be, lb, ub, [], options);


%plot temperature as function of time
figure(1)
plot(x(N + 1:2*N))

%final cost
cost_f = fval


%% 4

% (a) Modify (5) to include the temperature constraint and the extra cost. The resulting
% problem should be a quadratic problem.
% (b) Solve it in Matlab using quadprog. What is the total optimal cost? How much of
% that is destined to pay the terminal cost?

% new variables
pq = (1+EE2)/10; %price for quadratic difference
Tmax = 368; %maximum temperature
Tref = 323; %reference temperature

% quadratic problem formulation
Hq = zeros(2*N,2*N); %H matrix in quadratic formulation
Hq(end,end) = 2*pq;

cq = [ dt*read_prices.Price(1:N)*1E-6; zeros(N-1,1); - 2*pq*Tref];

ubq = [Qin_max*ones(1,N) Tmax*ones(1,N)]';

% solve quadratic problem

o = optimoptions('quadprog', 'Algorithm', 'interior-point-convex');
[xq,fvalq,exitflagq,outputq] = quadprog(Hq,cq,[],[],Ae,be,lb,ubq,[],o);

%plot temperature as function of time
figure(2)
plot(xq(N + 1:2*N));

%final costs, terminal cost
cost_fq = fvalq + Tref^2*pq; %total minimum costs
cost_t = (xq(end) - Tref)^2*pq %terminal cost (TN+1 - Tref)^2*pq
costs = [cost_f; cost_fq] %cost linear vs quadratic

cost_rel = cost_t/cost_fq %relative costs terminal versus total 









%% 5. Extra

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EXTRA: PLAYING AROUND WITH THE VARIABLE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%5.1 higher terminal cost
%5.2 quadratic costs for every time step

%% 5.1 higher terminal costs

fac = 15:2:40;
costs_e = zeros(length(fac),1);
pqe = zeros(length(fac),1);

for idx = 1:length(fac)
    
pqe(idx) = pq*fac(idx);

% quadratic problem formulation
Hqe = zeros(2*N,2*N); %H matrix in quadratic formulation
Hqe(end,end) = pqe(idx);

cqe = [ dt*read_prices.Price(1:N)*1E-6; zeros(N-1,1); - 2*pqe(idx)*Tref];

[xqe,fvalqe,exitflagqe,outputqe] = quadprog(Hqe,cqe,[],[],Ae,be,lb,ubq,[],o);

costs_e(idx) = fvalqe + Tref^2*pqe(idx);

%plot temperature versus time
figure(3); 
    plot(xqe(N+1:2*N))
    hold on
end

figure(4)
plot(costs,pqe)

%% 5.2 quadratic cost for every time step

% quadratic problem formulation
pq = 1000;
Hqt = pq*eye(2*N); %H matrix in quadratic formulation

cqt = [dt*read_prices.Price(1:N)*1E-6; -2*pq*Tref*ones(N,1)];

[xqt,fvalqt,exitflagqt,outputqt] = quadprog(Hqt,cqt,[],[],Ae,be,lb,ubq,[],o);

fvalqt
costs_t = fvalqt + N*pq*Tref^2

%plot temperature versus time
figure(5); 
    plot(xqt(N+1:2*N))

    



