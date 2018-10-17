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

%plotsettings
LineWidth = 1; %set default linewidth for plots
MarkerSize = 10; %set default markersize for plots
FontSize = 18; %set default fontsize for plots


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
o2 = optimoptions('quadprog', 'Algorithm', 'trust-region-reflective');
a = quadprog(H,c,A,b,Aeq,Beq,lb,ub,a0,o);
a2 = quadprog(H,c,A,b,Aeq,Beq,lb,ub,a0,o2)

%%%%%%%%%%%
% RESULTS %
%%%%%%%%%%%
% a =
% 
%    1.0e-06 *
% 
%     0.1349 [1/s]
%     0.0037 [K/J]


%% 3
N = 360; % 15 days horizon expressed in hours
t = 1:N;
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

facJ = (10^6*dt)^(-1);
price = read_prices.Price * facJ; % conversion from Eu/MWh --> Euro/J


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



%costs
cost_f = fval;
cost_at_t = P(1:N)'.*x(1:N); %Euro

for idx = 1:N;
    
    tot_cost_at_t(idx) = sum(cost_at_t(1:idx));
    
end

%%%%%%%%%%%
% RESULTS %
%%%%%%%%%%%

%plot temperature as function of time
figure(1);clf;
    
    subplot(3,1,1)
    
    plot(t, x(N + 1:2*N), t, ones(N,1)*Tmin,'-.r', LineWidth, LineWidth)
    grid on;
    ylabel('$T  [ ^{\circ}C]$', 'Interpreter','latex');
    ylim([310 345]);
    xlim([0 N]);
    xlabel('Time $[hour]$', 'Interpreter','latex');
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', FontSize,'FontName','cmr12')
    title('Hourly energy trade without upper bound on T and without final Temperature optimization', 'Interpreter','latex')
    subplot(3,1,2)
    
    yyaxis left
    plot(t, x(1:N)/1000,t, ones(N,1)*Qin_max/1000,'-.r',LineWidth, LineWidth)
    grid on;
    ylabel('$\dot{Q_{in}}$ [kW]', 'Interpreter','latex'); 
    xlabel('Time $[hour]$', 'Interpreter','latex');
    xlim([0 N]);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', FontSize,'FontName','cmr12')
    
    yyaxis right
    plot(t, price(1:N)*dt, LineWidth, LineWidth);
    ylabel('$\lambda_{in} [\frac{Euro}{Wh}]$', 'Interpreter','latex');
    ylim([0 1.3*dt*max(price(1:N))]);
    
    subplot(3,1,3)
    
    yyaxis left
    plot(t, cost_at_t, LineWidth, LineWidth)
    grid on;
    ylabel('Price invested $[\frac{Euro}{hour}]$', 'Interpreter','latex'); 
    xlabel('Time $[hour]$', 'Interpreter','latex');
    xlim([0 N]);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', FontSize,'FontName','cmr12');
    
    yyaxis right
    plot(t, tot_cost_at_t, LineWidth, LineWidth);
    ylabel('Total price $[Euro]$', 'Interpreter','latex');
    ylim([0 1.3*max(tot_cost_at_t(1:N))]);
% fval = 127.5761 Euro


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

cq = [ dt*read_prices.Price(1:N)*facJ; zeros(N-1,1); - 2*pq*Tref];

ubq = [Qin_max*ones(1,N) Tmax*ones(1,N)]';


% solve quadratic problem
o = optimoptions('quadprog', 'Algorithm', 'interior-point-convex');
[xq,fvalq,exitflagq,outputq] = quadprog(Hq,cq,[],[],Ae,be,lb,ubq,[],o);


%final costs, terminal cost
cost_fq = fvalq + Tref^2*pq; %total minimum costs
cost_t = (xq(end) - Tref)^2*pq; %terminal cost (TN+1 - Tref)^2*pq
costs = [cost_f; cost_fq]; %cost linear vs quadratic

cost_rel = cost_t/cost_fq; %relative costs terminal versus total 

cost_at_tq = [cq(1:N-1).*xq(1:N-1); cq(N).*xq(N) + (xq(end) - Tref)^2*pq]; %Euro/h

for idx = 1:N;
    
    tot_cost_at_tq(idx) = sum(cost_at_tq(1:idx));
    
end


%%%%%%%%%%%
% RESULTS %
%%%%%%%%%%%

%plotting values
figure(2)
    subplot(3,1,1)
    
    plot(t, xq(N + 1:2*N), t, ones(N,1)*Tmin,'-.r', t, x(N + 1:2*N), '--b' , LineWidth, LineWidth)
    grid on;
    ylabel('$T  [ ^{\circ}C]$', 'Interpreter','latex');
    ylim([310 345]);
    xlim([0 N]);
    xlabel('Time $[hour]$', 'Interpreter','latex');
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', FontSize,'FontName','cmr12')
    title('Hourly energy trade with upper bound on T and final Temperature optimization', 'Interpreter','latex')
    
    subplot(3,1,2)
    
    yyaxis left
    plot(t, xq(1:N)/1000, t, ones(N,1)*Qin_max/1000,'-.r',t, x(1:N)/1000, '--b', LineWidth, LineWidth)
    grid on;
    ylabel('$\dot{Q_{in}}$ [kW]', 'Interpreter','latex'); 
    xlabel('Time $[hour]$', 'Interpreter','latex');
    xlim([0 N]);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', FontSize,'FontName','cmr12')
    
    yyaxis right
    plot(t, price(1:N)*dt, LineWidth, LineWidth);
    ylabel('$\lambda_{in} [\frac{Euro}{Wh}]$', 'Interpreter','latex');
    ylim([0 1.3*dt*max(price(1:N))]);
    
    subplot(3,1,3)
    
    yyaxis left
    plot(t, cost_at_tq, t, cost_at_t, '--b', LineWidth, LineWidth)
    grid on;
    ylabel('Price invested $[\frac{Euro}{hour}]$', 'Interpreter','latex'); 
    xlabel('Time $[hour]$', 'Interpreter','latex');
    xlim([0 N]);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', FontSize,'FontName','cmr12');
    
    yyaxis right
    plot(t, tot_cost_at_tq, t, tot_cost_at_t, '-.', LineWidth, LineWidth);
    ylabel('Total price $[Euro]$', 'Interpreter','latex');
    ylim([0 1.3*max(tot_cost_at_tq(1:N))]);
    

% total cost (cost_fq) = 142.78 [Euro]
% terminal cost (cost_t) = 1.2075 [Euro]
% relative cost of terminal (cost_rel) = 0.0085 [-] = 0.85 [%]