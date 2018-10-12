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