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



%% 1.



%% 2.

%loading data
% read_meas = readtable('measurements.csv');
% read_demand = readtable('demand.csv');
% read_prices = readtable('inputprices.csv');
read_meas = readtable(strcat(pwd, '\Data\measurements.csv'));
read_demand = readtable(strcat(pwd, '\Data\heatDemand.csv'));
read_prices = readtable(strcat(pwd, '\Data\inputPrices.csv'));



























