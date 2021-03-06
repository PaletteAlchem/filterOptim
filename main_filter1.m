%%-----------------------------------------------------------------------
% This is the source code for "Improved Flower Pollination Algorithm for 
% Synthesis of Cross-coupled Filters";
% Copyright: only for reaserch studies;
% Author: Pei-Wen Shu, Qing-Xin Chu;
% Date:30/01/2021;
% Email: eespw@mail.scut.edu.cn;
% Anyone who wish to use/cite it could contact with us.

%%
clear all��
clc

% This is for filter1;
%% Parameters of the iFPA method; 
n = 15; % Population size;
max_iteration = 1000; % Maximum iterations steps;
p = 0.8; % Phase transition probability;

%% Search space for each design variables;
Lb = [0, 0, 0, -1.2 ,0]; % Lower bound;
Ub = [1.2, 1.2, 1.2, 1.2, 1.2]; % Upper bound;
dim = length(Lb); % The number of design variables;

% Initial guess for coupling matrix (without cross coupling);
initial_sol = [0.8, 0.7, 0.8, 0.2, 1.05];

% Optimization begins;
[fitness_history_ifpa, best_solution_ifpa, best_fitness_ifpa] = iFPA(@objFunc1, n, Lb, Ub, max_iteration, p, initial_sol);

%% Post-processing;
% Convergence curve;
subplot(1,2,1);
% we use logrithm scale (base 10);
plot(log10(fitness_history_ifpa),'LineWidth',2);
title('Convergence curve');
xlabel('Iteration');
ylabel('Objective function value (Log. scale)');
grid();
xlim([0, max_iteration]);
legend('iFPA');

% Response from optimized coupling matrix;
[w, S21_iFPA, S11_iFPA] = PlotResponse_filter1(best_solution_ifpa);
subplot(1,2,2);
plot(w, 20*log10(abs((S21_iFPA))), 'LineWidth',2);
hold on;
plot(w, 20*log10(abs((S11_iFPA))), 'LineWidth',2);
title('Synthesis Responses');
xlabel('Normalized Frequency');
ylabel('S-Parameters');
legend('dB(S_{21})', 'dB(S_{11})');
grid();