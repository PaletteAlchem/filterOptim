function K = objFunc1(x)
%%
% This is the objective function for filter1 

m12 = x(1);
m23 = x(2);
m34 = x(3);
m14 = x(4);
r11 = x(5);
rnn = x(5);

% Coupling matrix;
M = [0, m12, 0, m14;
     m12, 0, m23, 0;
     0, m23, 0, m34;
     m14, 0, m34, 0];
order = length(M);
R = zeros(order);
R(1,1) = r11;
R(order,order) = rnn;

% Computed from recursion;
ps = [-1.8, 1.8];
tz = [-0.9358, -0.4126, 0.4126, 0.9358]; 

% Return Loss in the passband;
RL = 20;
epsilon = (10^(RL/10)-1)^(-1/2);

% Different parts of objective function;
for i = 1:length(tz)    
    A = A_matrix(tz(i), M, R);
    [~,s11] = M2SP(A, R);
    K_part1_temp(i) = (abs(s11)^2);
end

for i = 1:length(ps)    
    A = A_matrix(ps(i), M, R);
    [s21,~] = M2SP(A, R);
    K_part2_temp(i) = (abs(s21)^2);
end

K_part1 = sum(K_part1_temp);

K_part2 = sum(K_part2_temp);

A = A_matrix(+1, M, R);
[~,s11] = M2SP(A, R);
K_part3 = (abs(s11) - epsilon/sqrt(1+epsilon^2))^2;

A = A_matrix(-1, M, R);
[~,s11] = M2SP(A, R);
K_part4 = (abs(s11) - epsilon/sqrt(1+epsilon^2))^2;

% Total cost function
K = K_part1 + K_part2 + K_part3 + K_part4;  

end

function A = A_matrix(w, M, R)
ONE = eye(length(M));
A = ONE * w - 1j * R + M;
end

function [s21,s11] = M2SP(A,R)
r11 = R(1,1);
rnn = R(end,end);
A_inv = pinv(A);
s21 = A_inv(length(A), 1) * -2j * sqrt(r11*rnn);
s11 = 1 + 2j * r11 * A_inv(1,1);
end
